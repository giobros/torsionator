# torsionfit/dihedral.py
import os
import glob
import numpy as np
import logging
from dataclasses import dataclass
import networkx as nx
from ase import neighborlist
from ase.io import read
from rdkit import Chem
import pandas as pd
from .io_utils import write_pdb, write_xyz_file
from .geometry import GeometryOptimizer, DihedralStepper, energy_eh
from .lj import lj_energy_between_groups, local_rg
from .constants import CONV_EH_TO_KCAL_MOL, VDW_RADII


class DihedralFinder:
    """Find rotatable dihedrals (four-atom tuples) from a PDB using RDKit."""
    @staticmethod
    def find_rotables(pdb_file: str):
        mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
        dihedral_indices = []
        uniq = {}

        for bond in mol.GetBonds():
            if bond.GetBondType() != Chem.BondType.SINGLE or bond.IsInRing():
                continue

            A = bond.GetBeginAtom()
            B = bond.GetEndAtom()
            AN = [a.GetIdx() for a in A.GetNeighbors() if a.GetIdx() != B.GetIdx()]
            BN = [a.GetIdx() for a in B.GetNeighbors() if a.GetIdx() != A.GetIdx()]

            if (A.GetDegree() == 1 and A.GetAtomicNum() == 6 and
                len([n.GetAtomicNum() for n in A.GetNeighbors()]) == 1):
                continue
            if (B.GetDegree() == 1 and B.GetAtomicNum() == 6 and
                len([n.GetAtomicNum() for n in B.GetNeighbors()]) == 1):
                continue
            if A.GetAtomicNum() == 1 or B.GetAtomicNum() == 1:
                continue
            if not AN or not BN:
                continue

            for i in AN:
                for l in BN:
                    if mol.GetAtomWithIdx(i).GetAtomicNum() > 1 and mol.GetAtomWithIdx(l).GetAtomicNum() > 1:
                        j = A.GetIdx()
                        k = B.GetIdx()
                        dihedral = (l, k, j, i)
                        pair = tuple(sorted([j, k]))
                        if pair not in uniq:
                            uniq[pair] = dihedral

        return list(uniq.values())


class GraphSplitter:
    """
    Split an ASE Atoms graph into two fragments (Group A / Group D) around a target dihedral (A-B-C-D).
    Robust to rings/fused systems. Call `split_by_dihedral(atoms, dihedral, method='auto')`.

    Methods:
      - method='auto'        : choose best strategy based on ring analysis
      - method='edge_removal': remove B–C edge, take components (fast; ok for acyclic)
      - method='node_removal': remove node B for A-side and node C for D-side (very robust in rings)
      - method='hybrid'      : try edge_removal; fallback to dijkstra if overlap/low coverage
      - method='dijkstra'    : assign by shortest-path distance; tie-break by paths avoiding B–C
    """

    # -------- Public API --------
    def split_by_dihedral(self, atoms, dihedral, method: str = "auto"):
        A, B, C, D = map(int, dihedral)
        G = self._build_connectivity_graph(atoms)

        if method == "auto":
            method = self._recommend_method(G, A, B, C, D)

        if method == "edge_removal":
            group_A, group_D = self._split_by_edge_removal(G, A, B, C, D)
        elif method == "node_removal":
            group_A, group_D = self._split_by_node_removal(G, A, B, C, D)
        elif method == "hybrid":
            group_A, group_D = self._split_hybrid(G, A, B, C, D)
        elif method == "dijkstra":
            group_A, group_D = self._split_by_dijkstra(G, A, B, C, D)
        else:
            raise ValueError(f"Unknown method: {method}")

        #  ensure disjoint and exclude dihedral atoms
        dihedral_set = {A, B, C, D}
        group_A = sorted(set(group_A) - dihedral_set)
        group_D = sorted(set(group_D) - dihedral_set)
        overlap = set(group_A) & set(group_D)
        if overlap:
            group_A = sorted(set(group_A) - overlap)
            group_D = sorted(set(group_D) - overlap)
        return group_A, group_D

    # -------- Internals --------
    def _build_connectivity_graph(self, atoms):
        """Build molecular connectivity from ASE neighbor list."""
        cutoffs = neighborlist.natural_cutoffs(atoms)
        nl = neighborlist.NeighborList(cutoffs, self_interaction=False, bothways=True)
        nl.update(atoms)

        G = nx.Graph()
        n = len(atoms)
        G.add_nodes_from(range(n))
        for i in range(n):
            neigh, _ = nl.get_neighbors(i)
            for j in neigh:
                if i < j:
                    G.add_edge(i, j)
        return G

    def _split_by_edge_removal(self, G, A, B, C, D):
        """Remove only the B–C edge, then take components for A and D."""
        Gs = G.copy()
        if Gs.has_edge(B, C):
            Gs.remove_edge(B, C)
        try:
            comp_A = set(nx.node_connected_component(Gs, A))
            comp_D = set(nx.node_connected_component(Gs, D))
        except nx.NetworkXError:
            return self._split_by_dijkstra(G, A, B, C, D)

        comp_A -= {A, B, C, D}
        comp_D -= {A, B, C, D}
        # Fallback if still overlapping
        if comp_A & comp_D:
            return self._split_by_dijkstra(G, A, B, C, D)
        return sorted(comp_A), sorted(comp_D)

    def _split_by_node_removal(self, G, A, B, C, D):
        """Remove node B to find A-side; remove node C to find D-side (robust in rings)."""
        GA = G.copy()
        if GA.has_node(B):
            GA.remove_node(B)
        comp_A = set(nx.node_connected_component(GA, A))

        GD = G.copy()
        if GD.has_node(C):
            GD.remove_node(C)
        comp_D = set(nx.node_connected_component(GD, D))

        comp_A -= {A, B, C, D}
        comp_D -= {A, B, C, D}
        overlap = comp_A & comp_D
        if overlap:
            comp_A -= overlap
            comp_D -= overlap
        return sorted(comp_A), sorted(comp_D)

    def _split_hybrid(self, G, A, B, C, D):
        """Try edge removal; validate coverage/overlap; fallback to dijkstra."""
        try:
            Gt = G.copy()
            if Gt.has_edge(B, C):
                Gt.remove_edge(B, C)
            comp_A = set(nx.node_connected_component(Gt, A))
            comp_D = set(nx.node_connected_component(Gt, D))
            overlap = comp_A & comp_D
            coverage = len(comp_A | comp_D) / len(G.nodes())
            if not overlap and coverage > 0.8:
                comp_A -= {A, B, C, D}
                comp_D -= {A, B, C, D}
                return sorted(comp_A), sorted(comp_D)
        except Exception:
            pass
        return self._split_by_dijkstra(G, A, B, C, D)

    def _split_by_dijkstra(self, G, A, B, C, D):
        """Assign atoms by shortest-path distance to A vs D; tie-break by paths avoiding B–C."""
        group_A, group_D = [], []
        dihedral_atoms = {A, B, C, D}

        # Precompute distances (unweighted) via BFS for efficiency
        dist_from_A = nx.single_source_shortest_path_length(G, A)
        dist_from_D = nx.single_source_shortest_path_length(G, D)

        for node in G.nodes():
            if node in dihedral_atoms:
                continue
            da = dist_from_A.get(node, None)
            dd = dist_from_D.get(node, None)
            if da is None or dd is None:
                continue
            if da < dd:
                group_A.append(node)
            elif dd < da:
                group_D.append(node)
            else:
                # Tie: prefer side with more shortest paths that avoid B–C
                paths_to_A = list(nx.all_shortest_paths(G, node, A))
                paths_to_D = list(nx.all_shortest_paths(G, node, D))
                a_avoid = sum(1 for p in paths_to_A if not self._path_uses_edge(p, B, C))
                d_avoid = sum(1 for p in paths_to_D if not self._path_uses_edge(p, B, C))
                if a_avoid > d_avoid:
                    group_A.append(node)
                elif d_avoid > a_avoid:
                    group_D.append(node)
                else:
                    group_A.append(node) 

        Aset, Dset = set(group_A), set(group_D)
        overlap = Aset & Dset
        if overlap:
            Aset -= overlap
            Dset -= overlap
        return sorted(Aset), sorted(Dset)

    def _path_uses_edge(self, path, u, v):
        """Return True if path uses edge (u, v) in either direction."""
        for i in range(len(path) - 1):
            if (path[i] == u and path[i+1] == v) or (path[i] == v and path[i+1] == u):
                return True
        return False

    # -------- Auto-selection helpers --------
    def _detect_ring_atoms(self, G):
        """Return set of atoms that belong to any cycle (ring)."""
        ring_atoms = set()
        try:
            for cyc in nx.cycle_basis(G):
                ring_atoms.update(cyc)
        except Exception:
            pass
        return ring_atoms

    def _recommend_method(self, G, A, B, C, D):
        """Choose a splitting method based on whether B/C lie in rings."""
        ring_atoms = self._detect_ring_atoms(G)
        B_ring = B in ring_atoms
        C_ring = C in ring_atoms
        if B_ring and C_ring:
            return "node_removal"
        elif B_ring or C_ring:
            return "hybrid"
        else:
            return "edge_removal"


@dataclass
class ClashScanResult:
    max_lj_kcal: float  # maximum LJ energy (kcal/mol) across the scan
    no_clash: bool      # True if no clashes were detected at any angle
    max_av_rg: float      # maximum average of local radii of gyration for the two groups


class ClashAnalyzer:
    """
    Scan a dihedral by rotating around it and:
      - compute the maximum Lennard-Jones energy between the split groups
      - detect steric clashes based on vdW radii (with a small buffer)
      - report a mean local radius of gyration
    Note: LJ energies here are already in kcal/mol given the parameter units.
    """
    def __init__(self, buffer=0.3):
        self.buffer = buffer

    def scan_check(self, atoms, dihedral, group_A, group_D, steps=37) -> ClashScanResult:
        init_angle = atoms.get_dihedral(*dihedral)

        max_lj = -1e9
        max_Rg = -1e9
        clash_any = False

        for i in range(steps):
            angle = (i / steps) * 360.0 + init_angle
            rot = atoms.copy()
            rot.set_dihedral(*dihedral, angle)

            # LJ energy in kcal/mol (parameters are kcal/mol)
            lj = lj_energy_between_groups(rot, group_A, group_D)
            rgA = local_rg(rot, group_A)
            rgD = local_rg(rot, group_D)
            avg_rg = 0.5 * (rgA + rgD)

            if lj > max_lj:
                max_lj = lj
            if avg_rg > max_Rg:
                max_Rg = avg_rg

            if self._has_clash(rot, group_A, group_D):
                clash_any = True

        return ClashScanResult(max_lj_kcal=max_lj, no_clash=not clash_any, max_av_rg=max_Rg)

    def _has_clash(self, atoms, group_A, group_D):
        for i in group_A:
            for j in group_D:
                d = np.linalg.norm(atoms.positions[i] - atoms.positions[j])
                r1 = VDW_RADII.get(atoms[i].symbol, 1.7)
                r2 = VDW_RADII.get(atoms[j].symbol, 1.7)
                threshold = r1 + r2 - self.buffer
                if d < threshold:
                    return True
        return False


class DihedralScanner:
    """
    Perform dihedral scans with constrained optimizations at a sequence of angles.
    Writes:
      - geometries.xyz (all frames)
      - scanned_pdbs/{angle}.pdb per snapshot
      - angles_vs_energies.txt (angles, energies in Eh)
      - angles_vs_energies_final.txt (angles sorted, relative Eh)
      - energies.dat (relative Eh column only)
    """
    def __init__(self, base_dir: str, optimizer: GeometryOptimizer):
        self.base_dir = base_dir
        self.optimizer = optimizer

    def scan(self, pdb_file, calc, dihedral_indices, method_label):
        atoms = read(pdb_file)
        dih_str = "_".join(map(str, dihedral_indices))
        out_dir = os.path.join(self.base_dir, "scanning", dih_str, method_label)
        xyz_file = os.path.join(out_dir, "geometries.xyz")
        txt_file = os.path.join(out_dir, "angles_vs_energies.txt")

        os.makedirs(out_dir, exist_ok=True)
        os.makedirs(os.path.join(out_dir, "scanned_pdbs"), exist_ok=True)

        start = atoms.get_dihedral(*dihedral_indices)
        angles = np.linspace(start, start + 360.0, 37)
        energies, true_angles = [], []
        a = atoms.copy()
        a.calc = calc
        stepper = DihedralStepper(self.optimizer.fmax, self.optimizer.steps)
        
        for angle in angles:

            a = stepper.optimize_with_dihedral(a, dihedral_indices, angle)

            e_eh = energy_eh(a)  # store in Hartree (Eh), consistent with original script
            tangle = a.get_dihedral(*dihedral_indices)

            energies.append(e_eh)
            true_angles.append(tangle)
            write_pdb(os.path.join(out_dir, "scanned_pdbs", f"{tangle:.0f}.pdb"), a)
            write_xyz_file(xyz_file, a)

        np.savetxt(txt_file, np.column_stack((true_angles, energies)))
        self._shift_files(out_dir)
    
    @staticmethod
    def _shift_files(output_folder):
        raw = os.path.join(output_folder, "angles_vs_energies.txt")
        final = os.path.join(output_folder, "angles_vs_energies_final.txt")
        edat = os.path.join(output_folder, "energies.dat")

        if not os.path.exists(raw):
            raise FileNotFoundError(raw)

        df = pd.read_csv(
            raw, delim_whitespace=True, header=None, names=["Angle", "Energy"]
        )

        df["Energy"] = df["Energy"] - df["Energy"].min()

        df["Energy"].to_csv(edat, sep=" ", index=False, header=False)
        df.sort_values("Angle").to_csv(final, sep=" ", index=False, header=False)


def plot_shifted_energy_profiles(base_dir: str, dihedral):
    """
    Plot relative energy profiles (kcal/mol) for obi and mace, if present.
    Reads .../angles_vs_energies_final.txt (energies in Eh) and converts to kcal/mol.
    """
    import matplotlib.pyplot as plt

    dih_str = "_".join(map(str, dihedral))
    base = os.path.join(base_dir, "scanning", dih_str)

    plt.figure(figsize=(10, 6))
    has_data = False

    for label in ("obi", "mace"):
        try:
            data = np.loadtxt(os.path.join(base, label, "angles_vs_energies_final.txt"))
            angles = data[:, 0]
            energies_kcal = data[:, 1] * CONV_EH_TO_KCAL_MOL  # Eh -> kcal/mol
            energies_kcal -= energies_kcal.min()
            plt.plot(angles, energies_kcal, 'o-', label=label)  
            has_data = True
        except Exception:
            pass

    if has_data:
        plt.xlabel("Dihedral Angle (degrees)")
        plt.ylabel("Relative Energy (kcal/mol)")
        plt.title(f"Dihedral Scan Energy Profile: {dih_str}")
        plt.legend()
        plt.xlim(0, 180)
        plt.grid(True)
        plt.tight_layout()
        out_png = os.path.join(base, f"{dih_str}.png")
        plt.savefig(out_png)
        plt.close()
        logging.info("Saved plot: %s", out_png)
