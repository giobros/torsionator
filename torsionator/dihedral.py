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
from collections import deque
from ase.data import covalent_radii

class GraphSplitter:
    """
    Split a molecule into two disjoint groups with respect to a dihedral (A,B,C,D),
    i.e., relative to the central bond B–C.
    """

    # ---------- Low-level helpers (ASE-only adjacency + BFS) ----------

    @staticmethod

    def split_two_groups(atoms, A, B, C, D, *, scale=1.10, pad=0.20, exclude_dihedral=False):
        """
        Build bonds by covalent radii and split by two-front BFS.
        - Edge (i,j) exists if distance(i,j) <= scale*(rc[i]+rc[j]) + pad  (rc = covalent radius).
        - Remove central edge B–C.
        - Seeds = first neighbors of B and of C (if empty -> A/D).
        - One queue with (node, label), label=0 for B-side, 1 for C-side.
        For each popped node, push all unvisited neighbors with the same label.
        - Any unvisited leftovers (rare) are assigned by graph-distance to B vs C.

        Returns (group_A, group_D) as sorted 0-based indices.
        """
        N = len(atoms)
        pos = atoms.get_positions()
        Z   = atoms.get_atomic_numbers()
        rc  = [covalent_radii[z] if z < len(covalent_radii) else 0.77 for z in Z]  # 0.77~C as fallback

        # --- 1) Ricostruisci i legami con criterio covalente ---
        adj = [set() for _ in range(N)]
        for i in range(N):
            ri = rc[i]
            pi = pos[i]
            for j in range(i+1, N):
                rj = rc[j]
                cutoff = scale * (ri + rj) + pad
                d = np.linalg.norm(pi - pos[j])
                # evita legami “zero” e collega se coerente con la somma dei rc
                if d > 0.25 and d <= cutoff:
                    adj[i].add(j); adj[j].add(i)

        # --- 2) Taglia B–C ---
        adj[B].discard(C); adj[C].discard(B)

        # --- 3) Seeds: primi vicini di B/C (fallback A/D) ---
        seedsB = list(adj[B]) or [A]
        seedsC = list(adj[C]) or [D]

        # --- 4) BFS a due fronti: “assegna se non ancora visitato” ---
        label = [-1] * N          # -1=unseen, 0=B-side (group_A), 1=C-side (group_D)
        dist  = [10**9] * N       # per eventuali tie/fallback

        q = deque()
        for s in seedsB:
            if label[s] == -1:
                label[s] = 0; dist[s] = 0; q.append(s)
        for s in seedsC:
            if label[s] == -1:
                label[s] = 1; dist[s] = 0; q.append(s)

        while q:
            u = q.popleft()
            lu = label[u]
            for v in adj[u]:
                if label[v] == -1:              # “vicino non ancora visitato”
                    label[v] = lu
                    dist[v] = dist[u] + 1
                    q.append(v)
                # se già visitato, NON lo rietichettiamo: “first arrival wins”

        # --- 5) Residui (grafo spezzato): assegna per distanza ai centri ---
        # BFS dai centri B e C sul grafo già senza B–C
        def bfs_single(src):
            d = {src: 0}; qq = deque([src])
            while qq:
                x = qq.popleft()
                for y in adj[x]:
                    if y not in d:
                        d[y] = d[x] + 1; qq.append(y)
            return d
        distB = bfs_single(B)
        distC = bfs_single(C)

        for v in range(N):
            if label[v] == -1:
                dB = distB.get(v, 10**9)
                dC = distC.get(v, 10**9)
                label[v] = 0 if dB <= dC else 1

        group_A = sorted(i for i in range(N) if label[i] == 0)
        group_D = sorted(i for i in range(N) if label[i] == 1)

        if exclude_dihedral:
            dihed = {A, B, C, D}
            group_A = [i for i in group_A if i not in dihed] or [A]
            group_D = [i for i in group_D if i not in dihed] or [D]

        # sanity
        assert not (set(group_A) & set(group_D))
        assert set(group_A) | set(group_D) == set(range(N))
        return group_A, group_D
    
    # ---------- Public API ----------

    def split_by_dihedral(self, atoms, dihedral, method: str = "neighbors",
                          *, exclude_dihedral: bool = False,
                             lock_cycle_to: str = "C",
                             mult: float = 1.0):
        """
        Entry point.

        Parameters
        ----------
        atoms : ase.Atoms
        dihedral : (int,int,int,int)
            (A,B,C,D) as 0-based ASE indices.
        method : str
            "neighbors"/"first_neighbors" or "cut_edge".
        exclude_dihedral : bool
            If True, remove {A,B,C,D} from groups (fallback [A]/[D] if a side ends empty).
        lock_cycle_to : str | None
            None, "B", or "C". Only used by the neighbors method.
        mult : float
            Multiplier for ASE natural cutoffs when building connectivity.

        Returns
        -------
        (group_A, group_D) : tuple[list[int], list[int]]
            Two disjoint lists whose union covers all atoms.
        """
        A, B, C, D = map(int, dihedral)

        if method in ("neighbors", "first_neighbors"):
            group_A, group_D = self.split_two_groups(
                atoms, A, B, C, D, scale=mult, exclude_dihedral=exclude_dihedral)
        elif method == "cut_edge":
            group_A, group_D = self._split_by_cut_edge_ase(
                atoms, A, B, C, D, mult=mult
            )
        else:
            raise ValueError(
                "Unknown method: {} (valid: neighbors/first_neighbors, cut_edge)".format(method)
            )

        if exclude_dihedral:
            dihed = {A, B, C, D}
            group_A = [i for i in group_A if i not in dihed] or [A]
            group_D = [i for i in group_D if i not in dihed] or [D]

        # final guarantees
        assert not (set(group_A) & set(group_D)), "Unexpected overlap between groups"
        assert set(group_A) | set(group_D) == set(range(len(atoms))), "Not all atoms were assigned"
        return group_A, group_D

    # ---------- Optional: keep a NetworkX graph builder if you need it elsewhere ----------

    def _build_connectivity_graph(self, atoms):
        """
        Build a networkx.Graph from ASE NeighborList (not used by the split itself).
        """
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


class DihedralFinder:
    """Find rotatable dihedrals (four-atom tuples) from a PDB using RDKit."""
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
