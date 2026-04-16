import logging
import os
from collections import deque
from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
from ase.data import covalent_radii
from ase.io import read
from rdkit import Chem

from .calculators import attach_calc
from .constants import CONV_EH_TO_KCAL_MOL, VDW_RADII
from .geometry import DihedralStepper, energy_eh
from .io_utils import write_pdb, write_xyz_file
from .lj import lj_energy_between_groups, local_rg


# =====================================================================
# Angle helpers
# =====================================================================

def _norm360(a: float) -> float:
    return ((a % 360.0) + 360.0) % 360.0


def _validate_step_size(step_size: int) -> int:
    try:
        s = int(step_size)
    except Exception as e:
        raise ValueError(f"step_size must be int-like, got {step_size!r}") from e
    if s <= 0:
        raise ValueError(f"step_size must be > 0, got {s}")
    if 360 % s != 0:
        raise ValueError(f"step_size must divide 360 exactly. Got {s}")
    return s


def _snap_to_step_nearest(a_deg: float, step: int) -> int:
    """Snap an angle to the nearest grid multiple of `step` in [0, 360)."""
    step = _validate_step_size(step)
    a = _norm360(a_deg)
    return int(np.floor(a / step + 0.5) * step) % 360


def _shift_scan_files(output_folder: str) -> None:
    """
    Post-process raw scan output:
      - angles_vs_energies.txt  → angles_vs_energies_final.txt (sorted by angle)
      - energies.dat            (energies shifted to min=0)
    """
    raw = os.path.join(output_folder, "angles_vs_energies.txt")
    final = os.path.join(output_folder, "angles_vs_energies_final.txt")
    edat = os.path.join(output_folder, "energies.dat")

    df = pd.read_csv(raw, sep=r"\s+", header=None, names=["Angle", "Energy"])
    df["Energy"] -= df["Energy"].min()
    df["Energy"].to_csv(edat, sep=" ", index=False, header=False)
    df.sort_values("Angle").to_csv(final, sep=" ", index=False, header=False)


# =====================================================================
# Graph splitter
# =====================================================================

class GraphSplitter:
    """
    Split a molecule into two disjoint groups relative to the central bond B–C
    of a dihedral (A, B, C, D).
    """

    @staticmethod
    def split_two_groups(
        atoms,
        A: int, B: int, C: int, D: int,
        *,
        scale: float = 1.10,
        pad: float = 0.20,
        exclude_dihedral: bool = False,
    ) -> Tuple[List[int], List[int]]:
        """Two-front BFS split after cutting the B–C bond."""
        N = len(atoms)
        pos = atoms.get_positions()
        Z = atoms.get_atomic_numbers()
        rc = [covalent_radii[z] if z < len(covalent_radii) else 0.77 for z in Z]

        # Build adjacency list from covalent radii
        adj = [set() for _ in range(N)]
        for i in range(N):
            for j in range(i + 1, N):
                cutoff = scale * (rc[i] + rc[j]) + pad
                d = np.linalg.norm(pos[i] - pos[j])
                if 0.25 < d <= cutoff:
                    adj[i].add(j)
                    adj[j].add(i)

        # Cut central bond
        adj[B].discard(C)
        adj[C].discard(B)

        # Two-front BFS
        seeds_b = list(adj[B]) or [A]
        seeds_c = list(adj[C]) or [D]

        label = [-1] * N
        dist = [10 ** 9] * N
        q = deque()

        for s in seeds_b:
            if label[s] == -1:
                label[s] = 0
                dist[s] = 0
                q.append(s)
        for s in seeds_c:
            if label[s] == -1:
                label[s] = 1
                dist[s] = 0
                q.append(s)

        while q:
            u = q.popleft()
            for v in adj[u]:
                if label[v] == -1:
                    label[v] = label[u]
                    dist[v] = dist[u] + 1
                    q.append(v)

        # Assign unlabelled atoms by graph distance to B vs C
        def _bfs_dist(src: int):
            dmap = {src: 0}
            qq = deque([src])
            while qq:
                x = qq.popleft()
                for y in adj[x]:
                    if y not in dmap:
                        dmap[y] = dmap[x] + 1
                        qq.append(y)
            return dmap

        dist_b = _bfs_dist(B)
        dist_c = _bfs_dist(C)

        for v in range(N):
            if label[v] == -1:
                label[v] = 0 if dist_b.get(v, 10 ** 9) <= dist_c.get(v, 10 ** 9) else 1

        group_A = sorted(i for i in range(N) if label[i] == 0)
        group_D = sorted(i for i in range(N) if label[i] == 1)

        if exclude_dihedral:
            dihed = {A, B, C, D}
            group_A = [i for i in group_A if i not in dihed] or [A]
            group_D = [i for i in group_D if i not in dihed] or [D]

        assert not (set(group_A) & set(group_D))
        assert set(group_A) | set(group_D) == set(range(N))
        return group_A, group_D

    def split_by_dihedral(
        self,
        atoms,
        dihedral: Tuple[int, int, int, int],
        method: str = "neighbors",
        *,
        exclude_dihedral: bool = False,
        mult: float = 1.0,
    ) -> Tuple[List[int], List[int]]:
        A, B, C, D = map(int, dihedral)
        if method not in ("neighbors", "first_neighbors"):
            raise ValueError(f"Unknown method: {method!r}")
        return self.split_two_groups(
            atoms, A, B, C, D, scale=mult, exclude_dihedral=exclude_dihedral
        )



# =====================================================================
# Dihedral finder (RDKit)
# =====================================================================

class DihedralFinder:
    """Find all rotatable dihedral four-atom tuples in a PDB using RDKit."""

    @staticmethod
    def find_rotables(pdb_file: str) -> List[Tuple[int, int, int, int]]:
        mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
        if mol is None:
            raise ValueError(f"RDKit could not read PDB: {pdb_file}")

        uniq = {}
        for bond in mol.GetBonds():
            if bond.GetBondType() != Chem.BondType.SINGLE or bond.IsInRing():
                continue

            A = bond.GetBeginAtom()
            B = bond.GetEndAtom()

            # Skip hydrogen-terminated or terminal methyl-like ends
            if A.GetAtomicNum() == 1 or B.GetAtomicNum() == 1:
                continue
            if A.GetDegree() == 1 and A.GetAtomicNum() == 6:
                continue
            if B.GetDegree() == 1 and B.GetAtomicNum() == 6:
                continue

            AN = [a.GetIdx() for a in A.GetNeighbors() if a.GetIdx() != B.GetIdx()]
            BN = [a.GetIdx() for a in B.GetNeighbors() if a.GetIdx() != A.GetIdx()]
            if not AN or not BN:
                continue

            for i in AN:
                for d in BN:
                    if (
                        mol.GetAtomWithIdx(i).GetAtomicNum() > 1
                        and mol.GetAtomWithIdx(d).GetAtomicNum() > 1
                    ):
                        j, k = A.GetIdx(), B.GetIdx()
                        pair = tuple(sorted([j, k]))
                        if pair not in uniq:
                            uniq[pair] = (d, k, j, i)

        return list(uniq.values())


# =====================================================================
# Dihedral validation
# =====================================================================

def _validate_dihedral_with_rdkit(
    pdb_file: str, a: int, b: int, c: int, d: int
) -> List[Tuple[int, int, int, int]]:
    mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
    if mol is None:
        raise ValueError(f"RDKit could not read PDB: {pdb_file}")

    n = mol.GetNumAtoms()
    for x in (a, b, c, d):
        if not (0 <= x < n):
            raise ValueError(f"Atom index {x} out of range [0, {n - 1}]")

    def bonded(i, j):
        return mol.GetBondBetweenAtoms(int(i), int(j)) is not None

    if not (bonded(a, b) and bonded(b, c) and bonded(c, d)):
        raise ValueError(
            f"Indices [{a},{b},{c},{d}] do not form a bonded dihedral "
            "(required bonds: a-b, b-c, c-d)."
        )
    return [(a, b, c, d)]


def resolve_dihedrals_arg(
    log, dihedral_arg: str, pdb_file: str
) -> Tuple[List[Tuple[int, int, int, int]], bool]:
    """
    Parse --dihedral argument.
    Returns (dihedrals, should_exit).
    """
    if dihedral_arg == "all":
        return DihedralFinder.find_rotables(pdb_file), False

    if dihedral_arg == "print":
        dihedrals = DihedralFinder.find_rotables(pdb_file)
        log.info("Dihedrals in the molecule:")
        for d in dihedrals:
            log.info("[%d,%d,%d,%d]", *d)
        return dihedrals, True

    if dihedral_arg.startswith("[") and dihedral_arg.endswith("]"):
        a, b, c, d = [int(x.strip()) for x in dihedral_arg.strip("[]").split(",")]
        return _validate_dihedral_with_rdkit(pdb_file, a, b, c, d), False

    raise ValueError("Invalid --dihedral argument. Use 'all', 'print', or '[a,b,c,d]'.")


# =====================================================================
# Clash analyzer
# =====================================================================

@dataclass
class ClashScanResult:
    max_lj_kcal: float
    no_clash: bool
    max_av_rg: float


class ClashAnalyzer:
    def __init__(self, buffer: float = 0.3):
        self.buffer = buffer

    def scan_check(
        self,
        atoms,
        dihedral: Tuple[int, int, int, int],
        group_A: List[int],
        group_D: List[int],
        steps: int = 36,
    ) -> ClashScanResult:
        init_angle = atoms.get_dihedral(*dihedral)
        max_lj = -np.inf
        max_rg = -np.inf
        clash_any = False

        for i in range(steps):
            angle = (i / steps) * 360.0 + init_angle
            rot = atoms.copy()
            rot.set_dihedral(*dihedral, angle)

            lj = lj_energy_between_groups(rot, group_A, group_D)
            avg_rg = 0.5 * (local_rg(rot, group_A) + local_rg(rot, group_D))

            max_lj = max(max_lj, lj)
            max_rg = max(max_rg, avg_rg)

            if self._has_clash(rot, group_A, group_D):
                clash_any = True

        return ClashScanResult(
            max_lj_kcal=float(max_lj),
            no_clash=not clash_any,
            max_av_rg=float(max_rg),
        )

    def _has_clash(self, atoms, group_A: List[int], group_D: List[int]) -> bool:
        for i in group_A:
            for j in group_D:
                d = np.linalg.norm(atoms.positions[i] - atoms.positions[j])
                r1 = VDW_RADII.get(atoms[i].symbol, 1.7)
                r2 = VDW_RADII.get(atoms[j].symbol, 1.7)
                if d < (r1 + r2 - self.buffer):
                    return True
        return False


# =====================================================================
# Dihedral scanners
# =====================================================================

class DihedralScanner:
    """
    Standard (non-MCS) dihedral scan.

    Snaps start angle to the nearest grid point, then sweeps a full 360°
    in step_size increments with constrained geometry optimization at each step.

    Writes per `out_dir`:
      - geometries.xyz
      - scanned_pdb/<angle>.pdb
      - angles_vs_energies.txt
      - angles_vs_energies_final.txt  (sorted by angle, energy shifted to min=0)
      - energies.dat                  (shifted energies only)
    """

    def __init__(self, base_dir: str, optimizer):
        self.base_dir = base_dir
        self.optimizer = optimizer

    def scan(
        self,
        pdb_file: str,
        calc,
        dihedral_indices: Tuple[int, int, int, int],
        method_label: str,
        step_size: int = 10,
    ) -> None:
        step_size = _validate_step_size(step_size)
        dih_str = "_".join(map(str, dihedral_indices))
        out_dir = os.path.join(self.base_dir, "scanning", dih_str, method_label)
        os.makedirs(out_dir, exist_ok=True)
        os.makedirs(os.path.join(out_dir, "scanned_pdb"), exist_ok=True)

        atoms = read(pdb_file)
        angle0 = _snap_to_step_nearest(float(atoms.get_dihedral(*dihedral_indices)), step_size)
        n_steps = 360 // step_size
        angles = [int((angle0 + i * step_size) % 360) for i in range(n_steps)]

        a = atoms.copy()
        attach_calc(a, calc)
        stepper = DihedralStepper(self.optimizer.fmax, self.optimizer.steps)

        a.set_dihedral(*dihedral_indices, angle0)
        a = stepper.optimize_with_dihedral(a, dihedral_indices, angle0)

        energies, saved_angles = [], []
        xyz_file = os.path.join(out_dir, "geometries.xyz")
        for angle in angles:
            a = stepper.optimize_with_dihedral(a, dihedral_indices, angle)
            energies.append(energy_eh(a))
            saved_angles.append(angle)
            write_pdb(os.path.join(out_dir, "scanned_pdb", f"{angle}.pdb"), a)
            write_xyz_file(xyz_file, a)

        np.savetxt(
            os.path.join(out_dir, "angles_vs_energies.txt"),
            np.column_stack((saved_angles, energies)),
        )
        _shift_scan_files(out_dir)


class DihedralScannerMCS:
    """
    MCS (multi-conformer) dihedral scan.

    Same as DihedralScanner but supports a scan direction ('cw' or 'ccw')
    and also writes scanned_xyz/<angle>.xyz alongside scanned_pdb/.
    """

    def __init__(self, base_dir: str, optimizer):
        self.base_dir = base_dir
        self.optimizer = optimizer

    def scan(
        self,
        pdb_file: str,
        calc,
        dihedral_indices: Tuple[int, int, int, int],
        method_label: str,
        direction: str = "cw",
        step_size: int = 10,
    ) -> None:
        step_size = _validate_step_size(step_size)
        if direction not in ("cw", "ccw"):
            raise ValueError(f"Invalid direction={direction!r}; expected 'cw' or 'ccw'")

        dih_str = "_".join(map(str, dihedral_indices))
        out_dir = os.path.join(self.base_dir, "scanning", dih_str, method_label)
        os.makedirs(out_dir, exist_ok=True)
        os.makedirs(os.path.join(out_dir, "scanned_pdb"), exist_ok=True)
        os.makedirs(os.path.join(out_dir, "scanned_xyz"), exist_ok=True)

        atoms = read(pdb_file)
        angle0 = _snap_to_step_nearest(float(atoms.get_dihedral(*dihedral_indices)), step_size)
        sign = -1 if direction == "cw" else 1
        n_steps = 360 // step_size
        angles = [int((angle0 + sign * i * step_size) % 360) for i in range(n_steps)]

        a = atoms.copy()
        attach_calc(a, calc)
        stepper = DihedralStepper(self.optimizer.fmax, self.optimizer.steps)

        a.set_dihedral(*dihedral_indices, angle0)
        a = stepper.optimize_with_dihedral(a, dihedral_indices, angle0)

        energies, saved_angles = [], []
        for angle in angles:
            a = stepper.optimize_with_dihedral(a, dihedral_indices, angle)
            energies.append(energy_eh(a))
            saved_angles.append(angle)
            write_pdb(os.path.join(out_dir, "scanned_pdb", f"{angle}.pdb"), a)
            write_xyz_file(os.path.join(out_dir, "scanned_xyz", f"{angle}.xyz"), a)

        np.savetxt(
            os.path.join(out_dir, "angles_vs_energies.txt"),
            np.column_stack((saved_angles, energies)),
        )
        _shift_scan_files(out_dir)


# =====================================================================
# Plotting
# =====================================================================

def plot_shifted_energy_profiles(
    base_dir: str,
    dihedral: Tuple[int, int, int, int],
    method: str,
    best_profile: Optional[str] = None,
) -> None:
    """Plot standard (non-MCS) energy profiles for a single method: NN/ML, GAFF2-old, GAFF2-best."""
    import matplotlib.pyplot as plt

    dih_str = "_".join(map(str, dihedral))
    base = os.path.join(base_dir, "scanning", dih_str)

    plt.figure(figsize=(10, 6))
    has_data = False

    qm_path = os.path.join(base, method, "angles_vs_energies_final.txt")
    if os.path.isfile(qm_path):
        data = np.loadtxt(qm_path)
        plt.plot(data[:, 0], data[:, 1] * CONV_EH_TO_KCAL_MOL, "o-", label=f"{method}")
        has_data = True

    gaff2_old_path = os.path.join(base, "GAFF2", "old", method, "dd_ee_shifted")
    if os.path.isfile(gaff2_old_path):
        data = np.loadtxt(gaff2_old_path)
        plt.plot(data[:, 0], data[:, 1], "s--", label= "GAFF2")
        has_data = True

    if best_profile and os.path.isfile(best_profile):
        data = np.loadtxt(best_profile)
        plt.plot(data[:, 0], data[:, 1], "^-.", label="GAFF2 rep")
        has_data = True

    if has_data:
        i, j, k, l = dihedral
        plt.xlabel("Dihedral Angle (degrees)")
        plt.ylabel("Relative Energy (kcal/mol)")
        plt.title(f"Dihedral Scan Energy Profile [{method}]: {i}-{j}-{k}-{l}")
        plt.legend()
        plt.xlim(0, 360)
        plt.grid(True)
        plt.tight_layout()
        out_png = os.path.join(base, f"{dih_str}_{method}.png")
        plt.savefig(out_png)
        logging.info("Saved plot: %s", out_png)
    plt.close()


def plot_shifted_energy_profiles_mcs(
    base_dir: str,
    dihedral: Tuple[int, int, int, int],
    method: str,
    best_profile: Optional[str] = None,
) -> None:
    """Plot MCS energy profiles for a single method: NN/ML-MCS, GAFF2-old-MCS, GAFF2-best."""
    import matplotlib.pyplot as plt

    dih_str = "_".join(map(str, dihedral))
    base = os.path.join(base_dir, "scanning", dih_str)

    plt.figure(figsize=(10, 6))
    has_data = False

    mcs_path = os.path.join(base, method, "MCS", "angles_vs_energies_final.txt")
    if os.path.isfile(mcs_path):
        data = np.loadtxt(mcs_path)
        plt.plot(data[:, 0], data[:, 1] * CONV_EH_TO_KCAL_MOL, "o-", label=f"{method}-MCS")
        has_data = True

    gaff2_old_path = os.path.join(base, "GAFF2", "old", method, "dd_ee_shifted")
    if os.path.isfile(gaff2_old_path):
        data = np.loadtxt(gaff2_old_path)
        plt.plot(data[:, 0], data[:, 1], "s--", label = "GAFF2")
        has_data = True

    legacy_merged = os.path.join(base, "GAFF2", "old", "dd_ee_merged_shifted")
    if os.path.isfile(legacy_merged):
        data = np.loadtxt(legacy_merged)
        plt.plot(data[:, 0], data[:, 1], "d--", label="GAFF2_MCS")
        has_data = True

    if best_profile and os.path.isfile(best_profile):
        data = np.loadtxt(best_profile)
        plt.plot(data[:, 0], data[:, 1], "^-.", label="GAFF2 rep")
        has_data = True

    if has_data:
        i, j, k, l = dihedral
        plt.xlabel("Dihedral Angle (degrees)")
        plt.ylabel("Relative Energy (kcal/mol)")
        plt.title(f"MCS Dihedral Scan Energy Profile [{method}]: {i}-{j}-{k}-{l}")
        plt.legend()
        plt.xlim(0, 360)
        plt.grid(True)
        plt.tight_layout()
        out_png = os.path.join(base, f"{dih_str}_{method}_MCS.png")
        plt.savefig(out_png)
        logging.info("Saved MCS plot: %s", out_png)
    plt.close()
