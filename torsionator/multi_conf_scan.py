"""
Multi-Conformer Scan (MCS) helpers.

These three functions implement the three steps of the MCS workflow:
  1. _generate_and_minimize_conformers  – generate, screen, and minimize conformers
  2. _scan_all_conformers               – run dihedral scans on each conformer
  3. _MCS_merge_scans                   – collapse all scans to the per-angle minimum
"""

import os
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Tuple, Union

import numpy as np
from ase.io import read

from .dihedral import DihedralScannerMCS, _norm360
from .conformers import ConformerGenerator, ConformerScreen
from .io_utils import write_xyz_file, banner


class _BestEntry(NamedTuple):
    e_compare: float
    xyz_path: str
    angle_true: float
    tag: str
    e_abs: float


def _round_to_step(a: float, step: int) -> int:
    return int(round(_norm360(a) / step) * step) % 360


def _load_angles_energies(path_txt: str):
    data = np.loadtxt(path_txt)
    if data.ndim == 1:
        if data.size == 0:
            return np.array([]), np.array([])
        if data.size == 2:
            return data[[0]], data[[1]]
        raise ValueError(f"Unexpected data format in {path_txt}: shape={data.shape}")
    return data[:, 0], data[:, 1]


# =====================================================================
# Step 1: generate & minimize conformers
# =====================================================================

def _generate_and_minimize_conformers(
    self,
    pdb_file: str,
    methods: List[str],
    n_confs: int,
    threshold: float,
) -> Dict[str, List[str]]:
    banner(self.log, "STEP 1: generate & minimize (all conformers)")

    gen = ConformerGenerator(os.path.join(self.cfg.base_dir, "conformers"))
    conf_dir = gen.generate(pdb_file, n_confs)

    screen = ConformerScreen()
    screen.screen_by_rmsd(conf_dir, threshold)
    screen.ensure_ref_if_unique(conf_dir, pdb_file, threshold)
    gen.normalize_conformer(conf_dir)

    def _sort_key(p: Path) -> Tuple[int, Union[int, str]]:
        return (0, int(p.stem)) if p.stem.isdigit() else (1, p.stem)

    conf_paths = sorted(Path(conf_dir).glob("*.pdb"), key=_sort_key)
    per_method_conf_dirs: Dict[str, List[str]] = {m: [] for m in methods}

    for conf_pdb in conf_paths:
        conf_id = conf_pdb.stem
        for m in methods:
            calc, _ = self.calc_factory.get(m)
            target_dir = os.path.join(self.cfg.base_dir, "conformers", m, conf_id)
            os.makedirs(target_dir, exist_ok=True)
            dst = os.path.join(target_dir, "minimized.pdb")
            self.log.info("[%s][%s] minimizing -> %s", m, conf_id, os.path.relpath(dst, self.cfg.base_dir))
            self.opt.minimize_and_write_pdb(str(conf_pdb), calc, dst)
            per_method_conf_dirs[m].append(target_dir)

    for m in methods:
        self.log.info("[%s] MCS: %d conformers prepared", m, len(per_method_conf_dirs[m]))

    return per_method_conf_dirs


# =====================================================================
# Step 2: dihedral scan on every conformer
# =====================================================================

def _scan_all_conformers(
    self,
    per_method_conf_dirs: Dict[str, List[str]],
    methods: List[str],
    dihedrals_to_scan: List[Tuple[int, int, int, int]],
    step_size: int = 10,
    double_rotation: bool = False,
) -> Dict[str, Dict[Tuple[int, int, int, int], List[str]]]:
    banner(self.log, "STEP 2: Dihedral scan (all conformers)")

    scanner = DihedralScannerMCS(self.cfg.base_dir, self.opt)
    scanned_tags: Dict[str, Dict[Tuple[int, int, int, int], List[str]]] = {
        m: defaultdict(list) for m in methods
    }

    def _dir_sort_key(path: str) -> Tuple[int, Union[int, str]]:
        name = os.path.basename(path)
        return (0, int(name)) if name.isdigit() else (1, name)

    for m in methods:
        calc, label = self.calc_factory.get(m)
        conf_dirs = sorted(per_method_conf_dirs.get(m, []), key=_dir_sort_key)
        self.log.info("[%s] scanning %d conformers", m, len(conf_dirs))

        for d_tup in dihedrals_to_scan:
            for cdir in conf_dirs:
                pdb_src = os.path.join(cdir, "minimized.pdb")
                if not os.path.exists(pdb_src):
                    self.log.warning("[%s][%s] missing %s", m, d_tup, pdb_src)
                    continue

                conf_id = os.path.basename(cdir)
                tag = f"{label}/{conf_id}"
                self.log.info("[%s][%s] SCAN (%s) input=%s", m, d_tup, tag, pdb_src)

                t0 = time.time()
                scanner.scan(pdb_src, calc, d_tup, tag, direction="cw", step_size=step_size)
                self.log.info("[%s][%s] SCAN done in %.2f s", m, d_tup, time.time() - t0)
                scanned_tags[m][d_tup].append(tag)

                if double_rotation:
                    tag_ccw = f"{label}/{conf_id}_ccw"
                    self.log.info("[%s][%s] SCAN (ccw) (%s) input=%s", m, d_tup, tag_ccw, pdb_src)
                    t1 = time.time()
                    scanner.scan(pdb_src, calc, d_tup, tag_ccw, direction="ccw", step_size=step_size)
                    self.log.info("[%s][%s] SCAN (ccw) done in %.2f s", m, d_tup, time.time() - t1)
                    scanned_tags[m][d_tup].append(tag_ccw)

    return scanned_tags


# =====================================================================
# Step 3: merge all conformer scans → single MCS profile
# =====================================================================

def _MCS_merge_scans(
    self,
    base_dir: str,
    dihedral_indices: Tuple[int, int, int, int],
    out_method_name: str,
    scanned_tags: List[str],
    step_size: int = 10,
    mode: str = "absolute",
) -> str:
    """
    Per angle bin, keep the geometry with the lowest energy across all conformer scans.

    Output directory: scanning/<a_b_c_d>/<out_method_name>/MCS/
    Files written:
      - geometries.xyz
      - energies.dat                   (absolute energies shifted to min=0)
      - angles_vs_energies.txt         (angle  absolute_energy)
      - angles_vs_energies_final.txt   (angle  relative_energy)
    """
    assert mode in ("absolute", "shifted"), f"Unsupported mode: {mode!r}"
    banner(self.log, "STEP 3: Merge MCS scans")

    dih_str = "_".join(map(str, dihedral_indices))
    scan_root = os.path.join(base_dir, "scanning", dih_str)
    out_dir = os.path.join(scan_root, out_method_name, "MCS")
    os.makedirs(out_dir, exist_ok=True)

    angle_bins = list(range(0, 360, step_size))
    best: Dict[int, Optional[_BestEntry]] = {b: None for b in angle_bins}

    for tag in scanned_tags:
        tag_dir = os.path.join(scan_root, tag)
        txt_path = os.path.join(tag_dir, "angles_vs_energies.txt")
        xyz_dir = os.path.join(tag_dir, "scanned_xyz")

        if not (os.path.exists(txt_path) and os.path.isdir(xyz_dir)):
            continue

        angles_deg, energies_eh = _load_angles_energies(txt_path)
        if energies_eh.size == 0:
            continue

        if mode == "shifted":
            energies_compare = energies_eh - energies_eh.min()
        else:
            energies_compare = energies_eh

        for ang, E_cmp, E_abs in zip(angles_deg, energies_compare, energies_eh):
            ang_norm = _norm360(ang)
            bin_angle = _round_to_step(ang_norm, step_size)
            deg_name = int(round(ang_norm)) % 360

            # Locate the xyz file, allowing ±3° tolerance
            xyz_path = os.path.join(xyz_dir, f"{deg_name}.xyz")
            if not os.path.exists(xyz_path):
                xyz_path = next(
                    (
                        os.path.join(xyz_dir, f"{(deg_name + d) % 360}.xyz")
                        for d in (1, -1, 2, -2, 3, -3)
                        if os.path.exists(os.path.join(xyz_dir, f"{(deg_name + d) % 360}.xyz"))
                    ),
                    None,
                )
                if xyz_path is None:
                    continue

            cur = best[bin_angle]
            if cur is None or E_cmp < cur.e_compare - 1e-15:
                best[bin_angle] = _BestEntry(E_cmp, xyz_path, ang_norm, tag, E_abs)

    chosen = [(b, v) for b, v in sorted(best.items()) if v is not None]
    if not chosen:
        return out_dir

    xyz_out = os.path.join(out_dir, "geometries.xyz")
    energies_abs: List[float] = []
    angles_kept: List[int] = []

    for b, entry in chosen:
        frame = read(entry.xyz_path)
        header = f"angle={entry.angle_true:.2f} E={entry.e_abs:.8f} conformer={entry.tag}"
        write_xyz_file(xyz_out, frame, mode="a", sig=16, header=header)
        energies_abs.append(entry.e_abs)
        angles_kept.append(b)

    energies_arr = np.array(energies_abs)
    energies_rel = energies_arr - energies_arr.min()

    np.savetxt(os.path.join(out_dir, "energies.dat"), energies_rel, fmt="%.12f")
    np.savetxt(
        os.path.join(out_dir, "angles_vs_energies.txt"),
        np.column_stack([angles_kept, energies_arr]),
        fmt=["%d", "%.12f"],
    )
    np.savetxt(
        os.path.join(out_dir, "angles_vs_energies_final.txt"),
        np.column_stack([angles_kept, energies_rel]),
        fmt=["%d", "%.12f"],
    )

    return out_dir
