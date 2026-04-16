"""
AMBER/GAFF2 refinement utilities.

Public API:
  - find_best_hrst_by_RMSE
  - get_atom_types_from_mol2
  - param_input_file
  - update_frcmod
  - expand_dihedral_in_frcmod
  - Gaff2SingleDihedralRefinement
  - Gaff2NewDihedralRefinement
"""

import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
from ase.io import read, write

from ..config import Config


# ============================================================================
# RMSD / energy helpers
# ============================================================================

def calculate_rmsd(energy_nn: np.ndarray, energy_mm: np.ndarray) -> float:
    if len(energy_nn) != len(energy_mm):
        raise ValueError(
            f"Energy arrays have different lengths: {len(energy_nn)} vs {len(energy_mm)}"
        )
    return float(np.sqrt(np.mean((energy_nn - energy_mm) ** 2)))


def load_energy_file(filepath: str) -> Tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(filepath)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data[:, 0], data[:, 1]


def find_best_hrst_by_RMSE(
    nn_energy_file: str,
    gaff2_new_dir: str,
    log,
) -> Tuple[float, str, float, str]:
    """
    Compare NN energies against each hrst scan in GAFF2/new/<method>/<hrst>/.
    Returns (best_hrst, path_best_frcmod, best_rmsd, best_kcal_profile).
    """
    nn_angles, nn_energies = load_energy_file(nn_energy_file)

    hrst_dirs = []
    for item in os.listdir(gaff2_new_dir):
        item_path = os.path.join(gaff2_new_dir, item)
        if os.path.isdir(item_path):
            try:
                hrst_dirs.append((float(item), item_path))
            except ValueError:
                continue

    if not hrst_dirs:
        raise FileNotFoundError(f"No hrst directories found in {gaff2_new_dir}")

    hrst_dirs.sort()
    rmsd_results: Dict[float, dict] = {}

    for hrst_val, hrst_path in hrst_dirs:
        hartree_file = os.path.join(hrst_path, "dd_ee_shifted_hartree")
        frcmod_file = os.path.join(hrst_path, "antechamber", "MOL.frcmod")

        mm_angles, mm_energies = load_energy_file(hartree_file)
        if len(mm_energies) != len(nn_energies):
            mm_energies = np.interp(nn_angles, mm_angles, mm_energies)

        rmsd_results[hrst_val] = {
            "rmsd": calculate_rmsd(nn_energies, mm_energies),
            "frcmod": frcmod_file,
            "kcal_profile": os.path.join(hrst_path, "dd_ee_shifted"),
        }

    best_hrst = min(rmsd_results, key=lambda x: rmsd_results[x]["rmsd"])
    log.info("[RMSD] BEST hrst=%s RMSD=%.6f", best_hrst, rmsd_results[best_hrst]["rmsd"])

    return (
        best_hrst,
        rmsd_results[best_hrst]["frcmod"],
        float(rmsd_results[best_hrst]["rmsd"]),
        rmsd_results[best_hrst]["kcal_profile"],
    )


# ============================================================================
# Angle grid helpers
# ============================================================================

def _validate_step_size(step_size: int) -> int:
    try:
        s = int(step_size)
    except Exception:
        raise ValueError(f"step_size must be int-like, got {step_size!r}")
    if s <= 0:
        raise ValueError(f"step_size must be > 0, got {s}")
    if 360 % s != 0:
        raise ValueError(f"step_size must divide 360 exactly. Got {s}.")
    return s


def _norm360(angle: float) -> float:
    return float(((angle % 360.0) + 360.0) % 360.0)


def nearest_multiple(angle_deg: float, step: int) -> int:
    step = _validate_step_size(step)
    return int(round(_norm360(angle_deg) / step) * step) % 360


def angles_0_350(step: int = 10) -> List[int]:
    return list(range(0, 360, _validate_step_size(step)))


def angles_from_angle0(angle0: int, step: int) -> List[int]:
    step = _validate_step_size(step)
    a0 = int(angle0) % 360
    return [int((a0 + k * step) % 360) for k in range(360 // step)]


# ============================================================================
# Small helpers
# ============================================================================

def strip_conect_records(pdb_in: str) -> str:
    with open(pdb_in) as f:
        lines = f.readlines()
    with open(pdb_in, "w") as f:
        for line in lines:
            if not line.startswith("CONECT"):
                f.write(line)
    return pdb_in


def get_atom_types_from_mol2(mol2_file: str) -> Dict[int, str]:
    atom_types: Dict[int, str] = {}
    with open(mol2_file) as f:
        lines = f.readlines()
    in_atoms = False
    for line in lines:
        if line.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            continue
        if line.startswith("@<TRIPOS>BOND"):
            break
        if in_atoms:
            parts = line.split()
            if len(parts) >= 6:
                atom_types[int(parts[0]) - 1] = parts[5]
    return atom_types


# ============================================================================
# MCS: xyz → per-angle folders
# ============================================================================

def write_mcs_frames_to_angle_folders(
    mcs_xyz: str,
    out_root: str,
    angles: Optional[Iterable[int]] = None,
    pdb_name: str = "angle.pdb",
) -> List[int]:
    """Split a multi-frame XYZ into out_root/<angle>/angle.pdb."""
    os.makedirs(out_root, exist_ok=True)
    angs = list(angles) if angles is not None else angles_0_350(10)
    frames = read(mcs_xyz, index=":")

    if len(frames) < len(angs):
        raise RuntimeError(
            f"MCS xyz has {len(frames)} frames but expected at least {len(angs)}"
        )

    for k, ang in enumerate(angs):
        a_dir = os.path.join(out_root, str(int(ang)))
        os.makedirs(a_dir, exist_ok=True)
        write(os.path.join(a_dir, pdb_name), frames[k], format="proteindatabank")

    return [int(a) for a in angs]


def cleanup_mcs_angle_outputs(root_dir: str, step: int = 10, log=None) -> None:
    """Remove angle folders and dd_ee* files from a previous MCS run."""
    step = _validate_step_size(step)
    for fn in ("dd_ee", "dd_ee_shifted", "dd_ee_shifted_hartree"):
        p = os.path.join(root_dir, fn)
        if os.path.exists(p):
            try:
                os.remove(p)
            except Exception as e:
                if log:
                    log.warning("[CLEAN][MCS] cannot remove %s: %s", p, e)
    for ang in angles_0_350(step):
        d = os.path.join(root_dir, str(int(ang)))
        if os.path.isdir(d):
            try:
                shutil.rmtree(d)
            except Exception as e:
                if log:
                    log.warning("[CLEAN][MCS] cannot remove dir %s: %s", d, e)


# ============================================================================
# Sander output parsers
# ============================================================================


def extract_sp_energy_from_nstep(out_file: str) -> float:
    """Parse ENERGY from the table following NSTEP/ENERGY header (MCS single-point path)."""
    with open(out_file) as f:
        lines = f.readlines()

    last_header = next(
        (i for i in reversed(range(len(lines))) if "NSTEP" in lines[i] and "ENERGY" in lines[i]),
        None,
    )
    if last_header is None:
        raise RuntimeError(f"No NSTEP/ENERGY table found in {out_file}")

    j = last_header + 1
    while j < len(lines) and not lines[j].strip():
        j += 1
    if j >= len(lines):
        raise RuntimeError(f"NSTEP header found but no data line in {out_file}")

    parts = lines[j].split()
    if len(parts) < 2:
        raise RuntimeError(f"Cannot parse ENERGY from line: {lines[j].rstrip()}")

    try:
        return float(parts[1].replace("D", "E"))
    except ValueError as e:
        raise RuntimeError(f"Could not convert ENERGY='{parts[1]}' to float in {out_file}") from e


def write_shifted_profiles(dd_ee_path: str, hartree_factor: float = 0.0015936) -> Tuple[str, str]:
    """Given dd_ee (angle, energy), write dd_ee_shifted (kcal) and dd_ee_shifted_hartree."""
    data = np.loadtxt(dd_ee_path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    ang = data[:, 0]
    shifted = data[:, 1] - data[:, 1].min()

    out_dir = os.path.dirname(dd_ee_path)
    kcal_path = os.path.join(out_dir, "dd_ee_shifted")
    eh_path = os.path.join(out_dir, "dd_ee_shifted_hartree")
    np.savetxt(kcal_path, np.column_stack([ang, shifted]), fmt=["%g", "%.8f"])
    np.savetxt(eh_path, np.column_stack([ang, shifted * hartree_factor]), fmt=["%g", "%.10f"])
    return kcal_path, eh_path


# ============================================================================
# frcmod: dihedral expansion (basis terms)
# ============================================================================

_DIH_LINE_RE = re.compile(
    r"""^(?P<a1>.{2})-(?P<a2>.{2})-(?P<a3>.{2})-(?P<a4>.{2})
        (?P<div>\s+[-\d\.Ee+]+)
        (?P<bar>\s+[-\d\.Ee+]+)
        (?P<ph>\s+[-\d\.Ee+]+)
        (?P<per>\s+[-\d\.Ee+]+)
        (?P<rest>.*)$""",
    re.VERBOSE,
)


def _parse_frcmod_dihe_line(line: str) -> Optional[dict]:
    m = _DIH_LINE_RE.match(line.rstrip("\n"))
    if not m:
        return None
    key = f"{m.group('a1')}-{m.group('a2')}-{m.group('a3')}-{m.group('a4')}"
    try:
        return {
            "key": key,
            "divisor": float(m.group("div")),
            "barrier": float(m.group("bar")),
            "phase": float(m.group("ph")),
            "periodicity": float(m.group("per")),
            "rest": m.group("rest"),
            "raw": line if line.endswith("\n") else line + "\n",
        }
    except ValueError:
        return None


def _format_frcmod_dihe_line(
    key: str, divisor: float, barrier: float, phase: float, periodicity: float, rest: str = ""
) -> str:
    line = f"{key:<11}{divisor:>4.0f}{barrier:>9.3f}{phase:>14.3f}{periodicity:>16.3f}"
    if rest:
        line += rest.rstrip("\n")
    return line + "\n"


def _canon_per(per: float) -> int:
    p = int(round(float(per)))
    if p == 0:
        return 0
    if abs(p) == 1:
        return 1
    return -abs(p)


def _std2(t: str) -> str:
    return f"{str(t).strip():<2}"[:2]


def _reverse_key(key: str) -> str:
    return "-".join(reversed(key.split("-")))


def _find_dihe_section(lines: List[str]) -> Tuple[Optional[int], Optional[int]]:
    _SECTION_ENDS = {"IMPROPER", "IMPR", "NONB", "MASS", "BOND", "ANGLE", "HBON"}
    for i, ln in enumerate(lines):
        if ln.strip().upper() == "DIHE":
            j = i + 1
            while j < len(lines) and lines[j].strip().upper() not in _SECTION_ENDS:
                j += 1
            return i, j
    return None, None


def _index_frcmod_dihe_keys(lines: List[str]) -> Dict[str, List[int]]:
    start, end = _find_dihe_section(lines)
    if start is None:
        return {}
    key_to_indices: Dict[str, List[int]] = {}
    for idx, ln in enumerate(lines[start + 1: end]):
        pr = _parse_frcmod_dihe_line(ln)
        if pr:
            key_to_indices.setdefault(str(pr["key"]), []).append(idx)
    return key_to_indices


def dihe_candidate_keys_from_atoms(mol2_path: str, d_tup: Tuple[int, int, int, int]) -> List[str]:
    atom_types = get_atom_types_from_mol2(mol2_path)
    a1, a2, a3, a4 = (_std2(atom_types[i]) for i in d_tup)
    X = _std2("X")
    candidates = [
        f"{a1}-{a2}-{a3}-{a4}",
        f"{a4}-{a3}-{a2}-{a1}",
        f"{X}-{a2}-{a3}-{X}",
        f"{X}-{a3}-{a2}-{X}",
        f"{X}-{a2}-{a3}-{a4}",
        f"{a1}-{a2}-{a3}-{X}",
        _reverse_key(f"{X}-{a2}-{a3}-{a4}"),
        _reverse_key(f"{a1}-{a2}-{a3}-{X}"),
    ]
    # deduplicate preserving order
    seen, out = set(), []
    for k in candidates:
        if k not in seen:
            seen.add(k)
            out.append(k)
    return out


def find_dihedral_key_in_frcmod(
    frcmod_path: str, candidates: List[str], log=None
) -> Optional[str]:
    with open(frcmod_path) as f:
        lines = f.readlines()
    key_to_indices = _index_frcmod_dihe_keys(lines)
    if not key_to_indices:
        if log:
            log.warning("[FIND-DIHE] No DIHE section in %s", frcmod_path)
        return None
    for k in candidates:
        if k in key_to_indices:
            if log:
                log.info("[FIND-DIHE] dihedral=%r", k)
            return k
    for k in candidates:
        rk = _reverse_key(k)
        if rk in key_to_indices:
            if log:
                log.info("[FIND-DIHE] dihedral=%r (reversed from %r)", rk, k)
            return rk
    if log:
        log.warning(
            "[FIND-DIHE] NO MATCH. candidates=%s | sample keys=%s",
            [repr(x) for x in candidates],
            [repr(x) for x in list(key_to_indices.keys())[:20]],
        )
    return None


def read_basis_from_frcmod(
    frcmod_path: str,
    mol2_path: str,
    d_tup: Tuple[int, int, int, int],
    log=None,
) -> Tuple[int, ...]:
    """Return the periodicity terms present in *frcmod_path* for the dihedral
    defined by *d_tup* (atom-type lookup via *mol2_path*).

    The returned tuple uses the sign convention of ``_basis_from_max_mult``:
    negative for p > 1 (indicating more terms follow in AMBER convention),
    positive for p = 1.  Terms are ordered by descending absolute value with
    p = 1 last, matching the layout produced by ``_basis_from_max_mult``.
    """
    candidates = dihe_candidate_keys_from_atoms(mol2_path, d_tup)
    key = find_dihedral_key_in_frcmod(frcmod_path, candidates, log=log)
    if key is None:
        raise RuntimeError(
            f"Cannot detect basis periodicities: dihedral {d_tup} not found in {frcmod_path}"
        )
    with open(frcmod_path) as f:
        lines = f.readlines()
    start, end = _find_dihe_section(lines)
    if start is None:
        raise RuntimeError(f"No DIHE section in {frcmod_path}")
    periodicities: set = set()
    for ln in lines[start + 1: end]:
        pr = _parse_frcmod_dihe_line(ln)
        if pr and pr["key"] == key:
            periodicities.add(_canon_per(pr["periodicity"]))
    if not periodicities:
        raise RuntimeError(
            f"No periodicity terms found for key {key!r} in {frcmod_path}"
        )
    non_one = sorted((p for p in periodicities if p != 1), key=abs, reverse=True)
    has_one = 1 in periodicities
    return tuple(non_one + ([1] if has_one else []))


def expand_dihedral_in_frcmod(
    frcmod_path: str,
    mol2_path: str,
    d_tup: Tuple[int, int, int, int],
    basis_periodicities: Tuple[int, ...],
    frcmod_out: Optional[str] = None,
    log=None,
) -> str:
    if frcmod_out is None:
        frcmod_out = frcmod_path

    candidates = dihe_candidate_keys_from_atoms(mol2_path, d_tup)
    key_in_frcmod = find_dihedral_key_in_frcmod(frcmod_path, candidates, log=log)
    if key_in_frcmod is None:
        raise RuntimeError(
            f"Cannot find dihedral for atoms {d_tup} in {frcmod_path}. Tried: {candidates}"
        )

    # Deduplicate and canonicalize periodicities
    seen, basis_canon = set(), []
    for b in basis_periodicities:
        bi = _canon_per(b)
        if bi not in seen:
            basis_canon.append(bi)
            seen.add(bi)

    with open(frcmod_path) as f:
        orig = f.readlines()

    start, end = _find_dihe_section(orig)
    if start is None:
        raise RuntimeError(f"No DIHE section in {frcmod_path}")

    pre = orig[: start + 1]
    body = list(orig[start + 1: end])
    post = orig[end:]

    # Parse existing entries for this key
    key_idxs, parsed = [], []
    for idx, ln in enumerate(body):
        pr = _parse_frcmod_dihe_line(ln)
        parsed.append(pr)
        if pr and str(pr["key"]) == key_in_frcmod:
            key_idxs.append(idx)

    per_to_params: Dict[int, dict] = {}
    seed_pr = None
    for idx in key_idxs:
        pr = parsed[idx]
        if pr:
            if seed_pr is None:
                seed_pr = pr
            per_to_params[_canon_per(pr["periodicity"])] = pr

    seed_div = float(seed_pr["divisor"]) if seed_pr else 1.0
    seed_phase = float(seed_pr["phase"]) if seed_pr else 0.0
    seed_rest = str(seed_pr.get("rest", "")) if seed_pr else ""

    # Rebuild body without the old key lines
    new_body = [ln for i, ln in enumerate(body) if i not in set(key_idxs)]

    rebuilt = []
    for per_i in basis_canon:
        if per_i in per_to_params:
            pr = per_to_params[per_i]
            rebuilt.append(_format_frcmod_dihe_line(
                key_in_frcmod, float(pr["divisor"]), float(pr["barrier"]),
                float(pr["phase"]), float(per_i), str(pr.get("rest", ""))
            ))
        else:
            rebuilt.append(_format_frcmod_dihe_line(
                key_in_frcmod, seed_div, 0.0, seed_phase, float(per_i), seed_rest
            ))

    insert_at = min(key_idxs)
    new_body[insert_at:insert_at] = rebuilt

    with open(frcmod_out, "w") as f:
        f.writelines(pre + new_body + post)

    return frcmod_out


# ============================================================================
# MDGX param input
# ============================================================================

def param_input_file(
    atom_types: Dict[int, str],
    method: str,
    dihedral_indices: Tuple[int, int, int, int],
    scan_dir: str,
    method_subdir: str,
    hrst_val: float,
    use_old_expanded: bool = True,
) -> Tuple[str, str, str]:
    dstr = "_".join(map(str, dihedral_indices))
    dih_types = " ".join(atom_types[i] for i in dihedral_indices)

    if use_old_expanded:
        leap_top = os.path.abspath(os.path.join(scan_dir, f"GAFF2/old/{method}/leap_expanded/MOL.top"))
        fallback_frcmod = os.path.abspath(os.path.join(scan_dir, f"GAFF2/old/{method}/antechamber_expanded/MOL.frcmod"))
    else:
        leap_top = os.path.abspath(os.path.join(scan_dir, f"GAFF2/old/{method}/leap/MOL.top"))
        fallback_frcmod = os.path.abspath(os.path.join(scan_dir, f"GAFF2/old/{method}/antechamber/MOL.frcmod"))

    traj = os.path.abspath(os.path.join(scan_dir, method_subdir, "geometries.nc"))
    energy_file = os.path.abspath(os.path.join(scan_dir, method_subdir, "energies.dat"))
    output_dir = os.path.join(scan_dir, method_subdir, "output", str(hrst_val))
    os.makedirs(output_dir, exist_ok=True)

    local_frcmod = os.path.join(output_dir, "MOL.frcmod")
    frcmod = os.path.abspath(local_frcmod if os.path.isfile(local_frcmod) else fallback_frcmod)

    for path, name in [(energy_file, "energies.dat"), (traj, "geometries.nc"), (frcmod, "frcmod"), (leap_top, "MOL.top")]:
        if not os.path.exists(path):
            raise FileNotFoundError(f"Missing {name} at {path}")

    out_dat = os.path.join(output_dir, f"{hrst_val}.dat")
    out_out = os.path.join(output_dir, f"{hrst_val}.out")
    pfile = os.path.join(output_dir, f"parametrization_{method}_dihedral_{dstr}.in")

    with open(pfile, "w") as f:
        f.write("&files\n")
        f.write(f"  -p      {leap_top}\n")
        f.write(f"  -parm   {frcmod}\n")
        f.write(f"  -d      {out_dat}\n")
        f.write(f"  -o      {out_out}\n")
        f.write("&end\n\n&param\n")
        f.write(f"System {leap_top} {traj} {energy_file}\n\n")
        f.write("  ParmOutput frcmod\n")
        f.write("  eunits    hartree,\n\n")
        f.write("  verbose    1,\n\n")
        f.write("  % Torsion fitting input\n")
        f.write(f"  fith {dih_types}\n")
        f.write(f"  hrst       {hrst_val},\n&end\n")

    return pfile, out_dat, out_out


# ============================================================================
# frcmod merge (mdgx .dat → frcmod)
# ============================================================================

def _key_from_dat(line: str) -> Optional[str]:
    m = re.match(
        r"^([A-Za-z0-9]{1,2})\s*-\s*([A-Za-z0-9]{1,2})\s*-\s*([A-Za-z0-9]{1,2})\s*-\s*([A-Za-z0-9]{1,2})",
        line.strip(),
    )
    if not m:
        return None
    atoms = [f"{m.group(i):<2}" for i in range(1, 5)]
    return "-".join(atoms)


def _parse_dat_line(line: str) -> Tuple[Optional[str], Optional[dict]]:
    key = _key_from_dat(line)
    if not key:
        return None, None
    parts = line.strip().split()
    p0 = next((i for i, p in enumerate(parts) if _is_float(p)), None)
    if p0 is None or len(parts) < p0 + 4:
        return None, None
    return key, {
        "divisor": float(parts[p0]),
        "barrier": float(parts[p0 + 1]),
        "phase": float(parts[p0 + 2]),
        "periodicity": float(parts[p0 + 3]),
    }


def _is_float(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


def _format_frcmod_line(key: str, params: dict) -> str:
    return (
        f"{key:<11}"
        f"{params['divisor']:>4.0f}"
        f"{params['barrier']:>9.3f}"
        f"{params['phase']:>14.3f}"
        f"{params['periodicity']:>16.3f}"
    )


def update_frcmod(
    frcmod_path: str,
    dat_path: str,
    method: str,
    output_folder: str,
    replace_dat_keys: bool = False,
) -> str:
    """Merge mdgx .dat dihedral terms into an frcmod file."""
    grouped: Dict[str, Dict[int, dict]] = {}
    key_order: List[str] = []

    with open(dat_path) as f:
        dat_lines = f.readlines()

    in_dihe = False
    for line in dat_lines:
        s = line.strip()
        if not s:
            continue
        if s.upper() == "DIHE":
            in_dihe = True
            continue
        if in_dihe and s.upper().split()[0] in ("IMPROPER", "IMPR", "NONB", "MASS", "BOND", "ANGLE", "HBON"):
            in_dihe = False
            continue
        if not in_dihe:
            continue
        k, p = _parse_dat_line(line)
        if k and p:
            per_i = int(round(float(p["periodicity"])))
            if k not in grouped:
                grouped[k] = {}
                key_order.append(k)
            grouped[k][per_i] = p

    os.makedirs(output_folder, exist_ok=True)
    inplace = os.path.dirname(os.path.abspath(frcmod_path)) == os.path.abspath(output_folder)
    out_file = frcmod_path if inplace else os.path.join(output_folder, f"new_{method}.frcmod")

    if not grouped:
        if not inplace:
            shutil.copy2(frcmod_path, out_file)
        return out_file

    with open(frcmod_path) as f:
        orig = f.readlines()

    start, end = _find_dihe_section(orig)

    if start is None:
        # Append a new DIHE section
        contents = list(orig)
        if contents and not contents[-1].endswith("\n"):
            contents[-1] += "\n"
        contents.append("DIHE\n")
        for k in key_order:
            for per_i, p in grouped[k].items():
                contents.append(_format_frcmod_line(k, p) + "\n")
        with open(out_file, "w") as f:
            f.writelines(contents)
        return out_file

    pre = orig[: start + 1]
    body = list(orig[start + 1: end])
    post = orig[end:]

    parsed_body = [_parse_frcmod_dihe_line(ln) for ln in body]
    key_to_indices: Dict[str, List[int]] = {}
    for idx, pr in enumerate(parsed_body):
        if pr:
            key_to_indices.setdefault(str(pr["key"]), []).append(idx)

    def _update_key(key_fwd: str) -> None:
        key_rev = _reverse_key(key_fwd)
        key_in_frcmod = (
            key_fwd if key_fwd in key_to_indices
            else key_rev if key_rev in key_to_indices
            else None
        )
        terms = grouped[key_fwd]

        if key_in_frcmod is None:
            for per_i, p in terms.items():
                new_line = _format_frcmod_line(key_fwd, p) + "\n"
                body.append(new_line)
                pr = _parse_frcmod_dihe_line(new_line)
                parsed_body.append(pr)
                if pr:
                    key_to_indices.setdefault(str(pr["key"]), []).append(len(body) - 1)
            return

        if replace_dat_keys:
            idxs_to_remove = sorted(key_to_indices[key_in_frcmod], reverse=True)
            insert_pos = min(key_to_indices[key_in_frcmod])
            for bi in idxs_to_remove:
                body.pop(bi)
                parsed_body.pop(bi)
            key_to_indices.clear()
            for idx, pr in enumerate(parsed_body):
                if pr:
                    key_to_indices.setdefault(str(pr["key"]), []).append(idx)
            for offset, (per_i, p) in enumerate(sorted(terms.items())):
                new_line = _format_frcmod_line(key_fwd, p) + "\n"
                pos = insert_pos + offset
                body.insert(pos, new_line)
                parsed_body.insert(pos, _parse_frcmod_dihe_line(new_line))
            key_to_indices.clear()
            for idx, pr in enumerate(parsed_body):
                if pr:
                    key_to_indices.setdefault(str(pr["key"]), []).append(idx)
            return

        per_to_idx = {_canon_per(parsed_body[bi]["periodicity"]): bi
                      for bi in key_to_indices.get(key_in_frcmod, [])
                      if parsed_body[bi]}

        for per_i, p in terms.items():
            if per_i in per_to_idx:
                bi = per_to_idx[per_i]
                rest = str(parsed_body[bi].get("rest", "")) if parsed_body[bi] else ""
                body[bi] = _format_frcmod_dihe_line(
                    key_in_frcmod, p["divisor"], p["barrier"], p["phase"], p["periodicity"], rest
                )
                parsed_body[bi] = _parse_frcmod_dihe_line(body[bi])
            else:
                insert_at = max(key_to_indices[key_in_frcmod]) + 1
                new_line = _format_frcmod_line(key_in_frcmod, p) + "\n"
                body.insert(insert_at, new_line)
                parsed_body.insert(insert_at, _parse_frcmod_dihe_line(new_line))
                for k2 in key_to_indices:
                    key_to_indices[k2] = [x + 1 if x >= insert_at else x for x in key_to_indices[k2]]
                key_to_indices.setdefault(key_in_frcmod, []).append(insert_at)

    for k in key_order:
        _update_key(k)

    with open(out_file, "w") as f:
        f.writelines(pre + body + post)

    return out_file


# ============================================================================
# Shared base for Gaff2*DihedralRefinement
# ============================================================================

_MIN_TEMPLATE = """\
minimization
 &cntrl
   imin=1,
   ntx=1, irest=0,
   ntxo=2, ntpr=10, ntwr=10000000, ntwx=10, ioutfm=1,
   maxcyc=10000, ntmin=3,
   ntc=1, ntf=1, ntb=0, cut=999.,
   igb=6, rgbmax=999.,
   nmropt=1,
/
&wt type='DUMPFREQ', istep1=1 /
&wt type='END' /
DISANG=restraints-XXX
DUMPAVE=dihedral-XXX.dat
"""

_MIN0_TEMPLATE = """\
minimization
 &cntrl
   imin=1,
   ntx=1, irest=0,
   ntxo=2, ntpr=10, ntwr=10000000, ntwx=10, ioutfm=1,
   maxcyc=10000, ntmin=3,
   ntc=1, ntf=1, ntb=0, cut=999.,
   igb=6, rgbmax=999.,
/
"""

_SP_TEMPLATE = """\
sp
&cntrl
  imin=1,
  maxcyc=0,
  ntmin=1,
  ntb=0,
  cut=999.0,
  igb=6,
  ntpr=1,
/
"""

_SCR_EN_TEMPLATE = """\
#!/usr/bin/env bash
set -e

ANGLES=({ang_list})

for a in "${{ANGLES[@]}}"; do
  tail -n 1 "dihedral-$a.dat"
done > dd

for a in "${{ANGLES[@]}}"; do
  grep EAMBER "min-$a.out" | tail -n 1
done > ee

paste dd ee | awk '{{print $2, $5}}' > dd_ee
rm dd ee

awk '
NR==1 {{min=$2}}
NR>1  {{if ($2 < min) min=$2}}
{{angles[NR]=$1; energies[NR]=$2}}
END {{
  for (i=1;i<=NR;i++) {{
    printf "%g %f\\n", angles[i], energies[i]-min
  }}
}}
' dd_ee > dd_ee_shifted

awk '
NR==1 {{min=$2}}
NR>1  {{if ($2 < min) min=$2}}
{{angles[NR]=$1; energies[NR]=$2}}
END {{
  factor = 0.0015936
  for (i=1;i<=NR;i++) {{
    printf "%g %f\\n", angles[i], (energies[i]-min)*factor
  }}
}}
' dd_ee > dd_ee_shifted_hartree
"""


@dataclass
class _AmberRefinementBase:
    cfg: Config
    log: object
    tleap_cmd: str = "tleap"
    sander_cmd: str = "sander"

    # ---- file writers ----

    def _write_min_template(self, gdir: str) -> None:
        with open(os.path.join(gdir, "min.in.template"), "w") as f:
            f.write(_MIN_TEMPLATE)

    def _write_restrain_template(self, gdir: str, d_tup: Tuple[int, int, int, int]) -> None:
        a, b, c, d = (x + 1 for x in d_tup)
        content = f"# bias potential\n&rst  iat={a},{b},{c},{d}\n    r1=-1000,  r2=XXX,  r3=XXX,  r4=1000.,  rk2=200.,  rk3=200., /\n"
        with open(os.path.join(gdir, "restraints.template"), "w") as f:
            f.write(content)

    def _write_min0_file(self, gdir: str) -> None:
        with open(os.path.join(gdir, "min.in"), "w") as f:
            f.write(_MIN0_TEMPLATE)

    def _write_sp_in(self, out_dir: str) -> None:
        with open(os.path.join(out_dir, "sp.in"), "w") as f:
            f.write(_SP_TEMPLATE)

    def _write_scr_min(self, gdir: str) -> None:
        path = os.path.join(gdir, "scr_min.sh")
        with open(path, "w") as f:
            f.write(f"#!/usr/bin/env bash\nset -e\n{self.sander_cmd} -O -p leap/MOL.top -i min.in -c leap/MOL.rst -o min.out -x min.crd -r min.rst\n")
        os.chmod(path, 0o755)

    def _write_scr_min_sp(self, angle_dir: str) -> None:
        path = os.path.join(angle_dir, "scr_min.sh")
        with open(path, "w") as f:
            f.write(f"#!/usr/bin/env bash\nset -e\n{self.sander_cmd} -O -p fragment.top -i sp.in -c fragment.rst -o sp.out -x sp.crd -r sp.rst\n")
        os.chmod(path, 0o755)

    def _write_scr_input(self, gdir: str, angles: List[int]) -> None:
        ang_list = " ".join(str(int(a)) for a in angles)
        path = os.path.join(gdir, "scr_input.sh")
        with open(path, "w") as f:
            f.write(f"""#!/usr/bin/env bash
set -e

ANGLES=({ang_list})

for a in "${{ANGLES[@]}}"; do
  sed "s/XXX/$a/g" restraints.template > "restraints-$a"
  sed "s/XXX/$a/g" min.in.template    > "min.in-$a"
done

a0="${{ANGLES[0]}}"
{self.sander_cmd} -O -p leap/MOL.top -i "min.in-$a0" -c min.rst -o "min-$a0.out" -x "min-$a0.crd" -r "min-$a0.rst"

for ((idx=1; idx<${{#ANGLES[@]}}; idx++)); do
  a="${{ANGLES[$idx]}}"
  prev="${{ANGLES[$((idx-1))]}}"
  {self.sander_cmd} -O -p leap/MOL.top -i "min.in-$a" -c "min-$prev.rst" -o "min-$a.out" -x "min-$a.crd" -r "min-$a.rst"
done
""")
        os.chmod(path, 0o755)

    def _write_scr_en(self, gdir: str, angles: List[int]) -> None:
        ang_list = " ".join(str(int(a)) for a in angles)
        path = os.path.join(gdir, "scr_en.sh")
        with open(path, "w") as f:
            f.write(_SCR_EN_TEMPLATE.format(ang_list=ang_list))
        os.chmod(path, 0o755)

    def _pdb_to_rst_with_cpptraj(self, top_path: str, pdb_path: str, rst_out: str, cwd: str) -> None:
        inp_content = f"parm {top_path}\ntrajin {pdb_path}\ntrajout {rst_out} restart\nrun\nquit\n"
        tmp_in = os.path.join(cwd, "pdb_to_rst.in")
        with open(tmp_in, "w") as f:
            f.write(inp_content)
        subprocess.run(["cpptraj", "-i", tmp_in], check=True, cwd=cwd,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    def _merge_mcs_angle_energies_sp(self, angle_root: str, label: str) -> None:
        angles, energies = [], []
        for item in sorted(os.listdir(angle_root), key=lambda x: int(x) if x.isdigit() else 10 ** 9):
            a_dir = os.path.join(angle_root, item)
            if not os.path.isdir(a_dir) or not item.isdigit():
                continue
            out_file = os.path.join(a_dir, "sp.out")
            if not os.path.exists(out_file):
                self.log.warning("[%s][MCS-SP] missing %s", label, out_file)
                continue
            try:
                e = extract_sp_energy_from_nstep(out_file)
                angles.append(float(item))
                energies.append(float(e))
            except Exception as ex:
                self.log.warning("[%s][MCS-SP] cannot parse ENERGY for angle %s: %s", label, item, ex)

        if not angles:
            raise RuntimeError(f"No MM energies computed in {label} MCS single-point mode")

        dd_ee_path = os.path.join(angle_root, "dd_ee")
        arr = np.column_stack([np.array(angles), np.array(energies)])
        arr = arr[arr[:, 0].argsort()]
        np.savetxt(dd_ee_path, arr, fmt=["%g", "%.10f"])
        write_shifted_profiles(dd_ee_path)

    def _resolve_angles(self, input_pdb, d_tup, step_size, start_from_nearest) -> List[int]:
        if start_from_nearest:
            try:
                atoms = read(input_pdb)
                a0 = nearest_multiple(float(atoms.get_dihedral(*d_tup)), step_size)
                return angles_from_angle0(a0, step_size)
            except Exception as e:
                self.log.warning("[%s] start_from_nearest failed (%s); fallback to 0", type(self).__name__, e)
        return angles_0_350(step_size)

    def _run_cmd(self, cmd, cwd: Optional[str] = None, step_name: str = "") -> None:
        res = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if res.returncode != 0:
            self.log.error("[AMBER] %s failed: %s", step_name, res.stderr.decode(errors="ignore"))
            raise RuntimeError(f"Command failed: {step_name}")

    def _run_mcs_sp_loop(
        self,
        angle_root: str,
        angles: List[int],
        top_path: str,
        label: str,
    ) -> None:
        """Run single-point sander for each angle folder."""
        self._write_sp_in(angle_root)
        for ang in angles:
            a_dir = os.path.join(angle_root, str(ang))
            if not os.path.isdir(a_dir):
                continue
            pdb_path = os.path.join(a_dir, "angle.pdb")
            if not os.path.exists(pdb_path):
                self.log.warning("[%s][MCS-SP] missing %s", label, pdb_path)
                continue
            strip_conect_records(pdb_path)
            try:
                shutil.copy2(top_path, os.path.join(a_dir, "fragment.top"))
                shutil.copy2(os.path.join(angle_root, "sp.in"), os.path.join(a_dir, "sp.in"))
                self._pdb_to_rst_with_cpptraj(
                    top_path=os.path.join(a_dir, "fragment.top"),
                    pdb_path="angle.pdb",
                    rst_out="fragment.rst",
                    cwd=a_dir,
                )
                self._write_scr_min_sp(a_dir)
                self._run_cmd(["bash", "scr_min.sh"], cwd=a_dir, step_name=f"sp-{ang}")
            except Exception as e:
                self.log.error("[%s][MCS-SP][%s] SP failed: %s", label, ang, e)

    def postscan_expand_old_and_releap(
        self,
        method: str,
        d_tup: Tuple[int, int, int, int],
        basis_periodicities: Tuple[int, ...],
        expanded_tag: str = "expanded",
    ) -> Tuple[str, str]:
        d_str = "_".join(str(x) for x in d_tup)
        scan_root = os.path.join(self.cfg.base_dir, "scanning", d_str)
        gaff2_old_method = os.path.join(scan_root, "GAFF2", "old", method)

        ante_src = os.path.join(gaff2_old_method, "antechamber")
        ante_exp = os.path.join(gaff2_old_method, f"antechamber_{expanded_tag}")
        leap_exp = os.path.join(gaff2_old_method, f"leap_{expanded_tag}")
        os.makedirs(ante_exp, exist_ok=True)
        os.makedirs(leap_exp, exist_ok=True)

        mol2_dst = os.path.join(ante_exp, "MOL.mol2")
        frcmod_dst = os.path.join(ante_exp, "MOL.frcmod")
        shutil.copy2(os.path.join(ante_src, "MOL.mol2"), mol2_dst)
        shutil.copy2(os.path.join(ante_src, "MOL.frcmod"), frcmod_dst)

        expand_dihedral_in_frcmod(
            frcmod_path=frcmod_dst,
            mol2_path=mol2_dst,
            d_tup=d_tup,
            basis_periodicities=basis_periodicities,
            frcmod_out=frcmod_dst,
            log=self.log,
        )

        leap_in = (
            f"source leaprc.gaff2\n"
            f"loadamberparams MOL.frcmod\n"
            f"MOL = loadmol2 MOL.mol2\n"
            f"saveamberparm MOL {os.path.join(leap_exp, 'MOL.top')} {os.path.join(leap_exp, 'MOL.rst')}\n"
            f"quit\n"
        )
        with open(os.path.join(ante_exp, "leap.in"), "w") as f:
            f.write(leap_in)
        self._run_cmd([self.tleap_cmd, "-f", "leap.in"], cwd=ante_exp, step_name="tleap-postscan-expanded")

        return os.path.join(leap_exp, "MOL.top"), frcmod_dst


# ============================================================================
# Gaff2SingleDihedralRefinement   (GAFF2 old)
# ============================================================================

@dataclass
class Gaff2SingleDihedralRefinement(_AmberRefinementBase):

    def refine_old(
        self,
        root: str,
        method: str,
        d_tup: Tuple[int, int, int, int],
        basis_periodicities: Tuple[int, ...],
        input_pdb: str,
        net_charge: int = 0,
        spin: int = 1,
        residue_name: str = "MOL",
        MCS: bool = False,
        mcs_xyz: Optional[str] = None,
        postscan_make_expanded_for_mdgx: bool = True,
        step_size: int = 10,
        start_from_nearest: bool = False,
    ) -> None:
        step_size = _validate_step_size(step_size)
        d_str = "_".join(str(x) for x in d_tup)
        scan_root = os.path.join(self.cfg.base_dir, "scanning", d_str)
        gaff2_old_method = os.path.join(scan_root, "GAFF2", "old", method)

        ante_dir = os.path.join(gaff2_old_method, "antechamber")
        leap_dir = os.path.join(gaff2_old_method, "leap")
        os.makedirs(ante_dir, exist_ok=True)
        os.makedirs(leap_dir, exist_ok=True)

        staged_pdb = os.path.join(ante_dir, "input.pdb")
        shutil.copy2(input_pdb, staged_pdb)
        strip_conect_records(staged_pdb)

        ante_args = [
            "antechamber", "-i", "input.pdb", "-fi", "pdb",
            "-o", "MOL.mol2", "-fo", "mol2",
            "-c", "bcc", "-nc", str(net_charge), "-m", str(spin),
            "-rn", residue_name, "-at", "gaff2",
        ]
        self._run_cmd(ante_args, cwd=ante_dir, step_name="antechamber")

        if not MCS:
            self._run_cmd(
                ["antechamber", "-i", "sqm.pdb", "-fi", "pdb",
                 "-o", "MOL.mol2", "-fo", "mol2",
                 "-c", "bcc", "-nc", str(net_charge), "-m", str(spin),
                 "-rn", residue_name, "-at", "gaff2"],
                cwd=ante_dir, step_name="antechamber_sqm",
            )

        self._run_cmd(
            ["parmchk2", "-i", "MOL.mol2", "-f", "mol2", "-o", "MOL.frcmod", "-s", "2", "-a", "Y"],
            cwd=ante_dir, step_name="parmchk2",
        )

        leap_in = "source leaprc.gaff2\nloadamberparams MOL.frcmod\nMOL = loadmol2 MOL.mol2\nsaveamberparm MOL ../leap/MOL.top ../leap/MOL.rst\nquit\n"
        with open(os.path.join(ante_dir, "leap.in"), "w") as f:
            f.write(leap_in)
        self._run_cmd([self.tleap_cmd, "-f", "leap.in"], cwd=ante_dir, step_name="tleap")

        top_path = os.path.join(leap_dir, "MOL.top")
        rst0_path = os.path.join(leap_dir, "MOL.rst")
        if not (os.path.isfile(top_path) and os.path.isfile(rst0_path)):
            raise FileNotFoundError(f"Missing MOL.top/MOL.rst in {leap_dir}")

        angles = self._resolve_angles(input_pdb, d_tup, step_size, start_from_nearest)

        if not MCS:
            self._write_min_template(gaff2_old_method)
            self._write_restrain_template(gaff2_old_method, d_tup)
            self._write_min0_file(gaff2_old_method)
            self._write_scr_min(gaff2_old_method)
            self._run_cmd(["bash", "scr_min.sh"], cwd=gaff2_old_method, step_name="scr_min")
            self._write_scr_input(gaff2_old_method, angles=angles)
            self._write_scr_en(gaff2_old_method, angles=angles)
            self._run_cmd(["bash", "scr_input.sh"], cwd=gaff2_old_method, step_name="scr_input")
            self._run_cmd(["bash", "scr_en.sh"], cwd=gaff2_old_method, step_name="scr_en")
            if postscan_make_expanded_for_mdgx:
                self.postscan_expand_old_and_releap(method, d_tup, basis_periodicities)
            return

        # MCS path
        if not mcs_xyz:
            raise ValueError("MCS=True requires mcs_xyz=<path to MCS/geometries.xyz>")
        cleanup_mcs_angle_outputs(gaff2_old_method, step=step_size, log=self.log)
        write_mcs_frames_to_angle_folders(mcs_xyz, gaff2_old_method, angles=angles)
        self._run_mcs_sp_loop(gaff2_old_method, angles, top_path, label="GAFF2/old")
        self._merge_mcs_angle_energies_sp(gaff2_old_method, label="GAFF2/old")


# ============================================================================
# Gaff2NewDihedralRefinement   (GAFF2 new, fitted params)
# ============================================================================

@dataclass
class Gaff2NewDihedralRefinement(_AmberRefinementBase):

    def refine_new(
        self,
        root: str,
        method: str,
        d_tup: Tuple[int, int, int, int],
        input_pdb: str,
        param_results: Dict[float, Dict[str, str]],
        net_charge: int = 0,
        spin: int = 1,
        residue_name: str = "fragment",
        MCS: bool = False,
        mcs_xyz: Optional[str] = None,
        step_size: int = 10,
        start_from_nearest: bool = False,
    ) -> None:
        step_size = _validate_step_size(step_size)
        d_str = "_".join(str(x) for x in d_tup)
        scan_root = os.path.join(self.cfg.base_dir, "scanning", d_str)
        gaff2_new_method = os.path.join(scan_root, "GAFF2", "new", method)
        os.makedirs(gaff2_new_method, exist_ok=True)

        mol2_src = os.path.join(scan_root, "GAFF2", "old", method, "antechamber", "MOL.mol2")
        if not os.path.isfile(mol2_src):
            raise FileNotFoundError(f"[GAFF2/new][{method}][{d_str}] Missing MOL.mol2 at {mol2_src}")
        if MCS and not mcs_xyz:
            raise ValueError("MCS=True requires mcs_xyz=<path to MCS/geometries.xyz>")

        angles = self._resolve_angles(input_pdb, d_tup, step_size, start_from_nearest)

        for hrst_val, paths in sorted(param_results.items()):
            frcmod_src = paths.get("frcmod")
            if not frcmod_src or not os.path.isfile(frcmod_src):
                self.log.warning("[GAFF2/new][%s][%s] hrst=%s: missing frcmod -> skip", method, d_tup, hrst_val)
                continue

            hrst_str = str(hrst_val)
            hrst_dir = os.path.join(gaff2_new_method, hrst_str)
            ante_dir = os.path.join(hrst_dir, "antechamber")
            leap_dir = os.path.join(hrst_dir, "leap")
            os.makedirs(ante_dir, exist_ok=True)
            os.makedirs(leap_dir, exist_ok=True)

            mol2_dst = os.path.join(ante_dir, "MOL.mol2")
            frcmod_dst = os.path.join(ante_dir, "MOL.frcmod")
            shutil.copy2(mol2_src, mol2_dst)
            shutil.copy2(frcmod_src, frcmod_dst)

            leap_in = (
                f"source leaprc.gaff2\n"
                f"loadamberparams MOL.frcmod\n"
                f"MOL = loadmol2 MOL.mol2\n"
                f"saveamberparm MOL ../leap/MOL.top ../leap/MOL.rst\n"
                f"quit\n"
            )
            with open(os.path.join(ante_dir, "leap.in"), "w") as f:
                f.write(leap_in)
            self._run_cmd([self.tleap_cmd, "-f", "leap.in"], cwd=ante_dir, step_name=f"tleap-{hrst_str}")

            top_path = os.path.join(leap_dir, "MOL.top")
            rst_path = os.path.join(leap_dir, "MOL.rst")
            if not (os.path.isfile(top_path) and os.path.isfile(rst_path)):
                self.log.error("[GAFF2/new][%s][%s] hrst=%s: tleap failed", method, d_tup, hrst_str)
                continue

            if not MCS:
                self._write_min_template(hrst_dir)
                self._write_restrain_template(hrst_dir, d_tup)
                self._write_min0_file(hrst_dir)
                self._write_scr_min(hrst_dir)
                self._run_cmd(["bash", "scr_min.sh"], cwd=hrst_dir, step_name=f"scr_min-{hrst_str}")
                self._write_scr_input(hrst_dir, angles=angles)
                self._write_scr_en(hrst_dir, angles=angles)
                self._run_cmd(["bash", "scr_input.sh"], cwd=hrst_dir, step_name=f"scr_input-{hrst_str}")
                self._run_cmd(["bash", "scr_en.sh"], cwd=hrst_dir, step_name=f"scr_en-{hrst_str}")
                continue

            # MCS path
            cleanup_mcs_angle_outputs(hrst_dir, step=step_size, log=self.log)
            write_mcs_frames_to_angle_folders(mcs_xyz, hrst_dir, angles=angles)
            self._run_mcs_sp_loop(hrst_dir, angles, top_path, label=f"GAFF2/new-{hrst_str}")
            self._merge_mcs_angle_energies_sp(hrst_dir, label=f"GAFF2/new-{hrst_str}")

