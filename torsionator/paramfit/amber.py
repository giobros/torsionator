# paramfit/amber.py

import os
import re
import numpy as np
from typing import Dict, List, Tuple, Optional

__all__ = [
    "get_atom_types_from_mol2",
    "param_input_file",
    "update_frcmod",
]

# ------------------------------
# AMBER / MDGX input preparation
# ------------------------------

def get_atom_types_from_mol2(mol2_file: str) -> Dict[int, str]:
    """
    Read Tripos MOL2 and return 0-based atom index -> atom-type (e.g., 'c3', 'os')
    """
    atom_types: Dict[int, str] = {}
    with open(mol2_file) as f:
        lines = f.readlines()
    atom_section = False
    for line in lines:
        if line.startswith("@<TRIPOS>ATOM"):
            atom_section = True
            continue
        if line.startswith("@<TRIPOS>BOND"):
            break
        if atom_section:
            parts = line.split()
            if len(parts) >= 6:
                atom_types[int(parts[0]) - 1] = parts[5]
    return atom_types


def param_input_file(
    atom_types: Dict[int, str],
    method: str,
    dihedral_indices: Tuple[int, int, int, int],
    scan_dir: str,
    method_label: str,
) -> None:
    """
    Create an mdgx parameterization input for a single dihedral scan.
    Writes: <scan_dir>/<method_label>/parametrization_<label>_dihedral_<a_b_c_d>.in
    """
    dih_types = " ".join([atom_types[i] for i in dihedral_indices])
    leap_top   = os.path.abspath(os.path.join(scan_dir, f"leap_{method_label}", "MOL.top"))
    frcmod     = os.path.abspath(os.path.join(scan_dir, f"antechamber_{method_label}", "MOL.frcmod"))
    out_dat    = os.path.abspath(os.path.join(scan_dir, method_label, f"{method.lower()}.dat"))
    out_out    = os.path.abspath(os.path.join(scan_dir, method_label, f"{method.lower()}.out"))
    traj       = os.path.abspath(os.path.join(scan_dir, method_label, "geometries.nc"))
    energy_file= os.path.abspath(os.path.join(scan_dir, method_label, "energies.dat"))

    # Load to ensure energies.dat exists & is readable (mdgx references it)
    _ = np.loadtxt(energy_file)  # noqa: F841

    pfile = os.path.join(
        scan_dir,
        method_label,
        f"parametrization_{method_label}_dihedral_{'_'.join(map(str, dihedral_indices))}.in",
    )
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
        f.write("  hrst       0.0002,\n&end\n")


# ------------------------------
# frcmod updating helpers
# ------------------------------

def _key_from_frcmod(line: str) -> Optional[str]:
    """
    Extract a 2-char padded DIHE key from an frcmod line: 'AA-BB-CC-DD'
    """
    m = re.match(r'^(.{2})-(.{2})-(.{2})-(.{2})', line)
    if m:
        return f"{m.group(1)}-{m.group(2)}-{m.group(3)}-{m.group(4)}"
    return None


def _key_from_dat(line: str) -> Optional[str]:
    """
    Extract a 2-char padded key from an mdgx .dat torsion line.
    Accepts tokens of length 1–2 and pads with spaces on the right to width 2.
    """
    s = line.strip()
    m = re.match(
        r'^([A-Za-z0-9]{1,2})\s*-\s*([A-Za-z0-9]{1,2})\s*-\s*([A-Za-z0-9]{1,2})\s*-\s*([A-Za-z0-9]{1,2})',
        s,
    )
    if not m:
        return None
    atoms = [f"{m.group(i):<2}" for i in range(1, 5)]
    return f"{atoms[0]}-{atoms[1]}-{atoms[2]}-{atoms[3]}"


def _parse_dat_line(line: str) -> Tuple[Optional[str], Optional[Dict[str, float]]]:
    """
    Parse one mdgx torsion line and return:
      key: 'AA-BB-CC-DD'
      params: dict(divisor, barrier, phase, periodicity)
    """
    key = _key_from_dat(line)
    if not key:
        return None, None

    parts = line.strip().split()
    p0 = None
    for i, p in enumerate(parts):
        try:
            float(p)
            p0 = i
            break
        except ValueError:
            pass
    if p0 is None or len(parts) < p0 + 4:
        return None, None

    return key, {
        "divisor": float(parts[p0]),
        "barrier": float(parts[p0 + 1]),
        "phase": float(parts[p0 + 2]),
        "periodicity": float(parts[p0 + 3]),
    }


def _format_frcmod_line(key: str, params: Dict[str, float]) -> str:
    """
    Format an frcmod DIHE line using your original widths (spacing preserved).
    """
    key = f"{key:<11}"
    return (
        f"{key}"
        f"{params['divisor']:>4.0f}"
        f"{params['barrier']:>9.3f}"
        f"{params['phase']:>14.3f}"
        f"{params['periodicity']:>16.3f}"
    )


def _reverse_key(key: str) -> str:
    parts = key.split("-")
    return "-".join(reversed(parts))


# ------------------------------
# Main merge function
# ------------------------------

def update_frcmod(
    frcmod_path: str,
    dat_path: str,
    method: str,
    output_folder: str,
) -> str:
    """
    Authoritative multi-term merge for torsion keys present in mdgx .dat:

    """
    # --- 1) Parse mdgx .dat: collect all torsion lines, grouped by key (keep order) ---
    grouped: Dict[str, List[Dict[str, float]]] = {}
    key_order: List[str] = []

    with open(dat_path, "r") as f:
        dat_lines = f.readlines()

    in_dihe = False
    for line in dat_lines:
        s = line.strip()
        # accept either explicit "DIHE" section or free torsion lines
        if s.upper() == "DIHE":
            in_dihe = True
            continue
        if in_dihe and (
            s == "" or s.upper().startswith(("IMPROPER", "IMPR", "NONB", "MASS", "BOND", "ANGLE", "HBON"))
        ):
            in_dihe = False
            continue

        k, p = _parse_dat_line(line)
        if k and p:
            if k not in grouped:
                grouped[k] = []
                key_order.append(k)
            grouped[k].append(p)

    # Decide output mode
    os.makedirs(output_folder, exist_ok=True)
    inplace = (os.path.dirname(os.path.abspath(frcmod_path)) == os.path.abspath(output_folder))
    out_file = frcmod_path if inplace else os.path.join(output_folder, f"new_{method}.frcmod")

    # If nothing to merge, just copy/keep as expected
    if not grouped:
        if not inplace:
            with open(frcmod_path, "r") as fin, open(out_file, "w") as fout:
                fout.write(fin.read())
        return out_file

    mdgx_keys = set(key_order)
    mdgx_keys_rev = {_reverse_key(k) for k in mdgx_keys}

    # --- 2) Read seed frcmod and locate DIHE section ---
    with open(frcmod_path, "r") as f:
        orig = f.readlines()

    dihe_start = None
    dihe_end = None
    for i, line in enumerate(orig):
        if line.strip().upper() == "DIHE":
            dihe_start = i
            j = i + 1
            while j < len(orig):
                t = orig[j].strip().upper()
                if t in ("IMPROPER", "IMPR", "NONB", "MASS", "BOND", "ANGLE", "HBON"):
                    break
                j += 1
            dihe_end = j
            break

    def _format_group(key: str) -> List[str]:
        return [_format_frcmod_line(key, p) + "\n" for p in grouped[key]]

    # --- 3) Build new contents ---
    if dihe_start is None:
        # No DIHE section in seed → append a new one
        new_contents = orig[:]
        if len(new_contents) and not new_contents[-1].endswith("\n"):
            new_contents[-1] = new_contents[-1] + "\n"
        new_contents.append("DIHE\n")
        for k in key_order:
            new_contents.extend(_format_group(k))
    else:
        pre  = orig[: dihe_start + 1]         # includes "DIHE"
        body = orig[dihe_start + 1 : dihe_end]
        post = orig[dihe_end :]

        # Drop any existing entries for these keys (forward or reverse)
        kept: List[str] = []
        for ln in body:
            key = _key_from_frcmod(ln)
            if not key:
                kept.append(ln)  # comments/blank
                continue
            if key in mdgx_keys or key in mdgx_keys_rev:
                continue         # drop old entries for these keys
            kept.append(ln)

        # Append authoritative mdgx lines
        # Drop trailing empty/comment-only lines before appending
        while kept and kept[-1].strip() == "":
            kept.pop()

        # Now append new authoritative lines
        new_dihe = kept[:]
        for k in key_order:
            new_dihe.extend(_format_group(k))
      

        new_contents = pre + new_dihe + post

    # --- 4) Write (in place for cumulative; temp file in method_dir for per-dihedral) ---
    with open(out_file, "w") as f:
        f.writelines(new_contents)

    return out_file
