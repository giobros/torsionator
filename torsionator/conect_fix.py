import os
import subprocess


def _strip_generic_header(pdb_path: str, in_place: bool = True) -> bool:
    """Remove any leading lines before the first ATOM/HETATM record."""
    with open(pdb_path) as f:
        lines = f.readlines()

    first_idx = next(
        (i for i, ln in enumerate(lines) if ln[:6] in ("ATOM  ", "HETATM")),
        None,
    )

    if first_idx in (None, 0):
        return False

    trimmed = lines[first_idx:]
    out_path = pdb_path if in_place else os.path.splitext(pdb_path)[0] + "_trimmed.pdb"

    tmp = out_path + ".tmp"
    with open(tmp, "w") as f:
        f.writelines(trimmed)
        if not trimmed or not trimmed[-1].endswith("\n"):
            f.write("\n")
    os.replace(tmp, out_path)
    return True


def _pdb_has_conect(pdb_path: str) -> bool:
    with open(pdb_path) as f:
        return any(ln.startswith("CONECT") for ln in f)


def ensure_conect_with_obabel(pdb_path: str, in_place: bool = True) -> bool:
    """Add CONECT records via OpenBabel if they are missing. Returns True if modified."""
    if _pdb_has_conect(pdb_path):
        return False

    _strip_generic_header(pdb_path, in_place=in_place)
    result = subprocess.run(
        ["obabel", "-ipdb", pdb_path, "-opdb", "-O", pdb_path, "-xp"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"OpenBabel failed: {result.stderr}")
    return True
