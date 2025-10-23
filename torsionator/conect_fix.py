import os, shutil, subprocess, tempfile

def strip_generic_header(pdb_path: str, in_place: bool = True, *, start_tags=("ATOM  ", "HETATM")) -> bool:
    """
    Remove any leading lines before the first structural record.
    By default, the first record must be 'ATOM  ' or 'HETATM' (exact PDB columns 1–6).
    """
    with open(pdb_path, "r") as f:
        lines = f.readlines()

    # find first true ATOM/HETATM record (column-accurate)
    first_idx = None
    for i, ln in enumerate(lines):
        tag = ln[:6]  # PDB record name field
        if tag in start_tags:
            first_idx = i
            break

    # nothing to trim or no atoms found
    if first_idx in (None, 0):
        return False

    trimmed = lines[first_idx:]
    out_path = pdb_path if in_place else os.path.splitext(pdb_path)[0] + "_trimmed.pdb"

    tmp = out_path + ".tmp"
    with open(tmp, "w") as w:
        w.writelines(trimmed)
        if not trimmed or not trimmed[-1].endswith("\n"):
            w.write("\n")
    os.replace(tmp, out_path)
    return True

def pdb_has_conect(pdb_path: str) -> bool:
    with open(pdb_path, "r") as f:
        for ln in f:
            if ln.startswith("CONECT"):
                return True
    return False

def ensure_conect_with_obabel(pdb_path: str, sdf_template: str = None, in_place: bool = True) -> bool:
  
    if pdb_has_conect(pdb_path):
        return False  
    strip_generic_header(pdb_path, in_place=in_place)
    with tempfile.TemporaryDirectory() as td:
        cmd = ["obabel", "-ipdb", pdb_path, "-opdb", "-O", pdb_path, "-xp"]
        r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"OpenBabel failed: {r.stderr}")
    
    return True