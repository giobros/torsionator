import os
import subprocess

# Supported input formats; "pdb" needs no conversion, the rest are Open
# Babel format codes for non-PDB inputs. Extend as needed.
INPUT_FORMATS = ("pdb", "sdf", "mol2", "mol", "xyz", "cif", "gro", "pdbqt")


def add_input_format_args(parser) -> None:
    """Add one mutually-exclusive CLI flag per entry in INPUT_FORMATS.

    Exactly one of --pdb/--sdf/--mol2/... must be given; whichever is used
    also tells resolve_input_pdb() which Open Babel format code to convert
    from (if not already PDB).
    """
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--pdb", help="Input PDB file.")
    for fmt in INPUT_FORMATS[1:]:
        group.add_argument(
            f"--{fmt}",
            help=f"Input {fmt.upper()} file (converted to PDB internally via Open Babel).",
        )


def convert_to_pdb(input_path: str, fmt: str) -> str:
    """Convert *input_path* (in format *fmt*) to PDB via Open Babel.

    Returns *input_path* unchanged if fmt == "pdb"; otherwise writes a
    sibling file with the same stem and a .pdb extension and returns that
    path.
    """
    if fmt == "pdb":
        return input_path

    out_path = os.path.splitext(input_path)[0] + ".pdb"
    result = subprocess.run(
        ["obabel", f"-i{fmt}", input_path, "-opdb", "-O", out_path, "-xp"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"OpenBabel failed converting {input_path} ({fmt}) to PDB: {result.stderr}")
    return out_path


def resolve_input_pdb(args) -> str:
    """Given parsed CLI args carrying one of the INPUT_FORMATS flags, return a PDB path."""
    for fmt in INPUT_FORMATS:
        path = getattr(args, fmt, None)
        if path is not None:
            return convert_to_pdb(path, fmt)
    raise ValueError(f"No input structure file provided (expected one of --{', --'.join(INPUT_FORMATS)}).")


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
