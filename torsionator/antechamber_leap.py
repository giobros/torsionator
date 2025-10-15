import os
import sys
import shutil
import subprocess
from dataclasses import dataclass



def _prepend_to_path(path: str):
    if path and os.path.isdir(path):
        os.environ["PATH"] = f"{path}:{os.environ.get('PATH', '')}"


def _which_all(cmds):
    missing = [c for c in cmds if shutil.which(c) is None]
    return missing

def run_command(cmd_args, cwd=None):
    try:
        result = subprocess.run(cmd_args, cwd=cwd, check=True)
        return result.returncode
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Command failed ({e.returncode}): {' '.join(cmd_args)}") from e


@dataclass
class AntechamberConfig:
    net_charge: int = 0
    residue_name: str = "MOL"
    forcefield: str = "gaff2"  # leaprc.gaff2
    at_type: str = "gaff2"     # -at gaff2



def run_workflow(scan_dir: str, method: str, cfg: AntechamberConfig = AntechamberConfig()):
    """
    scan_dir: path to the dihedral scan directory (the one that contains antechamber_<METHOD>/)
    method: 'obi' or 'mace' (case-insensitive)
    cfg: options (net charge, residue name, etc.)
    """
    scan_dir = os.path.abspath(scan_dir)
    ante_dir = os.path.join(scan_dir, f"antechamber_{method}")
    leap_dir = os.path.join(scan_dir, f"leap_{method}")

    os.makedirs(ante_dir, exist_ok=True)
    os.makedirs(leap_dir, exist_ok=True)

    input_pdb = os.path.join(ante_dir, "input.pdb")
    if not os.path.exists(input_pdb):
        sys.exit(f"Error: expected input PDB at: {input_pdb}")

    # Ensure AmberTools is on PATH
    _prepend_to_path(os.environ.get("AMBERTOOLS_BIN", ""))
    _prepend_to_path("/opt/conda/bin")

    missing = _which_all(["antechamber", "parmchk2", "tleap"])
    if missing:
        sys.exit(f"AmberTools commands not found in PATH: {', '.join(missing)}")

    # === Run Antechamber ===
    ac_cmd = [
        "antechamber",
        "-i", "input.pdb",
        "-fi", "pdb",
        "-o", "MOL.mol2",
        "-fo", "mol2",
        "-c", "bcc",
        "-nc", str(cfg.net_charge),
        "-rn", cfg.residue_name,
        "-at", cfg.at_type
    ]
    run_command(ac_cmd, cwd=ante_dir)

    if os.path.isfile(os.path.join(ante_dir, "sqm.pdb")):
        ac2_cmd = [
            "antechamber",
            "-i", "sqm.pdb",
            "-fi", "pdb",
            "-o", "MOL.mol2",
            "-fo", "mol2",
            "-c", "bcc",
            "-nc", str(cfg.net_charge),
            "-rn", cfg.residue_name,
            "-at", cfg.at_type
        ]
        run_command(ac2_cmd, cwd=ante_dir)

    # parmchk2 → frcmod
    run_command(
        ["parmchk2", "-i", "MOL.mol2", "-f", "mol2", "-o", "MOL.frcmod", "-s", "2", "-a", "Y"],
        cwd=ante_dir
    )

    # === tleap ===
    leap_in = f"""source leaprc.{cfg.forcefield}
MOL = loadmol2 MOL.mol2
loadamberparams MOL.frcmod
saveamberparm MOL ../leap_{method}/MOL.top ../leap_{method}/MOL.rst
quit
"""
    with open(os.path.join(ante_dir, "leap.in"), "w") as f:
        f.write(leap_in)

    run_command(["tleap", "-f", "leap.in"], cwd=ante_dir)


# ---- CLI / module entry ------------------------------------------------------

def main():
    if len(sys.argv) < 3 or len(sys.argv) > 5:
        print("Usage: python -m torsionator.antechamber_leap <scan_dir> <METHOD> [net_charge] [res_name]")
        sys.exit(1)

    scan_dir = sys.argv[1]
    method = sys.argv[2]
    net_charge = int(sys.argv[3]) if len(sys.argv) >= 4 else 0
    res_name = sys.argv[4] if len(sys.argv) == 5 else "MOL"

    if not os.path.isdir(scan_dir):
        sys.exit(f"Error: directory not found: {scan_dir}")

    cfg = AntechamberConfig(net_charge=net_charge, residue_name=res_name)
    try:
        run_workflow(scan_dir, method, cfg)
    except Exception as e:
        sys.exit(str(e))

if __name__ == "__main__":
    main()
