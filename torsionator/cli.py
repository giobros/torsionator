import argparse, os
from .config import Config
from .workflow import Workflow
from .conect_fix import ensure_conect_with_obabel

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pdb", required=True)
    p.add_argument("--method", choices=["obi","mace","all"], default="obi")
    p.add_argument("--dihedral", type=str, default="all",
                  help="Use 'all' or a dihedral as [a,b,c,d], or 'print' to list them.")
    p.add_argument("--conf_analysis", choices=["true","false"], default="true",
                     help = "false → scan everything anyway (clashes ignored), " \
                     "true → generate conformers and check for lower energy clash-free one")
    p.add_argument("--force_scanning", choices=["true","false"], default="true",
                    help = "if no clash-free conformer is found" \
                    " false → stop the workflow, " \
                     "true → uses best-conformer for scanning (lower LJ)")
    args = p.parse_args()
    #use conect_fix to ensure CONECT records are present and correct
    ensure_conect_with_obabel(args.pdb, in_place=True)  # if no CONECT, uses Open Babel (-xp) and overwrites the PDB

    cfg = Config(base_dir="/data")
    wf = Workflow(cfg)
    wf.run(args.pdb, args.method, args.dihedral, args.conf_analysis, args.force_scanning)

if __name__ == "__main__":
    main()
