import argparse, os
from .config import Config
from .workflow import Workflow

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pdb", required=True)
    p.add_argument("--method", choices=["obi","mace","obi"], default="obi")
    p.add_argument("--dihedral", type=str, default="all",
                  help="Use 'all' or a dihedral as [a,b,c,d], or 'print' to list them.")
    p.add_argument("--conf_analysis", choices=["true","false"], default="true",
                     help = "false → scan everything anyway (clashes ignored), " \
                     "true → generate conformers and check for lower energy clash-free one")
    p.add_argument("--force_scan", choices=["true","false"], default="true",
                    help = "if no clash-free conformer is found" \
                    " false → stop the workflow, " \
                     "true → uses best-conformer for scanning (lower LJ)")
    args = p.parse_args()

    cfg = Config(base_dir="/data")
    wf = Workflow(cfg)
    wf.run(args.pdb, args.method, args.dihedral, args.conf_analysis, args.force_scan)

if __name__ == "__main__":
    main()
