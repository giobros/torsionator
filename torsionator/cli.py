import argparse
import os

from .config import Config
from .workflow import Workflow
from .conect_fix import ensure_conect_with_obabel
from .io_utils import backup_existing_outputs


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="torsionator: torsion-parameter refinement with NN potentials + GAFF2."
    )
    p.add_argument("--pdb", required=True, help="Input PDB file.")
    p.add_argument(
        "--method", choices=["obi", "mace", "uma"], default="obi",
        help="NN calculator to use (default: obi).",
    )
    p.add_argument(
        "--dihedral", type=str, default="all",
        help="'all' to scan all rotatable bonds, 'print' to list them, or '[a,b,c,d]' for a specific one.",
    )
    p.add_argument(
        "--conf_analysis", choices=["true", "false"], default="false",
        help=(
            "false → scan minimized input even if clashes exist; "
            "true  → generate conformers and use a clash-free starting geometry."
        ),
    )
    p.add_argument(
        "--BCS", choices=["true", "false"], default="false",
        help=(
            "Best-Conformer-Scan fallback: "
            "if no clash-free conformer is found, "
            "false → abort; true → use the conformer with lowest LJ energy."
        ),
    )
    p.add_argument(
        "--multiplicity", type=int, choices=[0, 1, 2, 3, 4, 6], default=6,
        help="Max torsion periodicity preset. 0 = keep frcmod multiplicities unchanged.",
    )
    p.add_argument("--RMSD", type=float, default=0.5, help="RMSD pruning threshold for conformer screening.")
    p.add_argument("--n_confs", type=int, default=20, help="Number of conformers to generate.")
    p.add_argument(
        "--MCS", choices=["true", "false"], default="true",
        help="Multi-Conformer Scan: use minimum-energy collapse per angle across all conformers.",
    )
    p.add_argument(
        "--double_rotation", choices=["true", "false"], default="true",
        help="(MCS only) Also scan counter-clockwise in addition to clockwise.",
    )
    p.add_argument(
        "--step_size", type=int, default=10,
        help="Dihedral scan step in degrees (must divide 360). Example: 5, 10, 15, 20.",
    )
    p.add_argument(
        "--net_charge", type=int, default=0,
        help="Net charge of the molecule for antechamber (default: 0).",
    )
    p.add_argument(
        "--spin", type=int, default=1,
        help="Spin multiplicity of the molecule (default: 1). Only UMA supports open-shell NNP calculations.",
    )
    return p


def main() -> None:
    args = _build_parser().parse_args()

    ensure_conect_with_obabel(args.pdb, in_place=True)

    cfg = Config(base_dir="/data")
    backup_existing_outputs(cfg.base_dir)
    wf = Workflow(cfg)
    wf.run(
        pdb_file=args.pdb,
        method=args.method,
        dihedral=args.dihedral,
        conf_analysis=args.conf_analysis,
        multiplicity=args.multiplicity,
        BCS=args.BCS,
        RMSD=args.RMSD,
        n_confs=args.n_confs,
        MCS=args.MCS,
        double_rotation=args.double_rotation,
        step_size=args.step_size,
        net_charge=args.net_charge,
        spin=args.spin,
    )


if __name__ == "__main__":
    main()
