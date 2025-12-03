#!/bin/bash
set -euo pipefail

unset PYTHONPATH

# ---------------- USER SETTINGS ----------------

ROOT="ROOT"                     # CHANGE THIS! name of your PDB, without .pdb
DATA_DIR="$HOME/your_folder"    # CHANGE THIS! folder that contains $ROOT.pdb

METHOD="all"                     # "all" | "mace" | "obi"
DIHEDRAL="all"                   # "all" | "[a,b,c,d]" | "print"  (0-based indices)
CONF_ANALYSIS="false"            # "true" | "false" | "none"
BCS="false"                      # "true" | "false" | "none"
MCS="false"                      # "true" | "false" | "none"
N_CONF=50                        # number of conformers
RMSD=0.5                         # RMSD pruning threshold

# ------------------------------------------------

apptainer exec \
  --nv \
  --bind "$HOME/torsionator/torsionator:/torsionator" \
  --bind "$DATA_DIR:/data" \
  --env PYTHONPATH=/ \
  --env HOME=/root \
  torsionator.sif \
  python3.9 -m torsionator.cli \
    --pdb "/data/${ROOT}.pdb" \
    --method "$METHOD" \
    --dihedral "$DIHEDRAL" \
    --conf_analysis "$CONF_ANALYSIS" \
    --BCS "$BCS" \
    --MCS "$MCS" \
    --N "$N_CONF" \
    --RMSD "$RMSD"
