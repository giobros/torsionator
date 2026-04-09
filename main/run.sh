#!/bin/bash
set -euo pipefail

# ---------------- USER SETTINGS ----------------

ROOT="fragment"                     # CHANGE THIS! name of your PDB, without .pdb
DATA_DIR="../tests/16_mcs_5"        # CHANGE THIS! folder that contains $ROOT.pdb

METHOD="mace"                       # "mace" | "obi" | "uma"
DIHEDRAL="[0,1,2,3]"               # "all" | "[a,b,c,d]" | "print"  (0-based indices)
CONF_ANALYSIS="false"               # "true" | "false" | "none"
BCS="false"                         # "true" | "false" | "none"
MCS="true"                          # "true" | "false" | "none"
N_CONF=1                            # number of conformers
RMSD=0.5                            # RMSD pruning threshold
MULTIPLICITY=6                      # max expansion multiplicity (0 to keep the GAFF2 original one)
STEP_SIZE=5                         # scan steps (5,10,15,20)
DOUBLE_ROTATION="true"              # "true" | "false" | "none"
NET_CHARGE=0
# ------------------------------------------------

# Resolve the torsionator package root (one level up from this script's directory)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TORSIONATOR_ROOT="$(dirname "$SCRIPT_DIR")"

# shellcheck source=/dev/null
source "$(conda info --base)/etc/profile.d/conda.sh"

case "$METHOD" in
  mace)
    conda activate torsionator_mace
    export PYTHONPATH="$TORSIONATOR_ROOT:${PYTHONPATH:-}"
    ;;
  obi)
    conda activate torsionator_obi
    OBIWAN_DIR="$HOME/OBIWAN"
    if [ ! -d "$OBIWAN_DIR" ]; then
        echo "[setup] Cloning OBIWAN from GitHub..."
        git clone https://github.com/virtualmartire/OBIWAN.git "$OBIWAN_DIR"
    fi
    ASE_CALC_DIR="$(python3 -c 'import ase; import os; print(os.path.join(os.path.dirname(ase.__file__), "calculators"))')"
    cp -f "$SCRIPT_DIR/obi.py" "$ASE_CALC_DIR/obi.py"
    export OBIWAN_MODEL_PATH="$OBIWAN_DIR/results/models"
    export PYTHONPATH="$TORSIONATOR_ROOT:$OBIWAN_DIR:${PYTHONPATH:-}"
    ;;
  uma)
    conda activate torsionator_uma
    export PYTHONPATH="$TORSIONATOR_ROOT:${PYTHONPATH:-}"
    export UMA_MODEL_PATH="$TORSIONATOR_ROOT/torsionator/uma-s-1p1.pt"
    ;;
  *)
    echo "ERROR: unknown method '$METHOD'. Choose mace, obi, or uma." >&2
    exit 1
    ;;
esac

python3 -m torsionator.${METHOD}_torsionator \
    --pdb "${DATA_DIR}/${ROOT}.pdb" \
    --dihedral "$DIHEDRAL" \
    --conf_analysis "$CONF_ANALYSIS" \
    --BCS "$BCS" \
    --MCS "$MCS" \
    --n_confs "$N_CONF" \
    --RMSD "$RMSD" \
    --multiplicity "$MULTIPLICITY" \
    --step_size "$STEP_SIZE" \
    --double_rotation "$DOUBLE_ROTATION" \
    --net_charge "$NET_CHARGE"
