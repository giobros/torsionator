#!/bin/bash
set -euo pipefail

# ---------------- USER SETTINGS ----------------

PDB_FILE_NAME_ROOT="NAME"        # CHANGE THIS! file NAME of your PDB file, without .pdb
PDB_FILE_DIR="your_folder_path"  # CHANGE THIS! folder (full path) that contains $PDB_FILE_NAME_ROOT.pdb
METHOD="all"                     # "all" | "mace" | "obi", NN calculator to use (default: obi)
DIHEDRAL="all"                   # "all" | "[a,b,c,d]" | "print" (0-based indices): "all" to scan all rotatable bonds; "print" to list them; "[a,b,c,d]" for a specific one.
CONF_ANALYSIS="false"            # "true" | "false" | "none", "false" → scan minimized input even if clashes exist; "true"  → generate conformers and use a clash-free starting geometry
BCS="false"                      # "true" | "false" | "none", "false" → abort; "true" → use the conformer with lowest LJ energy.
MCS="true"                       # "true" | "false" | "none", "true"→  find lower-energy conformations per angle 
N_CONF=20                        # number of conformers
RMSD=0.5                         # RMSD pruning threshold
MULTIPLICITY=6                   # max expantion multiplicity (0 to keep the GAFF2 original one)
STEP_SIZE=10                     # scan steps (5,10,15,20)
DOUBLE_ROTATION="true"           # "true" | "false" | "none", false" →  just clockwise (cw), "true" →  both clockwise (cw) counterclockwise (ccw) scan when MCS=true; "
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
