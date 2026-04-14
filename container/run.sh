#!/bin/bash
set -euo pipefail

# ---------------- USER SETTINGS ----------------

PDB_FILE_NAME_ROOT="NAME"        # CHANGE THIS! file NAME of your PDB file, without .pdb
PDB_FILE_DIR="your_folder_path"  # CHANGE THIS! folder (full path) that contains $PDB_FILE_NAME_ROOT.pdb
METHOD="obi"                     # "all" | "mace" | "obi", NN calculator to use (default: obi)
DIHEDRAL="all"                   # "all" | "[a,b,c,d]" | "print" (0-based indices): "all" to scan all rotatable bonds; "print" to list them; "[a,b,c,d]" for a specific one.
CONF_ANALYSIS="false"            # "true" | "false" | "none", "false" → scan minimized input even if clashes exist; "true"  → generate conformers, clash-free starting geometry
BCS="false"                      # "true" | "false" | "none", "false" → abort; "true" → use the conformer with lowest LJ energy.
MCS="true"                       # "true" | "false" | "none", "true"→  find lower-energy conformations per angle 
N_CONF=20                        # number of conformers
RMSD=0.5                         # RMSD pruning threshold
MULTIPLICITY=6                   # max expantion multiplicity (0 to keep the GAFF2 original one)
STEP_SIZE=10                     # scan steps (5,10,15,20)
DOUBLE_ROTATION="true"           # "true" | "false" | "none", false" →  just clockwise (cw), "true" →  both clockwise (cw) counterclockwise (ccw) scan when MCS=true; "
NET_CHARGE="0"                   # net molecule charge, (default 0, other charges raise warning with mace and obi)
SPIN="1"                         # molecule spin (2S+1) (default 1; other spins raise warning with mace and obi)
# ------------------------------------------------
# ─────────────────────────────────────────────────────────────────────────────
# Resolve absolute paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
# The Python package dir (torsionator/torsionator/) is what the .def binds
# as ../torsionator → /torsionator inside the container.
PKG_DIR="${REPO_ROOT}/torsionator"

# UMA model lives in the torsionator package folder on the host,
# which is bound to /torsionator inside the container.
UMA_MODEL_PATH_CONTAINER="/torsionator/uma-s-1p1.pt"

case "$METHOD" in
  obi|mace|uma) ;;
  *)
    echo "ERROR: unknown METHOD '${METHOD}'. Choose obi, mace, or uma." >&2
    exit 1
    ;;
esac

if [ "$METHOD" = "uma" ] && [ ! -f "${PKG_DIR}/uma-s-1p1.pt" ]; then
    echo "ERROR: uma-s-1p1.pt not found in ${PKG_DIR}/" >&2
    echo "       Place the UMA model file there before running with METHOD=uma." >&2
    exit 1
fi

# Resolve PDB_FILE_DIR to an absolute path (supports relative paths too)
PDB_FILE_DIR="$(cd "$PDB_FILE_DIR" && pwd)"

# Each NNP lives in its own conda env inside the container.
ENV_BIN="/opt/conda/envs/torsionator_${METHOD}/bin"
CONTAINER_PATH="${ENV_BIN}:/opt/conda/bin:/usr/local/cuda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"

# Build optional UMA env arg cleanly
UMA_ENV_ARG=""
if [ "$METHOD" = "uma" ]; then
    UMA_ENV_ARG="--env UMA_MODEL_PATH=${UMA_MODEL_PATH_CONTAINER}"
fi

apptainer exec \
  --nv \
  --bind "${PKG_DIR}:/torsionator" \
  --bind "${PDB_FILE_DIR}:/data" \
  --env "PATH=${CONTAINER_PATH}" \
  --env "PYTHONPATH=/" \
  --env "PYTHONUNBUFFERED=1" \
  --env "LD_LIBRARY_PATH=/usr/local/cuda/lib64:/usr/local/cuda-11.8/lib64" \
  --env "MPLBACKEND=Agg" \
  $UMA_ENV_ARG \
  "${SCRIPT_DIR}/torsionator.sif" \
  python3 -m "torsionator.${METHOD}_torsionator" \
    --pdb "/data/${PDB_FILE_NAME_ROOT}.pdb" \
    --dihedral "$DIHEDRAL" \
    --conf_analysis "$CONF_ANALYSIS" \
    --BCS "$BCS" \
    --MCS "$MCS" \
    --n_confs "$N_CONF" \
    --RMSD "$RMSD" \
    --multiplicity "$MULTIPLICITY" \
    --step_size "$STEP_SIZE" \
    --double_rotation "$DOUBLE_ROTATION" \
    --net_charge "$NET_CHARGE" \
    --spin "$SPIN"
