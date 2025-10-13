#!/bin/bash

unset PYTHONPATH

apptainer exec \
  --nv \
  --bind "$HOME/torsionator/torsionator:/torsionator" \
  --bind "$HOME/your_folder:/data" \
  --env PYTHONPATH=/ \
  --env HOME=/root \
  torsionator.sif \
  python3.9 -m torsionator.cli \
   --pdb /data/file.pdb \
   --method obi|mace|both \
   --dihedral all|[a,b,c,d]|print \
   --conf_analysis true|false \
   --conf_scanning true|false
