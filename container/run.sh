#!/bin/bash

unset PYTHONPATH

apptainer exec \
  --nv \
  --bind "$HOME/torsionator/torsionator:/torsionator" \
  --bind "$HOME/<your_folder>:/data" \     ## MODIFY <your_folder> WITH YOU ACTUAL FOLDER NAME
  --env PYTHONPATH=/ \
  --env HOME=/root \
  torsionator.sif \
  python3.9 -m torsionator.cli \
   --pdb /data/<ROOT>.pdb \   # MODIFY <ROOT> WITH YOU ACTUAL PDB NAME
   --method obi|mace|all   #default=all   
   --dihedral all|[a,b,c,d]|"print"  #default=all \   # INDICES 0-BASED \
   --conf_analysis true|false|none #default=false \
   --BCS true|false|none  # Best-Conformer-Scan default=false\
   --MCS true|false|none # Multi-Conformer-Scan  default=false \
   --N # N. of conformers to generate #default=50 \
   --RMSD # RMSD pruning threshold  #default=0.5 \


