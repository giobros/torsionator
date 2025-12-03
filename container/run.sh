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
   --pdb /data/<ROOT>.pdb \   ## MODIFY <ROOT> WITH YOU ACTUAL PDB NAME
   --method obi|mace|all \      
   --dihedral all|[a,b,c,d]|"print" \   # INDICES 0-BASED
   --conf_analysis true|false \
   --BCS true|false # Best-Conformer-Scan
   --MCS true|false # Multi-Conformer-Scan
   --N # N. of conformers to generate
   --RMSD # RMSD pruning threshold


