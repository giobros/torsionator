# torsionator
## 1. **Overview** <br>
Torsionator is an end‑to‑end pipeline for dihedral scans and torsion parameter fitting. It minimizes an input PDB using ML force fields (OBI/MACE), screens for steric clashes, optionally explores RDKit conformers, performs constrained scans, and fits torsional terms with AMBERTools' progam mdgx, finally writing an updated frcmod.


<img width="9428" height="3573" alt="Picture" src="https://github.com/user-attachments/assets/90083681-c90d-4aa8-ac53-22760053b18e" />
## 2. **Instalaltion**
   
**Requirements** <br>
- Apptainer ≥ 1.x installed on the host<br>
  (Docker if installation starts from the Dockerfile)<br>
- NVIDIA GPU (optional) and host NVIDIA drivers; use --nv if you want GPU acceleration <br>
- A 'your_name' folder on the host that will be bind‑mounted as '/data' inside the container <br>

***Clone the repository**
First clone the repo and then move into the top-level directory of the package.<br>
```
git clone https://github.com/giobros/torsionator.git
```

**Build the image**
All the dependencies can be loaded together using the torsionator.sif generated with the Dockerfile and Apptainer.
Enter the folder container and lunch the file .sh to create the image
```
cd container
docker build -t ubuntu22_cuda11.2 .
docker save -o cuda11.2.tar ubuntu22_cuda11.2
apptainer build torsionator.sif docker-archive://./cuda11.2.tar
``` 
## 3 **Prepare your host work directory**
Place your pdb input and script inside a host directory that you’ll bind to /data, e.g.:
```
$HOME/your_folder/
├── file.pdb # your input structure
```
4 ##**Run (detached, with GPU)**
The user can change the script run.sh to select which options apply to the scanning

```
#!/bin/bash

unset PYTHONPATH  # Prevent host leakage

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
```

## 5. **Where outputs are written**
By default the script uses BASE_DIR = "/data"
You will find results under the following directories on the host inside your bound folder:
```
/data/
├── conformers/
│   ├── OBI/
│   │   ├── initial_energies.txt
│   │   ├── optimized_energies.txt
│   │   ├── sorted_energies.txt
│   │   ├── min_energy.txt
│   │   └── *.pdb / *.xyz (minimized conformers)
│   └── MACE/ ... (same layout)
└── scanning/
    └── a_b_c_d/
        ├── OBI/
        │   ├── geometries.xyz
        │   ├── angles_vs_energies.txt
        │   ├── angles_vs_energies_final.txt   # sorted & min-shifted
        │   ├── energies.dat                    # single-column, Hartree, min=0 (MDGX)
        │   ├── scan_pdbs/*.pdb
        │   ├── obi.dat                         # MDGX torsion fit
        │   └── final_OBI.frcmod                # frcmod with updated DIHE lines
        ├── MACE/ ... (same layout)
        └── a_b_c_d.png                         # plotted profile (kcal/mol, min=0)

```

The workflow togheter with the rrors that stop it (e.g., clashes without --conf_scanning true) are written to:
```
/data/workflow.log
```
