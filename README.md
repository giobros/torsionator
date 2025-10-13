# torsionator
1. ##**Overview**
Torsionator is an end‑to‑end pipeline for dihedral scans and torsion parameter fitting. It minimizes an input PDB using ML force fields (OBI/MACE), screens for steric clashes, optionally explores RDKit conformers, performs constrained scans, and fits torsional terms with AMBERTools' progam mdgx, finally writing an updated frcmod.
<img width="9428" height="3573" alt="Picture" src="https://github.com/user-attachments/assets/90083681-c90d-4aa8-ac53-22760053b18e" />
2. ##**Instalaltion**
   
**Requirements**
Apptainer ≥ 1.x installed on the host (also Docker if installation starts from the Dockerfile)
NVIDIA GPU (optional) and host NVIDIA drivers; use --nv if you want GPU acceleration
A folder on the host that will be bind‑mounted as /data inside the container

***Clone the repository**
First clone the repo and then move into the top-level directory of the package.
$git clone https://github.com/giobros/torsionator.git

**Build the image**
All the dependencies can be loaded together using the image.sif generated with the Dockerfile and Apptainer.
Enter the folder container and lunch the file .sh to create the image
$cd container
$docker build -t ubuntu22_cuda11.2 .
$docker save -o cuda11.2.tar ubuntu22_cuda11.2
$apptainer build image.sif docker-archive://./cuda11.2.tar

3 ## **Prepare your host work directory**
Place your inputs and script inside a host directory that you’ll bind to /data, e.g.:
```
```bash
$HOME/TOOL_NEW_EMI/prova/
├── MKC.pdb # your input structure
├── tool.py # the main script (entry point)
└── ... # optional: antechamber_leap.py, model files, etc.
```
