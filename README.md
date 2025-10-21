# torsionator                                                                                                                                                                               <img width="100" height="100" alt="logo" src="https://github.com/user-attachments/assets/e8f09f16-519e-436d-bff2-4764a1b889be" />

## 1. **Overview** <br>
Torsionator is an end‑to‑end pipeline for dihedral scans and torsion parameter fitting. It minimizes an input PDB using ML force fields (OBI/MACE), screens for steric clashes, optionally explores RDKit conformers, performs constrained scans, and fits torsional terms with AMBERTools' progam mdgx, finally writing an updated frcmod.


<img width="9369" height="3539" alt="Picture" src="https://github.com/user-attachments/assets/873f31b3-c82b-430b-9db0-4fc55123c327" />

## 2. **Instalaltion**
   
**Requirements** <br>
- Apptainer ≥ 1.x installed on the host<br>
- NVIDIA GPU (optional) and host NVIDIA drivers; use --nv if you want GPU acceleration <br>
- A <your_name> folder on the host that will be bind‑mounted as '/data' inside the container <br>

**Clone the repository**<br>
First clone the repo and then move into the top-level directory of the package.<br>
```
git clone https://github.com/giobros/torsionator.git
```

**Build the image**<br>
All the dependencies can be loaded together using the torsionator.sif generated with the .def file and Apptainer.
Enter the folder container and lunch the file .sh to create the image
```
cd torsionator/container
apptainer build --fakeroot torsionator.sif torsionator.def
```

## 3 **Prepare your host work directory**<br>
Place your pdb input and script inside a host directory that you’ll bind to /data, e.g.:
```
/$HOME/<your_folder>/
└── <ROOT>.pdb # your input structure
```

## 4 **Run (detached, with GPU)**<br>
The user can change the script run.sh inside the container folder to select which options apply to the scanning
In particular the modification to do are:
 - change the folder name *your_folder* with the actual folder name in --bind "$HOME/<your_folder>:/data" and the pdb *ROOT* in flag --pdb /data/<ROOT>.pdb 
 - change the scanning options:
```
   --method obi|mace|all\
   --dihedral all|[a,b,c,d]|"print" \
   --conf_analysis true|false \
   --force_scanning true|false
```
Suggestion: print the dihedrals before passing the wanted one, the code should recognize your input but my suggestion is to check the code-preferred dihedral definition.

## 5. **Where outputs are written**<br>
By default the script uses BASE_DIR = "/data"
You will find results under the following directories on the host inside your bound folder:
```
/<your_folder>/
├── conformers/
│   ├── pdb/*.pdb 
│   ├── OBI/
│   │   ├── initial_energies.txt
│   │   ├── optimized_energies.txt
│   │   ├── sorted_energies.txt
│   │   ├── min_energy.txt
│   │   ├── xyz/
│   │   └── *.pdb  
│   └── MACE/ ... (same layout)
├── scanning/
│     └── a_b_c_d/
│       ├── OBI/
│       │   ├── geometries.xyz
│       │   ├── angles_vs_energies.txt
│       │   ├── angles_vs_energies_final.txt   # sorted & min-shifted
│       │   ├── energies.dat                   # single-column, Hartree, min=0 
│       │   ├── scan_pdbs/*.pdb
│       │   └── obi.dat                        # MDGX torsion fit
│       ├── MACE/ ... (same layout)
│       └── a_b_c_d.png                        # plotted profile (kcal/mol)
└── parameters/<ROOT>_<method>_a_b_c_d.frcmod    # frcmod with updated DIHE lines


```

The workflow togheter with the rrors that stop it (e.g., clashes without --conf_scanning true) are written to:
```
/data/workflow.log
```
