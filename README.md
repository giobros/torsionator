<p align="center">
<img width="520" height="200" alt="logo" src="https://github.com/user-attachments/assets/82db9a0e-aeaa-40e5-9831-a50ceeeb0c2c" />
</p>

## Torsionator  <br>

## 1. **Overview** <br>
Torsionator is an end‑to‑end pipeline for dihedral scans and torsion parameter fitting. It minimizes an input PDB using ML force fields (OBI/MACE/UMA), screens for steric clashes, optionally explores RDKit conformers, performs constrained scans, and fits torsional terms with AMBERTools' progam mdgx, finally writing an updated frcmod.<br>

<img width="9285" height="8002" alt="Picture3" src="https://github.com/user-attachments/assets/25984aaf-26c1-41f0-a0ab-389f4adeca56" />


## 2. **Installation**

**Requirements** <br>
- Apptainer ≥ 1.x installed <br>
- NVIDIA GPU (optional) and host NVIDIA drivers; use --nv if you want GPU acceleration <br>
- <your_folder_path> on the host that will be bind‑mounted as '/data' inside the container <br>

To install Apptainer with the correct privileged version look at https://apptainer.org/docs/admin/1.4/installation.html, if you have ubuntu you can look install like that:
```
wget https://github.com/apptainer/apptainer/releases/download/v1.4.5/apptainer_1.4.5_amd64.deb
sudo apt install -y ./apptainer_1.4.5_amd64.deb

wget https://github.com/apptainer/apptainer/releases/download/v1.4.5/apptainer-suid_1.4.5_amd64.deb
sudo dpkg -i ./apptainer-suid_1.4.5_amd64.deb

```
**Clone the repository**<br>
First clone the repo and then move into the top-level directory of the package.<br>
```
git clone https://github.com/giobros/torsionator.git
```
**Build the image**<br>
All the dependencies can be loaded together using the torsionator.sif generated with the .def file and Apptainer.
You can create the .sif in two ways:
  1. Enter the folder container and lunch the file .sh to create the image
  ```
  cd torsionator/container
  sudo apptainer build torsionator.sif torsionator.def
  ```
  2. To obtain the image faster and in a more reproducible way, you can simply concatenate the split SIF parts (available under the Releases section of this GitHub repository) instead of rebuilding it from the .def file: 
  ```
  cat torsionator.sif.part_a* > torsionator.sif
  ```

The NNPs available are: uma (small) / mace (mace_off23) /obi

For UMA: you need to install from https://huggingface.co/facebook/UMA the checkpoints "uma-s-1p1.pt", and put in folder "torsionator" of the repository, there the pipeline scripts are located.

## 3 **Prepare your host work directory**<br>
Place your pdb input a folder /<your_folder_path>/:
```
/<your_folder_path>/
└── <NAME>.pdb # your input structure
```
Torsionator uses Open Babel and RDKit, and they  rely on strict atom typing and element recognition, input PDB files were curated to remove ambiguities in atom naming and to explicitly define element fields; for instance, a chlorine atom originally labeled as ‘Cl’ in the atom name field was standardized by explicitly specifying ‘Cl’ in the element column to ensure correct parsing.

## 4 **Run**<br>
The user can change and use the script container/run.sh of the repository to select which options apply to the scanning.
To change the scanning options modify:

```

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
```

Look at container/run_example.sh for a "real" case example, where the input .pbd (fragment.pdb) used is located in folder example/, same folder where Torsionator generated its outputs . 

Suggestion: print the dihedrals before passing the wanted one, the code should recognize your input but my suggestion is to check the code-preferred dihedral definition.

## 5. **Where outputs are written**<br>

You will find results under the following directories on the host inside your PDB_FILE_DIR:

If BCS=true
```
/<your_folder_path>/
├── conformers/
│   ├── pdb/n.pdb  n= n_conf + 1 (input pdb)   
│   └── method/
│      ├── initial_energies.txt
│      ├── optimized_energies.txt
│      ├── sorted_energies.txt
│      ├── min_energy.txt
│      ├── xyz/
│      └── *.pdb  
│   
├── scanning/
│     └── a_b_c_d/
│       ├── method/
│       │   ├── geometries.xyz
│       │   ├── angles_vs_energies.txt
│       │   ├── angles_vs_energies_final.txt   # sorted & min-shifted
│       │   ├── energies.dat                   # single-column, Hartree, min=0 
│       │   ├── scan_pdbs/*.pdb
│       │   └── output                       # MDGX torsion fit
│       │         └──hrst/hrst.dat
│       ├── GAFF2/old/method/           ← old GAFF2 fit
│       ├── GAFF2/new/method/  
│       └── a_b_c_d.png                        # plotted profile (kcal/mol)
└── parameters/<NAME>_<method>.frcmod  # frcmod with updated DIHE lines


```
If MCS = true:
```
/<your_folder_path>/
├── conformers/
│   ├── pdb/*.pdb 
│   └──  method/      
│        └── n folders/minimized.pdb  n= n_conf + 1 (input pdb)
├── scanning/
│     └── a_b_c_d/
│       ├── method/                     
│       │   ├── n folders (+ n_ccw)/scan_pdbs/*.pdb    n= n_conf + 1 (input pdb)
│       │   └──  MCS/
│       │       ├── angles_vs_energies_final.txt   # sorted & min-shifted
│       │       ├── geometries.xyz
│       │       ├── energies.dat
│       │       └── output                       # MDGX torsion fit
│       │         └──hrst/hrst.dat
│       ├── GAFF2/old/method/           ← old GAFF2 fit
│       ├── GAFF2/new/method/  
│       └── a_b_c_d_MCS.png                        # plotted profile (kcal/mol)
└── parameters/<NAME>_<method>.frcmod  # frcmod with updated DIHE lines


```
The workflow togheter with the errors are written to:
```
/your_folder_path/workflow.log
```

Note: If mdgx has already generated .dat and .out files in the folder, delete them before rerunning. Otherwise, errors may occur and the pipeline could stop unexpectedly.
