<img width="2227" height="1856" alt="papero" src="https://github.com/user-attachments/assets/a5b8ab19-13fc-400d-a98a-80d6539d4581" /># torsionator_plus
<p align="center">
<img width="400" height="400" alt="logo" src="https://github.com/user-attachments/assets/40132121-0ba2-44f2-8177-a5e1b68d2f5b" />
</p>

## Torsionator Plus  <br>

## 1. **Overview** <br>
Torsionator is an end‑to‑end pipeline for dihedral scans and torsion parameter fitting. It minimizes an input PDB using ML force fields (OBI/MACE), screens for steric clashes, optionally explores RDKit conformers, performs constrained scans, and fits torsional terms with AMBERTools' progam mdgx, finally writing an updated frcmod.<br>

<img width="8285" height="6592"  alt="papero" src="https://github.com/user-attachments/assets/aa15a552-6d3e-485b-bf23-192146a3758a" />

## 2. **Instalaltion**

**Clone the repository**<br>
First clone the repo and then move into the top-level directory of the package.<br>
```
git clone https://github.com/giobros/torsionator_plus.git
```

**Create the enviroments**<br>

 Create conda env — one per backend (MACE / OBI / UMA), pointing to the right environment.yaml
```
conda env create -f <method>/enviroment.yaml

```
<method> available: uma / mace (mace_off23) /obi

## 3 **Prepare your host work directory**<br>
Place your pdb input and script inside the "main" folder insaide the repo:
```
/<your_folder_path>/
└── <NAME>.pdb # your input structure
```

## 4 **Run**<br>
The user can change and use the script run.sh inside the "main" folder to select which options apply to the scanning.
To change the scanning options modify:
```

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
NET_CHARGE="0"                   # net molecule charge, (default 0, other charges are not accepted by mace and obi) 
# ------------------------------------------------
```

Suggestion: print the dihedrals before passing the wanted one, the code should recognize your input but my suggestion is to check the code-preferred dihedral definition.

## 5. **Where outputs are written**<br>

You will find results under the following directories on the host inside your PDB_FILE_DIR:

If BCS=true
```
/<your_folder_path>/
├── conformers/
│   ├── pdb/*.pdb 
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
│        └── n_confs folders/minimized.pdb  
├── scanning/
│     └── a_b_c_d/
│       ├── method/                     
│       │   ├── n_confs folders (+ n_confs_ccw)/scan_pdbs/*.pdb
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
