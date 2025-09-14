import os, glob, shutil, numpy as np, logging
from dataclasses import dataclass
from rdkit import Chem
from rdkit.Chem import rdDistGeom, AllChem, rdMolAlign
from ase.io import read
from .io_utils import write_pdb, write_xyz_file
from .geometry import GeometryOptimizer
from .constants import CONV_EV_TO_EH

@dataclass
class ConformerGenerator:
    out_root: str

    def generate(self, input_pdb: str, numconf: int = 50) -> str:
        pdb_dir = os.path.join(self.out_root, "pdb"); os.makedirs(pdb_dir, exist_ok=True)
        mol = Chem.MolFromPDBFile(input_pdb, removeHs=False)
        cids = rdDistGeom.EmbedMultipleConfs(mol, numconf, rdDistGeom.ETKDGv2())
        props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
        AllChem.MMFFOptimizeMoleculeConfs(mol, mmffVariant="MMFF94s")

        conformers = [(cid, AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid).CalcEnergy()) for cid in cids]
        conformers.sort(key=lambda x: x[1])

        for idx, (cid, energy) in enumerate(conformers, 1):
            mol.SetProp("CID", str(cid)); mol.SetProp("Energy", str(energy))
            path = os.path.join(pdb_dir, f"{idx}.pdb")
            with Chem.PDBWriter(path) as w:
                w.write(mol, confId=cid)
        return pdb_dir

class ConformerScreen:
    @staticmethod
    def calc_rmsd(ref_pdb, target_pdb):
        ref = Chem.MolFromPDBFile(ref_pdb, removeHs=False)
        tgt = Chem.MolFromPDBFile(target_pdb, removeHs=False)
        ref_atoms = [a for a in ref.GetAtoms() if a.GetAtomicNum() > 1]
        tgt_atoms = [a for a in tgt.GetAtoms() if a.GetAtomicNum() > 1]
        if len(ref_atoms) != len(tgt_atoms):
            raise ValueError("Mismatch in heavy atom counts.")
        amap = [(ra.GetIdx(), ta.GetIdx()) for ra, ta in zip(ref_atoms, tgt_atoms)]
        return rdMolAlign.AlignMol(tgt, ref, atomMap=amap)

    def screen_by_rmsd(self, pdb_dir: str, threshold: float = 0.5):
        files = sorted([os.path.join(pdb_dir, f) for f in os.listdir(pdb_dir) if f.endswith(".pdb")],
                       key=lambda x: int(os.path.basename(x).split(".")[0]))
        kept = [files[0]]
        for cand in files[1:]:
            if all(self.calc_rmsd(ref, cand) >= threshold for ref in kept):
                kept.append(cand)
            else:
                os.remove(cand)

    def ensure_ref_if_unique(self, conformers_dir: str, reference_pdb: str, threshold: float = 0.5):
        pdbs = [os.path.join(conformers_dir, f) for f in os.listdir(conformers_dir) if f.endswith(".pdb")]
        rmsds = {p: self.calc_rmsd(reference_pdb, p) for p in pdbs}
        if rmsds and min(rmsds.values()) > threshold:
            shutil.copy(reference_pdb, os.path.join(conformers_dir, "0.pdb"))

@dataclass
class ConformerMinimizer:
    base_dir: str
    optimizer: GeometryOptimizer

    def minimize_folder(self, method: str, calc, pdb_dir: str) -> str:
        out_dir= os.path.join(self.base_dir, "conformers", method)
        out_dir_xyz = os.path.join(out_dir, "xyz")

        os.makedirs(out_dir, exist_ok=True)
        os.makedirs(out_dir_xyz, exist_ok=True)

        init_log = os.path.join(out_dir, "initial_energies.txt")
        opt_log  = os.path.join(out_dir, "optimized_energies.txt")

        files = sorted(glob.glob(os.path.join(pdb_dir, "*.pdb")),

                       key=lambda x: int(os.path.basename(x).split('.')[0]))

        min_energy = float("inf"); min_geom = None
        with open(init_log,"w") as fi, open(opt_log,"w") as fo:
            fi.write("File Initial_Energy in Eh\n")
            fo.write("File  Minimized_Energy in Eh\n")

            for pdb in files:
                atoms = read(pdb); atoms.calc = calc
                e0 = atoms.get_potential_energy() * CONV_EV_TO_EH
                fi.write(f"{os.path.basename(pdb)}\t{e0:.6f}\n")
                self.optimizer.minimize(atoms)
                emin = atoms.get_potential_energy() * CONV_EV_TO_EH
                fo.write(f"{os.path.basename(pdb)}\t{emin:.6f}\n")
                if emin < min_energy:
                    min_energy, min_geom = emin, atoms.copy()
                out_path = os.path.join(out_dir, os.path.basename(pdb))
                write_pdb(out_path, atoms)
                write_xyz_file(os.path.join(out_dir_xyz, os.path.basename(pdb).replace(".pdb",".xyz")), atoms, 'w')
        # sort energies for downstream iteration
        self._sort_energies(opt_log)
        return os.path.join(out_dir, "sorted_energies.txt")
    
    @staticmethod
    def _sort_energies(energy_file_path, output_sorted_file="sorted_energies.txt"):
        with open(energy_file_path,"r") as f:
            lines = f.readlines()
        entries = []
        for line in lines[1:]:
            parts = line.strip().split()
            if len(parts) >= 2:
                entries.append((parts[0], float(parts[1])))
        entries.sort(key=lambda x:x[1])
        out_path = os.path.join(os.path.dirname(energy_file_path), output_sorted_file)
        with open(out_path,"w") as fo:
            fo.write("FileName\tEnergy_Hartee\n")
            for fname, e in entries:
                fo.write(f"{fname}\t{e:.6f}\n")

def iter_energy_ordered_conformers(folder):
    sorted_file = os.path.join(folder, "sorted_energies.txt")
    if not os.path.exists(sorted_file):
        return
    with open(sorted_file) as f:
        for line in f.readlines()[1:]:
            fname = line.split()[0]
            path = os.path.join(folder, fname)
            if path.lower().endswith(".pdb") and os.path.exists(path):
                yield path
