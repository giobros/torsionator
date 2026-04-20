import glob
import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, List

import numpy as np
from ase.io import read
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom, rdMolAlign

from .calculators import attach_calc
from .constants import CONV_EV_TO_EH
from .geometry import GeometryOptimizer
from .io_utils import write_pdb, write_xyz_file, _int_stem_key


@dataclass
class ConformerGenerator:
    out_root: str

    def generate(self, input_pdb: str, n_confs: int) -> str:
        """
        Generate `n_confs` conformers with ETKDGv2 + MMFF94s, write sorted PDBs,
        and return the directory path containing them.
        """
        pdb_dir = os.path.join(self.out_root, "pdb")
        os.makedirs(pdb_dir, exist_ok=True)

        mol = Chem.MolFromPDBFile(input_pdb, removeHs=False, proximityBonding=False)
        cids = rdDistGeom.EmbedMultipleConfs(mol, n_confs, rdDistGeom.ETKDGv2())
        AllChem.MMFFOptimizeMoleculeConfs(mol, mmffVariant="MMFF94s")

        props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
        conformers = sorted(
            [(cid, AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid).CalcEnergy()) for cid in cids],
            key=lambda x: x[1],
        )

        for idx, (cid, energy) in enumerate(conformers, start=1):
            mol.SetProp("CID", str(cid))
            mol.SetProp("Energy", str(energy))
            with Chem.PDBWriter(os.path.join(pdb_dir, f"{idx}.pdb")) as w:
                w.write(mol, confId=cid)

        return pdb_dir

    @staticmethod
    def normalize_conformer(conf_dir: str) -> None:
        """
        Rename conformer PDBs in `conf_dir` so that:
          - 0.pdb (if present) is kept as-is
          - all others become 1.pdb, 2.pdb, ... (no gaps)
        """
        conf_dir_path = Path(conf_dir)

        def _sort_key(p: Path):
            return (0, int(p.stem)) if p.stem.isdigit() else (1, p.stem)

        all_pdbs = sorted(conf_dir_path.glob("*.pdb"), key=_sort_key)
        zero_path = next((p for p in all_pdbs if p.stem == "0"), None)
        others = [p for p in all_pdbs if p is not zero_path]

        # Rename to temp files to avoid collisions
        tmp_paths = []
        for j, p in enumerate(others, start=1):
            tmp = p.with_name(f"_tmp_{j}.pdb")
            p.rename(tmp)
            tmp_paths.append(tmp)

        for idx, tmp in enumerate(tmp_paths, start=1):
            tmp.rename(tmp.with_name(f"{idx}.pdb"))


class ConformerScreen:
    @staticmethod
    def _calc_rmsd(ref_pdb: str, target_pdb: str) -> float:
        ref = Chem.MolFromPDBFile(ref_pdb, removeHs=False, proximityBonding=False)
        tgt = Chem.MolFromPDBFile(target_pdb, removeHs=False, proximityBonding=False)
        ref_heavy = [a for a in ref.GetAtoms() if a.GetAtomicNum() > 1]
        tgt_heavy = [a for a in tgt.GetAtoms() if a.GetAtomicNum() > 1]
        if len(ref_heavy) != len(tgt_heavy):
            raise ValueError("Mismatch in heavy atom counts.")
        amap = [(r.GetIdx(), t.GetIdx()) for r, t in zip(ref_heavy, tgt_heavy)]
        return rdMolAlign.AlignMol(tgt, ref, atomMap=amap)

    def screen_by_rmsd(self, pdb_dir: str, threshold: float) -> None:
        """Remove conformers that are within `threshold` RMSD of any already-kept one."""
        files = sorted(
            [os.path.join(pdb_dir, f) for f in os.listdir(pdb_dir) if f.endswith(".pdb")],
            key=_int_stem_key,
        )
        kept = [files[0]]
        for cand in files[1:]:
            if all(self._calc_rmsd(ref, cand) >= threshold for ref in kept):
                kept.append(cand)
            else:
                os.remove(cand)

    def ensure_ref_if_unique(self, conformers_dir: str, reference_pdb: str, threshold: float) -> None:
        """If no conformer is within `threshold` of the reference, copy the reference as 0.pdb."""
        pdbs = [os.path.join(conformers_dir, f) for f in os.listdir(conformers_dir) if f.endswith(".pdb")]
        rmsds = {p: self._calc_rmsd(reference_pdb, p) for p in pdbs}
        if rmsds and min(rmsds.values()) > threshold:
            shutil.copy(reference_pdb, os.path.join(conformers_dir, "0.pdb"))


@dataclass
class ConformerMinimizer:
    base_dir: str
    optimizer: GeometryOptimizer

    def minimize_folder(self, method: str, calc, pdb_dir: str) -> str:
        """
        Minimize all PDBs in `pdb_dir` with `calc`, write results to
        `base_dir/conformers/<method>/`, and return the path to sorted_energies.txt.
        """
        out_dir = os.path.join(self.base_dir, "conformers", method)
        xyz_dir = os.path.join(out_dir, "xyz")
        os.makedirs(out_dir, exist_ok=True)
        os.makedirs(xyz_dir, exist_ok=True)

        init_log = os.path.join(out_dir, "initial_energies.txt")
        opt_log = os.path.join(out_dir, "optimized_energies.txt")

        pdb_files = sorted(glob.glob(os.path.join(pdb_dir, "*.pdb")), key=_int_stem_key)

        with open(init_log, "w") as fi, open(opt_log, "w") as fo:
            fi.write("File Initial_Energy in Eh\n")
            fo.write("File  Minimized_Energy in Eh\n")

            for pdb in pdb_files:
                atoms = read(pdb)
                attach_calc(atoms, calc)

                e0 = atoms.get_potential_energy() * CONV_EV_TO_EH
                fi.write(f"{os.path.basename(pdb)}\t{e0:.6f}\n")

                self.optimizer.minimize(atoms)
                e_min = atoms.get_potential_energy() * CONV_EV_TO_EH
                fo.write(f"{os.path.basename(pdb)}\t{e_min:.6f}\n")

                out_pdb = os.path.join(out_dir, os.path.basename(pdb))
                out_xyz = os.path.join(xyz_dir, os.path.basename(pdb).replace(".pdb", ".xyz"))
                write_pdb(out_pdb, atoms)
                write_xyz_file(out_xyz, atoms, mode="w")

        self._write_sorted_energies(opt_log)
        return os.path.join(out_dir, "sorted_energies.txt")

    @staticmethod
    def _write_sorted_energies(energy_file: str, output_name: str = "sorted_energies.txt") -> None:
        with open(energy_file) as f:
            lines = f.readlines()

        entries = []
        for line in lines[1:]:
            parts = line.strip().split()
            if len(parts) >= 2:
                entries.append((parts[0], float(parts[1])))
        entries.sort(key=lambda x: x[1])

        out_path = os.path.join(os.path.dirname(energy_file), output_name)
        with open(out_path, "w") as f:
            f.write("FileName\tEnergy_Hartee\n")
            for fname, e in entries:
                f.write(f"{fname}\t{e:.6f}\n")


def iter_energy_ordered_conformers(folder: str) -> Iterator[str]:
    """Yield PDB paths from `folder` in ascending energy order (via sorted_energies.txt)."""
    sorted_file = os.path.join(folder, "sorted_energies.txt")
    if not os.path.exists(sorted_file):
        return
    with open(sorted_file) as f:
        for line in f.readlines()[1:]:
            fname = line.split()[0]
            path = os.path.join(folder, fname)
            if path.lower().endswith(".pdb") and os.path.exists(path):
                yield path
