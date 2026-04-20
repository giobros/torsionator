import glob
import os
from dataclasses import dataclass
from typing import Dict, List, Tuple

from ase.io import read

from .dihedral import ClashAnalyzer
from .io_utils import _int_stem_key


@dataclass
class ConformerSelector:
    analyzer: ClashAnalyzer

    def select_best(
        self,
        conf_folder: str,
        dihedral: Tuple[int, int, int, int],
        group_A: List[int],
        group_D: List[int],
        steps: int = 37,
    ) -> Dict:
        """
        Scan all conformers in `conf_folder` and return the one with the
        lowest max LJ energy (best_lj) and lowest radius of gyration (best_rg).
        """
        pdb_files = sorted(glob.glob(os.path.join(conf_folder, "*.pdb")), key=_int_stem_key)

        results = []
        for path in pdb_files:
            atoms = read(path)
            res = self.analyzer.scan_check(atoms, dihedral, group_A, group_D, steps=steps)
            results.append({"path": path, "lj": res.max_lj_kcal, "rg": res.max_av_rg})

        best_lj = min(results, key=lambda r: r["lj"])
        best_rg = min(results, key=lambda r: r["rg"])
        return {"best_lj": best_lj, "best_rg": best_rg, "all_results": results}
