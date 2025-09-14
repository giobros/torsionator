import logging, glob, os
from dataclasses import dataclass
from .dihedral import ClashAnalyzer
from .lj import local_rg
from ase.io import read

@dataclass
class ConformerSelector:
    analyzer: ClashAnalyzer

    def select_best(self, conf_folder: str, dihedral, group_A, group_D, steps=37):
        results = []
        for p in sorted(glob.glob(os.path.join(conf_folder,"*.pdb")),
                        key=lambda x: int(os.path.basename(x).split('.')[0])):
            atoms = read(p)
            res = self.analyzer.scan_check(atoms, dihedral, group_A, group_D, steps=steps)
            results.append({"path": p, "lj": res.max_lj_kcal, "rg": res.max_av_rg})
            #logging.info("→ %s | Max LJ = %.3f kcal/mol | ⟨Rg⟩ = %.3f Å", os.path.basename(p), res.max_lj_kcal, res.mean_rg)
        best_lj = min(results, key=lambda r: r["lj"])
        best_rg = min(results, key=lambda r: r["rg"])
        return {"best_lj": best_lj, "best_rg": best_rg, "all_results": results}
