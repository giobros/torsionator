import os
import time
import shutil
import subprocess
from dataclasses import dataclass, field
from typing import List, Tuple, Dict, Set, Optional
from ase.io import read
from rdkit import Chem

from .config import Config
from .logging_utils import setup_logging
from .calculators import CalculatorFactory
from .geometry import GeometryOptimizer
from .conformers import (
    ConformerGenerator,
    ConformerScreen,
    ConformerMinimizer,
    iter_energy_ordered_conformers,
)
from .dihedral import (
    DihedralFinder,
    GraphSplitter,
    ClashAnalyzer,
    DihedralScanner,
    plot_shifted_energy_profiles,
)
from .selection import ConformerSelector
from .paramfit.cpptraj import write_xyz_to_nc_in
from .paramfit.amber import get_atom_types_from_mol2, param_input_file, update_frcmod


@dataclass
class Workflow:
    cfg: Config
    # Track a cumulative final frcmod per method 
    final_frcmod: Dict[str, str] = field(default_factory=dict)

    def __post_init__(self):
        self.log = setup_logging(self.cfg.log_file)
        self.opt = GeometryOptimizer()
        self.calc_factory = CalculatorFactory(self.cfg.obiwan_path)
        self.clash = ClashAnalyzer()

    # ---------- pretty logging helpers ----------
    def _banner(self, title: str, char: str = "="):
        line = char * max(60, len(title) + 10)
        self.log.info("%s", line)
        self.log.info("  %s", title)
        self.log.info("%s", line)

    def _sub(self, title: str, char: str = "-"):
        line = char * (len(title) + 4)
        self.log.info("%s", line)
        self.log.info("  %s", title)
        self.log.info("%s", line)

    def _choice(
        self,
        m: str,
        d_tup: Tuple[int, int, int, int],
        used: str,
        reason: str,
    ):
        """Standardized one-liner for 'what did we use and why'."""
        self.log.info("[%s][%s] USED: %s | reason: %s", m, d_tup, used, reason)

    # -------- RDKit validation helper --------
    def _validate_and_expand_dihedral_with_rdkit(
        self, pdb_file: str, a: int, b: int, c: int, d: int
    ) -> List[Tuple[int, int, int, int]]:
        """
        Validate (a,b,c,d) describes a real dihedral in the molecule using RDKit.
        Assumes the user provides **0-based indices** (Python style). Raises ValueError if invalid.
        """
        mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
        if mol is None:
            raise ValueError(f"RDKit could not read PDB file: {pdb_file}")
        n = mol.GetNumAtoms()
        user_idxs = [a, b, c, d]
        def bonded(i, j):
            return mol.GetBondBetweenAtoms(i, j) is not None
        if not (bonded(a, b) and bonded(b, c) and bonded(c, d)):
            raise ValueError("Indices do not form a dihedral angle")
        forward = tuple(user_idxs)
        return [forward]

    # -------- args parsing for dihedrals --------
    def resolve_dihedrals_arg(self, dihedral_arg: str, pdb_file: str):
        if dihedral_arg == "all":
            return DihedralFinder.find_rotables(pdb_file), False

        if dihedral_arg == "print":
            dihedrals = DihedralFinder.find_rotables(pdb_file)
            self.log.info("Dihedrals in the molecule:")
            for d in dihedrals:
                self.log.info("[%d,%d,%d,%d]", d[0], d[1], d[2], d[3])
            return dihedrals, True

        if dihedral_arg.startswith("[") and dihedral_arg.endswith("]"):
            a, b, c, d = [int(x.strip()) for x in dihedral_arg.strip("[]").split(",")]
            try:
                all_defs = self._validate_and_expand_dihedral_with_rdkit(
                    pdb_file, a, b, c, d
                )
            except Exception as e:
                self.log.error("Validation for --dihedral %s: %s", dihedral_arg, e)
                raise
            return all_defs, False

        raise ValueError("Invalid --dihedral argument.")

    # -------- minimization of input structure per method --------
    def minimize_input(self, pdb_file: str, methods: List[str]) -> Dict[str, str]:
        out: Dict[str, str] = {}
        for m in methods:
            self.log.info("→ Minimizing input with %s", m)
            calc, label = self.calc_factory.get(m)
            out_path = os.path.join(self.cfg.base_dir, f"minimized_{label}.pdb")
            self.opt.minimize_and_write_pdb(pdb_file, calc, out_path)
            out[m] = out_path
        return out

    # -------- parametrization for a single dihedral (per method) --------
    def parameterize_one_dihedral(
        self,
        ROOT: str,
        m: str,
        label: str,
        d_tup: Tuple[int, int, int, int],
        pdb_source: str,
    ):
        """
        For a single dihedral and method:
        - prepare scan folder
        - run antechamber+LEaP, cpptraj, mdgx
        - merge torsion params
        - write two outputs under parameters/:
            1) per-dihedral: {ROOT}_{m}_{a}_{b}_{c}_{d}.frcmod
            2) cumulative : {ROOT}_{m}.frcmod  (updated atomically each call)
        """
        dstr = "_".join(map(str, d_tup))
        scan_dir = os.path.join(self.cfg.base_dir, "scanning", dstr)
        os.makedirs(scan_dir, exist_ok=True)

        try:
            shutil.copy(pdb_source, scan_dir)
        except Exception as e:
            self.log.error("Failed to copy source PDB for %s → %s: %s", d_tup, scan_dir, e)
            return

        ante_dir = os.path.join(scan_dir, f"antechamber_{label}")
        os.makedirs(ante_dir, exist_ok=True)
        shutil.copy(pdb_source, os.path.join(ante_dir, "input.pdb"))

        # 1) antechamber + LEaP (module lives in this repo)
        try:
            subprocess.run(
                ["python3", "-m", "torsionator.antechamber_leap", scan_dir, label],
                cwd=scan_dir,
                check=True,
            )
        except subprocess.CalledProcessError as e:
            self.log.error("Antechamber/LEaP failed for %s: %s", d_tup, e)
            return

        # 2) xyz -> nc via cpptraj (uses files generated by scan in <scan_dir>/<label>)
        write_xyz_to_nc_in(scan_dir, label)
        method_dir = os.path.join(scan_dir, label)
        try:
            subprocess.run(["cpptraj", "-i", "xyz_to_nc.in"], cwd=method_dir, check=True)
        except subprocess.CalledProcessError as e:
            self.log.error("cpptraj failed for %s: %s", d_tup, e)
            return

        # 3) mdgx: build input and fit
        mol2_path = os.path.join(ante_dir, "MOL.mol2")
        if not os.path.exists(mol2_path):
            self.log.error("Missing MOL.mol2 for %s at %s", d_tup, mol2_path)
            return

        atom_types = get_atom_types_from_mol2(mol2_path)
        param_input_file(atom_types, m, d_tup, scan_dir, label)

        pfile = f"parametrization_{label}_dihedral_{dstr}.in"
        try:
            subprocess.run(["mdgx", "-i", pfile], cwd=method_dir, check=True)
        except subprocess.CalledProcessError as e:
            self.log.error("mdgx failed for %s: %s", d_tup, e)
            return

        # 4) merge updated torsional params back into frcmod (for this dihedral only)
        frcmod_path = os.path.join(ante_dir, "MOL.frcmod")           # seed frcmod from parmchk2
        dat_path = os.path.join(method_dir, f"{m.lower()}.dat")      # mdgx fit results

        # Produce a merged frcmod in method_dir (update_frcmod may write new frcmod)
        update_frcmod(frcmod_path, dat_path, m, method_dir)
        self.log.info("Per-dihedral merge done for %s (%s)", d_tup, label)

        # Normalize any temp name produced in method_dir
        maybe_new = os.path.join(method_dir, f"new.frcmod")
        method_frcmod = (
            os.path.join(method_dir, f"{m}.frcmod")
            if os.path.exists(os.path.join(method_dir, f"{m}.frcmod"))
            else os.path.join(method_dir, "frcmod")
        )
        if os.path.exists(maybe_new):
            try:
                # prefer a stable name '{m}.frcmod' inside method_dir
                os.replace(maybe_new, os.path.join(method_dir, f"{m}.frcmod"))
                method_frcmod = os.path.join(method_dir, f"{m}.frcmod")
            except Exception as e:
                self.log.warning("Could not normalize %s: %s", maybe_new, e)
                # if normalization fails, fall back to the temp file path
                method_frcmod = maybe_new

        # 5) write the two FINAL deliverables under parameters/
        final_dir = os.path.join(self.cfg.base_dir, "parameters")
        os.makedirs(final_dir, exist_ok=True)

        # 5a) PER-DIHEDRAL: {ROOT}_{m}_{a}_{b}_{c}_{d}.frcmod
        per_dihedral_final = os.path.join(final_dir, f"{ROOT}_{m}_{dstr}.frcmod")

        # choose best available source for the per-dihedral output
        candidates = [
            os.path.join(method_dir, f"{m}.frcmod"),
            os.path.join(method_dir, "frcmod"),
            os.path.join(method_dir, "MOL.frcmod"),
            frcmod_path,  # last resort: pre-merged parmchk2 output
        ]
        src_found = next((p for p in candidates if os.path.exists(p)), None)
        if src_found is None:
            self.log.error("No frcmod produced for %s; candidates missing in %s", d_tup, method_dir)
            return

        try:
            shutil.copy(src_found, per_dihedral_final)
            self.log.info("Saved per-dihedral final frcmod → %s", per_dihedral_final)
        except Exception as e:
            self.log.error("Failed to save per-dihedral frcmod %s: %s", per_dihedral_final, e)
            return

        # 5b) CUMULATIVE: {ROOT}_{m}.frcmod  (seed once, then merge in)
        cumulative_frcmod = self.final_frcmod.get(m)
        if cumulative_frcmod is None:
            cumulative_frcmod = os.path.join(final_dir, f"{ROOT}_{m}.frcmod")
            try:
                shutil.copy(per_dihedral_final, cumulative_frcmod)
            except Exception as e:
                self.log.error("Failed to seed cumulative frcmod for %s: %s", m, e)
            else:
                self.final_frcmod[m] = cumulative_frcmod

        # Merge the new .dat into the cumulative and atomically replace if a temp is produced
        if self.final_frcmod.get(m):
            try:
                update_frcmod(cumulative_frcmod, dat_path, m, final_dir)

                temp_new = os.path.join(final_dir, f"new_{m}.frcmod")
                if os.path.exists(temp_new):
                    try:
                        # Atomically replace cumulative and remove temp
                        os.replace(temp_new, cumulative_frcmod)
                    except Exception as e:
                        self.log.warning("Could not finalize cumulative update with %s: %s", temp_new, e)
                        # best-effort cleanup
                        try:
                            if os.path.exists(temp_new):
                                os.remove(temp_new)
                        except Exception:
                            pass
            except Exception as e:
                self.log.error("Failed to update cumulative frcmod for %s: %s", m, e)

    # -------- main workflow --------
    def run(
        self,
        pdb_file: str,
        method: str,
        dihedral_arg: str,
        conf_analysis: Optional[str] = None,
        force_scan: Optional[str] = None,
    ):
        """
        dihedral_arg:
          - 'all'      -> scan all rotables (as provided by DihedralFinder)
          - 'print'    -> just list and exit
          - '[a,b,c,d]' (0-based indices) 
        """
        _ = read(pdb_file)
        ROOT = os.path.splitext(os.path.basename(pdb_file))[0]

        def _norm_flag(x):
            if x is None:
                return None
            s = str(x).strip().lower()
            return None if s in ("", "none", "null") else s

        conf_norm = _norm_flag(conf_analysis)
        force_norm = _norm_flag(force_scan)

        # ---- headers & config summary
        self._banner(f"WORKFLOW START for {ROOT}")
        self.log.info("Methods: %s", "all (obi, mace)" if method == "all" else method)
        self.log.info("Flags: conf_analysis=%r, force_scan=%r", conf_norm, force_norm)

        dihedrals, should_exit = self.resolve_dihedrals_arg(dihedral_arg, pdb_file)
        if should_exit:
            self._banner("WORKFLOW END (print-only)")
            return

        methods = ["obi", "mace"] if method == "all" else [method]

        self._banner("STEP 1: Resolve target dihedrals")
        self.log.info("Target dihedrals (0-based): %s", sorted(map(tuple, dihedrals)))

        # step 1: input minimization
        self._banner("STEP 2: Input minimization")
        start = time.time()
        minimized = self.minimize_input(pdb_file, methods)
        self.log.info("Minimization completed in %.2f s", time.time() - start)
        for m in methods:
            self.log.info("[%s] Minimized PDB → %s", m, minimized[m])

        # step 2: clash detection over target dihedrals
        self._banner("STEP 3: Clash detection")
        splitter = GraphSplitter()
        all_diheds: Set[Tuple[int, int, int, int]] = set(map(tuple, dihedrals))
        clash_dihedrals: Set[Tuple[int, int, int, int]] = set()

        for d in dihedrals:
            d_tup = tuple(d)
            # Check clashes across methods; if any method clashes, mark as clash
            for m in methods:
                atoms_min = read(minimized[m])
                gA, gD = splitter.split_by_dihedral(atoms_min, d_tup, method="neighbors")
                res = self.clash.scan_check(atoms_min, d_tup, gA, gD)
                self.log.info("[%s][%s] clash: %s", m, d_tup, "NO" if res.no_clash else "YES")
                if not res.no_clash:
                    clash_dihedrals.add(d_tup)
                    break  # no need to check other methods for this dihedral

        non_clash_dihedrals = sorted(all_diheds - clash_dihedrals)
        self.log.info("Clash dihedrals: %s", sorted(clash_dihedrals))
        self.log.info("Non-clash dihedrals: %s", non_clash_dihedrals)

        use_conformers = False
        best_conf: Dict[str, Dict[Tuple[int, int, int, int], str]] = {m: {} for m in methods}

        # step 3: optionally rescue clashing torsions with conformers
        if clash_dihedrals:
            self._banner("STEP 4: Clash handling & conformer strategy")

            # If the flag is ABSENT (None after normalization) → explicit error so “nothing shows up” never happens
            if conf_norm is None:
                msg = (
                    f"Clashes detected for {sorted(clash_dihedrals)}. "
                    "No --conf_analysis provided. "
                    "Pass --conf_analysis true to attempt conformer rescue, "
                    "or --conf_analysis false to force scanning despite clashes."
                )
                self.log.error(msg)
                raise RuntimeError(msg)

            if conf_norm not in ("true", "false"):
                msg = f"Invalid --conf_analysis value: {conf_analysis!r}. Use 'true' or 'false'."
                self.log.error(msg)
                raise RuntimeError(msg)

            if conf_norm == "false":
                # Force scanning even if clashing
                self._sub("Conformer rescue disabled (conf_analysis=false) → Scanning anyway")
                dihedrals_to_scan = sorted(all_diheds)
            else:  # conf_norm == "true"
                use_conformers = True
                dihedrals_to_scan = sorted(all_diheds)

                # generate & minimize conformers once
                self._sub("Conformer rescue enabled (conf_analysis=true)")
                gen = ConformerGenerator(os.path.join(self.cfg.base_dir, "conformers"))
                conf_dir = gen.generate(pdb_file)
                self.log.info("Conformers directory: %s", conf_dir)

                screen = ConformerScreen()
                screen.screen_by_rmsd(conf_dir)
                screen.ensure_ref_if_unique(conf_dir, pdb_file)

                minim = ConformerMinimizer(self.cfg.base_dir, self.opt)
                for m in methods:
                    calc, _ = self.calc_factory.get(m)
                    minim.minimize_folder(m, calc, conf_dir)
                    self.log.info("[%s] Conformers minimized", m)

                selector = ConformerSelector(self.clash)

                # choose best conformer ONLY for clashing dihedrals
                for m in methods:
                    self.log.info("Selecting conformers for %s ...", m)
                    # Need gA/gD for selection criteria → compute from minimized structure
                    atoms_min = read(minimized[m])
                    for d in sorted(clash_dihedrals):
                        d_tup = tuple(d)
                        gA, gD = splitter.split_by_dihedral(atoms_min, d_tup)
                        conf_method_dir = os.path.join(self.cfg.base_dir, "conformers", m)
                        # call select_best once per dihedral
                        res = selector.select_best(conf_method_dir, d_tup, gA, gD)
                        if res.get("clash_free"):
                            best_conf[m][d_tup] = res["clash_free"]["path"]
                            self._choice(
                                m, d_tup,
                                used="clash-free conformer",
                                reason="lowest energy among clash-free",
                            )
                        else:
                            # Fallbacks: print Best-LJ first, then Best-Rg (if present)
                            if force_norm == "true" and res.get("best_lj"):
                                # Record the used conformer (Best-LJ)
                                best_conf[m][d_tup] = res["best_lj"]["path"]
                                self._choice(
                                    m, d_tup,
                                    used="best-LJ conformer",
                                    reason="no clash-free available; force_scan=true allows LJ fallback",
                                )
                                # Explicit candidate order (avoid repetition elsewhere):
                                self.log.info("[%s][%s] Candidate (Best-LJ): %s", m, d_tup, res["best_lj"]["path"])
                                if res.get("best_rg"):
                                    self.log.info("[%s][%s] Candidate (Best-Rg): %s", m, d_tup, res["best_rg"]["path"])
                            else:
                                self.log.error(
                                    "[%s][%s] No clash-free conformer. best-LJ=%s, best-Rg=%s. "
                                    "Pass --force_scan true to allow LJ fallback.",
                                    m, d_tup, bool(res.get("best_lj")), bool(res.get("best_rg"))
                                )
        else:
            dihedrals_to_scan = sorted(all_diheds)

        # Fail-fast guard: never allow a silent no-op
        if not dihedrals_to_scan:
            msg = (
                "No dihedrals selected for scanning. "
                f"clash_dihedrals={sorted(clash_dihedrals)}, non_clash_dihedrals={non_clash_dihedrals}, "
                f"conf_analysis(normalized)={conf_norm!r}"
            )
            self.log.error(msg)
            raise RuntimeError(msg)

        # step 4: scanning + parameterization (+ plotting)
        self._banner("STEP 5: Dihedral scans & parameterization")
        scanner = DihedralScanner(self.cfg.base_dir, self.opt)
        for m in methods:
            calc, label = self.calc_factory.get(m)
            for d in dihedrals_to_scan:
                d_tup = tuple(d)

                # select structure for scan
                if use_conformers and d_tup in clash_dihedrals:
                    pdb_for_scan = best_conf[m].get(d_tup)
                    if not pdb_for_scan:
                        self.log.warning("Skipping %s for %s: no suitable conformer.", d_tup, m)
                        continue
                else:
                    pdb_for_scan = minimized[m]

                # run scan
                self.log.info("[%s][%s] Input for scan: %s", m, d_tup, pdb_for_scan)
                self.log.info("[%s][%s] SCAN start (%s)", m, d_tup, label)
                start = time.time()
                scanner.scan(pdb_for_scan, calc, d_tup, label)
                self.log.info("[%s][%s] SCAN done in %.2f s", m, d_tup, time.time() - start)
                # After scanning is complete, plotting
   
                plot_shifted_energy_profiles(self.cfg.base_dir, d_tup)

                self.log.info("[%s][%s] PLOT OK (shifted energy profile)", m, d_tup)

                # parameterization for this (method, dihedral)
                try:
                    self.parameterize_one_dihedral(ROOT, m, label, d_tup, pdb_for_scan)
                    self.log.info("[%s][%s] PARAM OK", m, d_tup)
                except Exception as e:
                    self.log.error("[%s][%s] PARAM FAILED: %s", m, d_tup, e)
                    continue

        # summary
        self._banner("STEP 6: Summary")
        scanned_dirs = []
        for d in dihedrals_to_scan:
            d_tup = tuple(d)
            dstr = "_".join(map(str, d_tup))
            if os.path.isdir(os.path.join(self.cfg.base_dir, "scanning", dstr)):
                scanned_dirs.append(dstr)
        self.log.info("Completed dihedral folders: %s", scanned_dirs if scanned_dirs else "—")

        self._sub("Final cumulative frcmods")
        if self.final_frcmod:
            for m, path in self.final_frcmod.items():
                self.log.info("[%s] %s", m, path)
        else:
            self.log.warning("No cumulative frcmod produced (no parameterizations ran).")

        self._banner("WORKFLOW END")
