import os
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from ase.io import read

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
    ClashAnalyzer,
    DihedralScanner,
    GraphSplitter,
    plot_shifted_energy_profiles,
    plot_shifted_energy_profiles_mcs,
    resolve_dihedrals_arg,
)
from .multi_conf_scan import (
    _MCS_merge_scans,
    _generate_and_minimize_conformers,
    _scan_all_conformers,
)
from .selection import ConformerSelector
from . import io_utils


@dataclass
class Workflow:
    cfg: Config
    final_frcmod: Dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.log = setup_logging(self.cfg.log_file)
        self.opt = GeometryOptimizer()
        self.calc_factory = CalculatorFactory(self.cfg.obiwan_path)
        self.clash = ClashAnalyzer()

    # -------------------------------------------------------------------------
    # Multiplicity → basis tuple
    # -------------------------------------------------------------------------

    @staticmethod
    def _basis_from_max_mult(n: int) -> Tuple[int, ...]:
        presets = {
            1: (1,),
            2: (-2, 1),
            3: (-3, -2, 1),
            4: (-4, -3, -2, 1),
            6: (-6, -4, -3, -2, 1),
        }
        if n not in presets:
            raise ValueError(f"multiplicity must be one of 0, 1, 2, 3, 4, 6 (got {n})")
        return presets[n]

    # -------------------------------------------------------------------------
    # Clash helpers
    # -------------------------------------------------------------------------

    def _pick_lowest_energy_clash_free(
        self,
        conf_method_dir: str,
        d_tup: Tuple[int, int, int, int],
        group_A: List[int],
        group_D: List[int],
        steps: int = 37,
    ) -> Optional[str]:
        for pdb_path in iter_energy_ordered_conformers(conf_method_dir):
            atoms = read(pdb_path)
            res = self.clash.scan_check(atoms, d_tup, group_A, group_D, steps=steps)
            if res.no_clash:
                return pdb_path
        return None

    def _detect_clashes(
        self,
        minimized: Dict[str, str],
        dihedrals_to_scan: List[Tuple[int, int, int, int]],
        methods: List[str],
    ) -> Dict[str, Dict[Tuple[int, int, int, int], bool]]:
        """Return {method: {d_tup: has_clash}} for all method/dihedral combinations."""
        splitter = GraphSplitter()
        clash_by_method: Dict[str, Dict[Tuple[int, int, int, int], bool]] = {m: {} for m in methods}

        for method in methods:
            atoms = read(minimized[method])
            for d_tup in dihedrals_to_scan:
                group_A, group_D = splitter.split_by_dihedral(atoms, d_tup)
                res = self.clash.scan_check(atoms, d_tup, group_A, group_D, steps=37)
                has_clash = not res.no_clash
                clash_by_method[method][d_tup] = has_clash
                self.log.info("[%s][%s] initial clash=%s", method, d_tup, has_clash)

        return clash_by_method

    def _handle_conformers(
        self,
        pdb_file: str,
        minimized: Dict[str, str],
        clash_by_method: Dict[str, dict],
        methods: List[str],
        RMSD: float,
        n_confs: int,
        conf_norm: Optional[str],
        bcs_norm: Optional[str],
    ):
        """
        Conformer selection logic (non-MCS path):

        - No clashes anywhere      → skip conformer pipeline.
        - Clashes exist:
            conf_analysis=None     → raise error explaining the options.
            conf_analysis='false'  → scan minimized input as-is.
            conf_analysis='true'   → generate & minimize conformers, then for each
                                     clashing (method, dihedral):
                                       1. pick lowest-energy clash-free conformer
                                       2. if none: use best_lj if BCS=true, else raise.

        Returns (use_conformers, best_conf).
        """
        # No clashes → nothing to do
        if all(not any(clash_by_method[m].values()) for m in methods):
            return False, {m: {} for m in methods}

        # Clashes exist but user didn't specify --conf_analysis
        if conf_norm is None:
            msg = (
                "Clashes detected in at least one target dihedral.\n\n"
                "Set --conf_analysis to:\n"
                "  true  → generate conformers and try to find a clash-free starting geometry\n"
                "  false → scan the minimized input even if clashes are present\n\n"
                "Clashes can distort the torsion scan profile."
            )
            self.log.error(msg)
            raise RuntimeError(msg)

        if conf_norm not in ("true", "false"):
            raise RuntimeError(f"Invalid --conf_analysis value: {conf_norm!r}. Use 'true' or 'false'.")

        # User explicitly accepts clashes
        if conf_norm == "false":
            self.log.info("[CONF] conf_analysis=false -> scanning minimized input (clashes accepted).")
            return False, {m: {} for m in methods}

        # conf_analysis='true': full conformer pipeline
        self.log.info("[CONF] conf_analysis=true -> generating/minimizing conformers.")

        gen = ConformerGenerator(os.path.join(self.cfg.base_dir, "conformers"))
        conf_dir = gen.generate(pdb_file, n_confs)

        screen = ConformerScreen()
        screen.screen_by_rmsd(conf_dir, RMSD)
        screen.ensure_ref_if_unique(conf_dir, pdb_file, RMSD)
        gen.normalize_conformer(conf_dir)

        minim = ConformerMinimizer(self.cfg.base_dir, self.opt)
        for m in methods:
            calc, _ = self.calc_factory.get(m)
            minim.minimize_folder(m, calc, conf_dir)

        selector = ConformerSelector(self.clash)
        splitter = GraphSplitter()
        best_conf: Dict[str, Dict[Tuple[int, int, int, int], str]] = {m: {} for m in methods}

        for m in methods:
            atoms_min = read(minimized[m])
            conf_method_dir = os.path.join(self.cfg.base_dir, "conformers", m)

            for d_tup in sorted(clash_by_method[m]):
                gA, gD = splitter.split_by_dihedral(atoms_min, d_tup)

                # Try lowest-energy clash-free conformer
                clash_free = self._pick_lowest_energy_clash_free(conf_method_dir, d_tup, gA, gD)
                if clash_free:
                    best_conf[m][d_tup] = clash_free
                    self.log.info("[CONF][%s][%s] using clash-free conformer: %s", m, d_tup, clash_free)
                    continue

                # Fallback: best_lj if BCS=true
                if bcs_norm != "true":
                    msg = (
                        f"[CONF][{m}][{d_tup}] No clash-free conformer found. "
                        "Set --BCS true to allow LJ-based fallback, or revise the input."
                    )
                    self.log.error(msg)
                    raise RuntimeError(msg)

                res = selector.select_best(conf_method_dir, d_tup, gA, gD, steps=37)
                if res.get("best_lj"):
                    best_conf[m][d_tup] = res["best_lj"]["path"]
                    self.log.warning(
                        "[CONF][%s][%s] no clash-free found; BCS=true -> best_lj: %s (LJ=%.4f)",
                        m, d_tup, best_conf[m][d_tup], res["best_lj"]["lj"],
                    )
                    if res.get("best_rg"):
                        self.log.info(
                            "[CONF][%s][%s] best_rg: %s (RG=%.4f)",
                            m, d_tup, res["best_rg"]["path"], res["best_rg"]["rg"],
                        )
                    continue

                msg = (
                    f"[CONF][{m}][{d_tup}] No clash-free conformer found and BCS fallback yielded no best_lj. "
                    "Revise the input conformers."
                )
                self.log.error(msg)
                raise RuntimeError(msg)

        return True, best_conf

    # -------------------------------------------------------------------------
    # Main entry
    # -------------------------------------------------------------------------

    def run(
        self,
        pdb_file: str,
        method: str,
        dihedral: str,
        conf_analysis: Optional[str] = "false",
        multiplicity: int = 6,
        BCS: Optional[str] = "false",
        RMSD: float = 0.5,
        n_confs: int = 20,
        MCS: Optional[str] = "true",
        double_rotation: Optional[str] = "false",
        step_size: int = 10,
        net_charge: int = 0,
        spin: int = 1,
    ) -> None:
        t0_total = time.time()
        self.calc_factory.net_charge = net_charge
        self.calc_factory.spin = spin
        root_name = os.path.splitext(os.path.basename(pdb_file))[0]

        conf_norm = io_utils.norm_flag(conf_analysis)
        bcs_norm = io_utils.norm_flag(BCS)
        mcs_norm = io_utils.norm_flag(MCS)
        dbl_norm = io_utils.norm_flag(double_rotation)
        step_size = io_utils.validate_step_size(step_size)

        if multiplicity == 0:
            auto_basis = True
            expand_basis = True
            basis_tuple = None
        else:
            auto_basis = False
            expand_basis = True
            basis_tuple = self._basis_from_max_mult(int(multiplicity))

        io_utils.banner(self.log, f"WORKFLOW START for {root_name}")
        self.log.warning(
            "WARNING: Since during the minimization, the protonation state and connectivity "
            "of the molecule can change, we suggest looking at the antechamber sqm.pdb output "
            "and NNP geometries.xyz one to check if any the molecule underwent any changes "
            "during the pipeline execution."
        )
        methods = [method]

        if net_charge != 0:
            charge_incompatible = [m for m in methods if m in ("obi", "mace")]
            if charge_incompatible:
                msg = (
                    f"Method(s) {charge_incompatible} do not support charged molecules "
                    f"(net_charge={net_charge}). Use --method uma for charged systems."
                )
                self.log.error(msg)
                raise RuntimeError(msg)

        if spin != 1:
            spin_incompatible = [m for m in methods if m != "uma"]
            if spin_incompatible:
                msg = (
                    f"Method(s) {spin_incompatible} do not support non-singlet spin states "
                    f"(spin={spin}). Use --method uma for open-shell systems."
                )
                self.log.error(msg)
                raise RuntimeError(msg)

        dihedrals, should_exit = resolve_dihedrals_arg(self.log, dihedral, pdb_file)
        if should_exit:
            return

        dihedrals_to_scan = sorted({tuple(d) for d in dihedrals})

        # Validate user-provided dihedral indices
        if dihedral.strip().startswith("[") and dihedral.strip().endswith("]"):
            for d_tup in dihedrals_to_scan:
                io_utils.validate_user_dihedral_indices(pdb_file, d_tup)

        if mcs_norm == "true":
            self._run_mcs(
                pdb_file=pdb_file,
                methods=methods,
                dihedrals=dihedrals_to_scan,
                root_name=root_name,
                n_confs=n_confs,
                RMSD=RMSD,
                basis_tuple=basis_tuple,
                expand_basis=expand_basis,
                auto_basis=auto_basis,
                double_rotation=(dbl_norm == "true"),
                step_size=step_size,
                net_charge=net_charge,
                spin=spin,
            )
            io_utils.banner(self.log, "WORKFLOW END")
            return

        if dbl_norm == "true":
            self.log.warning("--double_rotation ignored because MCS=false")

        minimized = io_utils.minimize_input(
            self.cfg, self.opt, self.calc_factory, self.log, pdb_file, methods
        )
        clash_by_method = self._detect_clashes(minimized, dihedrals_to_scan, methods)
        use_conformers, best_conf = self._handle_conformers(
            pdb_file=pdb_file,
            minimized=minimized,
            clash_by_method=clash_by_method,
            methods=methods,
            RMSD=RMSD,
            n_confs=n_confs,
            conf_norm=conf_norm,
            bcs_norm=bcs_norm,
        )

        scanner = DihedralScanner(self.cfg.base_dir, self.opt)

        for m in methods:
            calc, label = self.calc_factory.get(m)

            for d_tup in dihedrals_to_scan:
                pdb_for_scan = best_conf[m].get(d_tup, minimized[m]) if use_conformers else minimized[m]

                io_utils.run_nn_scan(self.log, scanner, pdb_for_scan, calc, d_tup, label, m, step_size=step_size)

                io_utils.run_gaff2_old(
                    self.log, self.cfg, root_name, m, d_tup, pdb_for_scan,
                    basis_periodicities=basis_tuple or (1,),
                    postscan_make_expanded_for_mdgx=not auto_basis,
                    MCS=False, mcs_xyz=None, step_size=step_size,
                    net_charge=net_charge, spin=spin,
                )

                effective_basis = (
                    io_utils.detect_basis_from_antechamber(self.log, self.cfg, d_tup, m)
                    if auto_basis else basis_tuple
                )

                param_results = io_utils.run_parametrization(
                    self.log, self.cfg, m, d_tup,
                    MCS=False, basis_periodicities=effective_basis, expand_basis=expand_basis,
                )

                io_utils.run_gaff2_new(
                    self.log, self.cfg, root_name, m, d_tup, pdb_for_scan,
                    param_results, MCS=False, mcs_xyz=None, step_size=step_size,
                    net_charge=net_charge, spin=spin,
                )

                best_profile = io_utils.select_and_update_frcmod(
                    self.log, self.cfg, root_name, m, d_tup,
                    self.final_frcmod, param_results, tag="", MCS=False,
                    expand_basis=expand_basis,
                )

                try:
                    plot_shifted_energy_profiles(self.cfg.base_dir, d_tup, m, best_profile)
                except Exception as exc:
                    self.log.warning("[%s][%s] plot failed: %s", m, d_tup, exc)

        io_utils.banner(self.log, "WORKFLOW END")

    # -------------------------------------------------------------------------
    # MCS workflow
    # -------------------------------------------------------------------------

    def _run_mcs(
        self,
        *,
        pdb_file: str,
        methods: List[str],
        dihedrals: List[Tuple[int, int, int, int]],
        root_name: str,
        n_confs: int,
        RMSD: float,
        basis_tuple: Optional[Tuple[int, ...]],
        expand_basis: bool,
        auto_basis: bool = False,
        double_rotation: bool,
        step_size: int,
        net_charge: int = 0,
        spin: int = 1,
    ) -> None:
        io_utils.banner(self.log, "RUNNING MCS WORKFLOW")

        per_method_conf = _generate_and_minimize_conformers(self, pdb_file, methods, n_confs, RMSD)
        scanned_tags = _scan_all_conformers(
            self, per_method_conf, methods, dihedrals,
            step_size=step_size, double_rotation=double_rotation,
        )

        for m in methods:
            ref_pdb = os.path.join(self.cfg.base_dir, "conformers", m, "0", "minimized.pdb")

            for d_tup in dihedrals:
                out_dir = _MCS_merge_scans(
                    self,
                    self.cfg.base_dir,
                    d_tup,
                    m,
                    scanned_tags[m].get(d_tup, []),
                    step_size=step_size,
                    mode="absolute",
                )

                mcs_xyz = os.path.join(out_dir, "geometries.xyz")
                if not os.path.exists(mcs_xyz):
                    raise FileNotFoundError(f"[MCS][{m}][{d_tup}] missing merged xyz in {out_dir}")

                io_utils.run_gaff2_old(
                    self.log, self.cfg, root_name, m, d_tup, ref_pdb,
                    basis_periodicities=basis_tuple or (1,),
                    postscan_make_expanded_for_mdgx=not auto_basis,
                    MCS=True, mcs_xyz=mcs_xyz, step_size=step_size,
                    net_charge=net_charge, spin=spin,
                )

                effective_basis = (
                    io_utils.detect_basis_from_antechamber(self.log, self.cfg, d_tup, m)
                    if auto_basis else basis_tuple
                )

                param_results = io_utils.run_parametrization(
                    self.log, self.cfg, m, d_tup,
                    MCS=True, basis_periodicities=effective_basis, expand_basis=expand_basis,
                )

                io_utils.run_gaff2_new(
                    self.log, self.cfg, root_name, m, d_tup, ref_pdb,
                    param_results, MCS=True, mcs_xyz=mcs_xyz, step_size=step_size,
                    net_charge=net_charge, spin=spin,
                )

                best_profile = io_utils.select_and_update_frcmod(
                    self.log, self.cfg, root_name, m, d_tup,
                    self.final_frcmod, param_results, tag="_MCS", MCS=True,
                    expand_basis=expand_basis,
                )

                try:
                    plot_shifted_energy_profiles_mcs(
                        self.cfg.base_dir, d_tup, m, best_profile=best_profile,
                    )
                except Exception as exc:
                    self.log.warning("[MCS][%s][%s] plot failed: %s", m, d_tup, exc)
