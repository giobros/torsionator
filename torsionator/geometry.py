import numpy as np
from ase.optimize import BFGS
from ase.constraints import FixInternals
from ase.io import read

from .calculators import attach_calc
from .io_utils import write_pdb
from .constants import CONV_EV_TO_EH


def energy_eh(atoms) -> float:
    return atoms.get_potential_energy() * CONV_EV_TO_EH


class GeometryOptimizer:
    def __init__(self, fmax: float = 1e-4, steps: int = 1000):
        self.fmax = fmax
        self.steps = steps

    def minimize(self, atoms):
        """Gradient descent pre-relaxation followed by BFGS."""
        self._gradient_descent(atoms)
        opt = BFGS(atoms, trajectory=None, logfile=None)
        opt.run(fmax=self.fmax, steps=self.steps)
        return atoms

    def minimize_and_write_pdb(self, pdb_input: str, calc, out_path: str) -> str:
        atoms = read(pdb_input)
        attach_calc(atoms, calc)
        self.minimize(atoms)
        write_pdb(out_path, atoms)
        return out_path

    @staticmethod
    def _gradient_descent(atoms, learning_rate: float = 0.01, max_steps: int = 100, force_tol: float = 0.01) -> None:
        for _ in range(max_steps):
            forces = atoms.get_forces()
            if np.max(np.abs(forces)) < force_tol:
                break
            atoms.positions += learning_rate * forces


class DihedralStepper:
    """Optimize a structure with the dihedral fixed at successive target angles."""

    def __init__(self, fmax: float = 1e-4, steps: int = 1000):
        self.fmax = fmax
        self.steps = steps

    def optimize_with_dihedral(self, atoms, dihedral_indices, angle_deg):
        atoms.set_dihedral(*dihedral_indices, angle_deg)
        atoms.set_constraint(FixInternals(dihedrals_deg=[[angle_deg, dihedral_indices]]))
        opt = BFGS(atoms, trajectory=None, logfile=None)
        opt.run(fmax=self.fmax, steps=self.steps)
        return atoms
