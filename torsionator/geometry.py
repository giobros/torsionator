import numpy as np
from ase.optimize import BFGS
from ase.constraints import FixInternals
from .io_utils import write_pdb
from .constants import CONV_EV_TO_EH

class GeometryOptimizer:
    def __init__(self, fmax=1e-4, steps=1000):
        self.fmax = fmax
        self.steps = steps

    @staticmethod
    def gradient_descent(atoms, learning_rate=0.01, max_steps=100, force_tol=0.01):
        for _ in range(max_steps):
            forces = atoms.get_forces()
            if np.max(np.abs(forces)) < force_tol:
                break
            atoms.positions += learning_rate * forces
        return atoms

    def minimize(self, atoms):
        self.gradient_descent(atoms)
        opt = BFGS(atoms, trajectory=None,logfile=None)
        opt.run(fmax=self.fmax, steps=self.steps)
        return atoms

    def minimize_and_write_pdb(self, pdb_input: str, calc, out_path: str):
        from ase.io import read
        atoms = read(pdb_input)
        atoms.calc = calc
        self.minimize(atoms)
        write_pdb(out_path, atoms)
        return out_path

class DihedralStepper:
    """Utility to optimize with a dihedral fixed at successive angles."""
    def __init__(self, fmax=1e-4, steps=1000):
        self.fmax, self.steps = fmax, steps

    def optimize_with_dihedral(self, atoms, dihedral_indices, angle_deg):
        from ase.optimize import BFGS
        atoms.set_dihedral(*dihedral_indices, angle_deg)
        atoms.set_constraint(FixInternals(dihedrals_deg=[[angle_deg, dihedral_indices]]))
        opt = BFGS(atoms, trajectory=None,logfile=None)
        opt.run(fmax=self.fmax, steps=self.steps)
        return atoms

def energy_eh(atoms):
    return atoms.get_potential_energy() * CONV_EV_TO_EH
