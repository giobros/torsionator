import numpy as np
from ase import neighborlist
from .constants import LJ_PARAMETERS, VDW_RADII

def combine_lj(atom1: str, atom2: str):
    d = {'sigma': 3.4, 'epsilon': 0.359}
    s1 = LJ_PARAMETERS.get(atom1, d)['sigma']; e1 = LJ_PARAMETERS.get(atom1, d)['epsilon']
    s2 = LJ_PARAMETERS.get(atom2, d)['sigma']; e2 = LJ_PARAMETERS.get(atom2, d)['epsilon']
    return (s1 + s2) / 2.0, np.sqrt(e1 * e2)

def v_lj(r, sigma, epsilon):
    sr6 = (sigma / r) ** 6
    return 4 * epsilon * (sr6**2 - sr6)

def excluded_pairs(atoms):
    cut = neighborlist.natural_cutoffs(atoms)
    nl = neighborlist.NeighborList(cut, self_interaction=False, bothways=True)
    nl.update(atoms)
    bonded = {i: set(nl.get_neighbors(i)[0]) for i in range(len(atoms))}
    ex = set()
    # 1-2
    for i in bonded:
        for j in bonded[i]: ex.add(tuple(sorted((i, j))))
    # 1-3
    for i in bonded:
        for j in bonded[i]:
            for k in bonded[j]:
                if i != k: ex.add(tuple(sorted((i, k))))
    # 1-4
    for i in bonded:
        for j in bonded[i]:
            for k in bonded[j]:
                if k == i: continue
                for l in bonded[k]:
                    if l != j and l != i:
                        ex.add(tuple(sorted((i, l))))
    return ex

def lj_energy_between_groups(atoms, group_A, group_D, cutoff=12.0):
    pos = atoms.get_positions(); syms = atoms.get_chemical_symbols()
    ex = excluded_pairs(atoms)
    total = 0.0
    for i in group_A:
        for j in group_D:
            if i == j: continue
            ij = (i, j) if i < j else (j, i)
            if ij in ex: continue
            r = np.linalg.norm(pos[i]-pos[j])
            if r > cutoff: continue
            sigma, eps = combine_lj(syms[i], syms[j])
            total += v_lj(r, sigma, eps)
    return total

def local_rg(atoms, indices):
    pos = atoms.get_positions()[indices]
    masses = np.array([atoms[i].mass for i in indices])
    com = np.average(pos, axis=0, weights=masses)
    rg2 = np.sum(masses * np.sum((pos - com)**2, axis=1)) / np.sum(masses)
    return float(np.sqrt(rg2))
