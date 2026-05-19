import numpy as np
from ase import neighborlist

from .constants import LJ_PARAMETERS, VDW_RADII

_LJ_DEFAULT = {"sigma": 3.4, "epsilon": 0.359}


def _lj_params(symbol: str):
    p = LJ_PARAMETERS.get(symbol, _LJ_DEFAULT)
    return p["sigma"], p["epsilon"]


def _lorentz_berthelot(sym1: str, sym2: str):
    """Lorentz-Berthelot combining rules: arithmetic sigma, geometric epsilon."""
    s1, e1 = _lj_params(sym1)
    s2, e2 = _lj_params(sym2)
    return (s1 + s2) / 2.0, np.sqrt(e1 * e2)


def _v_lj(r: float, sigma: float, epsilon: float) -> float:
    sr6 = (sigma / r) ** 6
    return 4.0 * epsilon * (sr6 ** 2 - sr6)


def _excluded_pairs(atoms) -> set:
    """Return set of 1-2, 1-3, and 1-4 bonded pairs (sorted tuples)."""
    cutoffs = neighborlist.natural_cutoffs(atoms)
    nl = neighborlist.NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)

    bonded = {i: set(nl.get_neighbors(i)[0]) for i in range(len(atoms))}
    excluded = set()

    for i, neighbors_i in bonded.items():
        # 1-2
        for j in neighbors_i:
            excluded.add((min(i, j), max(i, j)))
        # 1-3
        for j in neighbors_i:
            for k in bonded[j]:
                if k != i:
                    excluded.add((min(i, k), max(i, k)))
        # 1-4
        for j in neighbors_i:
            for k in bonded[j]:
                if k == i:
                    continue
                for m in bonded[k]:
                    if m != j and m != i:
                        excluded.add((min(i, m), max(i, m)))

    return excluded


def lj_energy_between_groups(atoms, group_A, group_D, cutoff: float = 12.0) -> float:
    """Non-bonded LJ energy between two disjoint atom groups (1-4 excluded)."""
    pos = atoms.get_positions()
    syms = atoms.get_chemical_symbols()
    excluded = _excluded_pairs(atoms)

    total = 0.0
    for i in group_A:
        for j in group_D:
            if i == j:
                continue
            pair = (min(i, j), max(i, j))
            if pair in excluded:
                continue
            r = np.linalg.norm(pos[i] - pos[j])
            if r > cutoff:
                continue
            sigma, eps = _lorentz_berthelot(syms[i], syms[j])
            total += _v_lj(r, sigma, eps)

    return total


def local_rg(atoms, indices) -> float:
    """Mass-weighted radius of gyration for a subset of atoms."""
    pos = atoms.get_positions()[indices]
    masses = np.array([atoms.get_masses()[i] for i in indices])
    com = np.average(pos, axis=0, weights=masses)
    rg2 = np.sum(masses * np.sum((pos - com) ** 2, axis=1)) / np.sum(masses)
    return float(np.sqrt(rg2))
