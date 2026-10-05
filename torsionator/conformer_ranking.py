"""
MCS + all-dihedrals conformer ranking.

After the MCS merge, every scanned dihedral has a minimum-energy profile
(scanning/<a_b_c_d>/<method>/MCS/). This module takes the lowest-energy point
of each profile, collects those geometries into one conformer set and ranks
them by energy.

Relative energies are computed from the *absolute* NNP energies (Hartree)
written by _MCS_merge_scans, not from the per-dihedral shifted files
(energies.dat / angles_vs_energies_final.txt), whose minimum is always 0.
All energies come from the same calculator, molecule, charge and spin, so
their differences are directly comparable across dihedrals.

Note: each geometry is a constrained minimum (the scanned dihedral is fixed
on the step_size grid), not a fully relaxed stationary point.

Output directory: conformer_ranking/<method>/
  - ranking.txt                 all per-dihedral minima, ranked by energy
  - conformers_ranked.xyz       unique conformers (symmetry-aware heavy-atom RMSD), lowest energy first
  - conformers_ranked_all.xyz   every per-dihedral minimum, lowest energy first
"""

import os
from typing import List, NamedTuple, Optional, Tuple

import numpy as np
from rdkit import Chem

from .constants import CONV_EH_TO_KCAL_MOL
from .io_utils import banner


class _Minimum(NamedTuple):
    dihedral: Tuple[int, int, int, int]
    angle: float
    energy_eh: float
    source: str
    symbols: List[str]
    positions: np.ndarray


def _read_xyz_frames(path: str) -> List[Tuple[str, List[str], np.ndarray]]:
    """Plain multi-frame xyz reader returning (comment, symbols, positions)."""
    frames = []
    with open(path, encoding="utf-8") as f:
        lines = f.read().splitlines()
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i].split()[0])
        comment = lines[i + 1]
        symbols, coords = [], []
        for line in lines[i + 2:i + 2 + n]:
            parts = line.split()
            symbols.append(parts[0])
            coords.append([float(x) for x in parts[1:4]])
        frames.append((comment, symbols, np.array(coords)))
        i += 2 + n
    return frames


def _parse_header(comment: str) -> dict:
    out = {}
    for tok in comment.split():
        if "=" in tok:
            k, v = tok.split("=", 1)
            out[k] = v
    return out


def _kabsch_rmsd(p: np.ndarray, q: np.ndarray) -> float:
    p = p - p.mean(axis=0)
    q = q - q.mean(axis=0)
    h = p.T @ q
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    diff = (r @ p.T).T - q
    return float(np.sqrt((diff ** 2).sum() / len(p)))


def _heavy_atom_permutations(
    log, ref_pdb: Optional[str], symbols: List[str], max_matches: int,
) -> Tuple[np.ndarray, List[np.ndarray]]:
    """
    Heavy-atom indices and their symmetry-equivalent orderings (graph
    automorphisms of the heavy-atom skeleton, from the connectivity of
    `ref_pdb`). Each ordering is an array of full-molecule atom indices, to be
    compared against the heavy-atom indices in their original order.
    Falls back to the identity ordering.
    """
    heavy_idx = np.array([i for i, s in enumerate(symbols) if s.upper() not in ("H", "D")])
    identity = [heavy_idx]
    if not ref_pdb or not os.path.exists(ref_pdb):
        log.warning("[RANK] no reference PDB; RMSD without symmetry handling")
        return heavy_idx, identity
    mol = Chem.MolFromPDBFile(ref_pdb, removeHs=False, sanitize=False, proximityBonding=False)
    if mol is None or [a.GetSymbol().upper() for a in mol.GetAtoms()] != [s.upper() for s in symbols]:
        log.warning("[RANK] cannot use %s for symmetry; RMSD without symmetry handling", ref_pdb)
        return heavy_idx, identity
    if mol.GetNumBonds() == 0:
        # No CONECT records: every same-element atom would look equivalent.
        log.warning("[RANK] %s has no bonds; RMSD without symmetry handling", ref_pdb)
        return heavy_idx, identity

    # Heavy-atom skeleton; element + connectivity only (PDB bond orders are unreliable).
    # Removing atoms keeps the relative order, so skeleton index k <-> heavy_idx[k].
    query = Chem.RWMol(mol)
    for i in sorted((a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 1), reverse=True):
        query.RemoveAtom(int(i))
    for b in query.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
    query = query.GetMol()
    query.UpdatePropertyCache(strict=False)
    matches = query.GetSubstructMatches(
        query, uniquify=False, useChirality=False, maxMatches=max_matches,
    )
    if not matches:
        return heavy_idx, identity
    if len(matches) >= max_matches:
        log.warning("[RANK] symmetry permutations capped at %d", max_matches)
    return heavy_idx, [heavy_idx[np.array(m)] for m in matches]


def _symmetric_rmsd(
    p: np.ndarray, q: np.ndarray, heavy_idx: np.ndarray, perms: List[np.ndarray],
) -> float:
    """Lowest aligned heavy-atom RMSD over all symmetry-equivalent orderings of p."""
    q_heavy = q[heavy_idx]
    return min(_kabsch_rmsd(p[perm], q_heavy) for perm in perms)


def _load_profile_minimum(log, mcs_dir: str, d_tup) -> Optional[_Minimum]:
    xyz_path = os.path.join(mcs_dir, "geometries.xyz")
    if not os.path.exists(xyz_path):
        log.warning("[RANK][%s] missing %s", d_tup, xyz_path)
        return None
    frames = _read_xyz_frames(xyz_path)
    if not frames:
        log.warning("[RANK][%s] empty %s", d_tup, xyz_path)
        return None

    # Energies are stored in each frame header by _MCS_merge_scans (absolute, Eh).
    best = None
    for comment, symbols, pos in frames:
        hdr = _parse_header(comment)
        if "E" not in hdr:
            log.warning("[RANK][%s] frame without energy in %s", d_tup, xyz_path)
            return None
        e = float(hdr["E"])
        if best is None or e < best.energy_eh:
            best = _Minimum(
                dihedral=tuple(d_tup),
                angle=float(hdr.get("angle", "nan")),
                energy_eh=e,
                source=hdr.get("conformer", "?"),
                symbols=symbols,
                positions=pos,
            )
    return best


def _write_frame(fh, m: _Minimum, rank: int, rel_kcal: float) -> None:
    dih = "_".join(map(str, m.dihedral))
    fh.write(f"{len(m.symbols)}\n")
    fh.write(
        f"rank={rank} E={m.energy_eh:.8f} dE_kcal={rel_kcal:.4f} "
        f"dihedral={dih} angle={m.angle:.2f} conformer={m.source}\n"
    )
    for s, (x, y, z) in zip(m.symbols, m.positions):
        fh.write(f"{s} {x:.10f} {y:.10f} {z:.10f}\n")


def rank_mcs_profile_minima(
    log,
    base_dir: str,
    method: str,
    dihedrals: List[Tuple[int, int, int, int]],
    rmsd_threshold: float = 0.25,
    ref_pdb: Optional[str] = None,
    max_symmetry_matches: int = 1000,
) -> Optional[str]:
    """
    Collect the minimum of every MCS dihedral profile and rank them by energy.

    Minima whose heavy-atom RMSD (after Kabsch alignment) to a lower-energy one is
    below `rmsd_threshold` Å are flagged as duplicates. The RMSD is minimised over
    the symmetry-equivalent atom orderings derived from the connectivity of
    `ref_pdb` (e.g. phenyl flips, equivalent carboxylate oxygens); these are and excluded from
    conformers_ranked.xyz (they are kept in ranking.txt and conformers_ranked_all.xyz).

    Returns the output directory, or None if nothing could be ranked.
    """
    banner(log, f"MCS conformer ranking ({method})")

    minima: List[_Minimum] = []
    for d_tup in dihedrals:
        dih_str = "_".join(map(str, d_tup))
        mcs_dir = os.path.join(base_dir, "scanning", dih_str, method, "MCS")
        m = _load_profile_minimum(log, mcs_dir, d_tup)
        if m is not None:
            minima.append(m)

    if not minima:
        log.warning("[RANK][%s] no MCS profiles found; skipping ranking", method)
        return None

    minima.sort(key=lambda m: m.energy_eh)
    e0 = minima[0].energy_eh

    # Duplicate detection: compare each minimum with the unique ones below it.
    heavy_idx, perms = _heavy_atom_permutations(log, ref_pdb, minima[0].symbols, max_symmetry_matches)
    log.info("[RANK][%s] %d symmetry-equivalent atom orderings", method, len(perms))
    unique_idx: List[int] = []
    duplicate_of: List[Optional[int]] = []
    for i, m in enumerate(minima):
        dup = None
        for j in unique_idx:
            if _symmetric_rmsd(m.positions, minima[j].positions, heavy_idx, perms) < rmsd_threshold:
                dup = j
                break
        duplicate_of.append(dup)
        if dup is None:
            unique_idx.append(i)

    out_dir = os.path.join(base_dir, "conformer_ranking", method)
    os.makedirs(out_dir, exist_ok=True)

    unique_rank = {idx: r + 1 for r, idx in enumerate(unique_idx)}

    with open(os.path.join(out_dir, "ranking.txt"), "w", encoding="utf-8") as f:
        f.write(
            f"# Minima of MCS dihedral profiles ranked by absolute {method} energy\n"
            f"# dE relative to the global minimum; duplicate = symmetry-aware heavy-atom RMSD < "
            f"{rmsd_threshold} A (after alignment) to a lower-energy unique conformer\n"
        )
        f.write(
            f"{'rank':>4} {'unique':>6} {'dihedral':>16} {'angle':>8} "
            f"{'E_Eh':>18} {'dE_kcal/mol':>12} {'conformer':>16} {'duplicate_of':>12}\n"
        )
        for i, m in enumerate(minima):
            rel = (m.energy_eh - e0) * CONV_EH_TO_KCAL_MOL
            u = str(unique_rank[i]) if i in unique_rank else "-"
            dup = str(unique_rank[duplicate_of[i]]) if duplicate_of[i] is not None else "-"
            f.write(
                f"{i + 1:>4} {u:>6} {'_'.join(map(str, m.dihedral)):>16} {m.angle:>8.2f} "
                f"{m.energy_eh:>18.8f} {rel:>12.4f} {m.source:>16} {dup:>12}\n"
            )

    with open(os.path.join(out_dir, "conformers_ranked_all.xyz"), "w", encoding="utf-8") as f:
        for i, m in enumerate(minima):
            _write_frame(f, m, i + 1, (m.energy_eh - e0) * CONV_EH_TO_KCAL_MOL)

    with open(os.path.join(out_dir, "conformers_ranked.xyz"), "w", encoding="utf-8") as f:
        for i in unique_idx:
            m = minima[i]
            _write_frame(f, m, unique_rank[i], (m.energy_eh - e0) * CONV_EH_TO_KCAL_MOL)

    log.info(
        "[RANK][%s] %d profile minima, %d unique conformers -> %s",
        method, len(minima), len(unique_idx), os.path.relpath(out_dir, base_dir),
    )
    for i in unique_idx:
        m = minima[i]
        log.info(
            "[RANK][%s] #%d dE=%.3f kcal/mol dihedral=%s angle=%.1f conformer=%s",
            method, unique_rank[i], (m.energy_eh - e0) * CONV_EH_TO_KCAL_MOL,
            m.dihedral, m.angle, m.source,
        )
    return out_dir
