#!/usr/bin/env python3
"""
Convert DeepH-predicted wavefunction coefficients to ABACUS WFC format.

Only gamma-point (k=0,0,0), non-spin-polarized case.

Usage
-----
    python deeph_to_abacus_wfc.py <data_dir> <output_path> [--all-bands]

Example
-------
    python deeph_to_abacus_wfc.py ./001 ./WFC_NAO_K1.txt --all-bands
"""

import numpy as np
from pathlib import Path
import click

from deepx_dock.compute.eigen.hamiltonian import HamiltonianObj


# -------------------------------------------------------------------
# Empirical per-l-shell mapping from DeepH (OpenMX convention) to
# ABACUS WFC ordering.
#
# For each l-shell with size = 2*l+1:
#   perm_inv[i] = DeepH orbital index to place at ABACUS WFC position i
#   sign[i]     = sign factor for ABACUS WFC position i
#
# Verified on Si DZP (64 atoms, 832 orbitals) with majority voting.
#
# DeepH (OpenMX) order per atom:
#   s: [s]
#   p: [p_x, p_y, p_z]  (m = +1, -1, 0)
#   d: [d_z², d_x²-y², d_xy, d_xz, d_yz]  (m = 0, +2, -2, +1, -1)
#
# ABACUS WFC order per atom (m = -l ... l, with empirical signs):
#   s: [s]
#   p: [-p_y, +p_z, +p_x]   (m = -1, 0, +1)
#   d: [-d_xy, +d_xz, +d_x²-y², -d_yz, -d_z²]
#      (m = -2, +1, +2, -1, 0)
# -------------------------------------------------------------------
_DEEPH_TO_WFC = {
    0: (np.array([0]),          np.array([-1.0])),
    1: (np.array([1, 2, 0]),    np.array([-1.0, 1.0, 1.0])),
    2: (np.array([2, 3, 1, 4, 0]), np.array([-1.0, 1.0, 1.0, -1.0, -1.0])),
}


def build_permutation(elements_orbital_map, elements):
    """
    Build global perm + sign arrays for all atoms.

    Parameters
    ----------
    elements_orbital_map : dict
        e.g. {"Si": [0, 0, 1, 1, 2]}
    elements : list of str
        Element symbol for each atom, e.g. ["Si", "Si", ...]

    Returns
    -------
    perm : ndarray of shape (n_orb,)
        perm[abacus_idx] = deeph_idx
    sign : ndarray of shape (n_orb,)
        sign factor for each ABACUS WFC position
    """
    atom_num_orbits = [
        sum(2 * l + 1 for l in elements_orbital_map[el]) for el in elements
    ]
    atom_cumsum = np.insert(np.cumsum(atom_num_orbits), 0, 0)
    n_orb = atom_cumsum[-1]

    perm = np.empty(n_orb, dtype=np.intp)
    sign = np.empty(n_orb, dtype=np.float64)

    for i_atom, el in enumerate(elements):
        orbs = elements_orbital_map[el]
        base = atom_cumsum[i_atom]
        deeph_off = 0
        abacus_off = 0

        for l_val in orbs:
            l_perm_inv, l_sign = _DEEPH_TO_WFC[l_val]
            n_m = 2 * l_val + 1
            for j in range(n_m):
                idx = base + abacus_off + j
                perm[idx] = base + deeph_off + l_perm_inv[j]
                sign[idx] = l_sign[j]
            deeph_off += n_m
            abacus_off += n_m

    return perm, sign


def apply_transform(eigvecs, perm, sign):
    """
    Apply permutation + sign to convert eigenvectors from
    DeepH ordering to ABACUS WFC ordering.

    Parameters
    ----------
    eigvecs : ndarray, shape (n_orb, n_bands)
    perm : ndarray, shape (n_orb,)
    sign : ndarray, shape (n_orb,)

    Returns
    -------
    wfc : ndarray, shape (n_orb, n_bands)
    """
    return sign[:, None] * eigvecs[perm, :]


def write_wfc(output_path, eigvals_ry, wfc_coeffs):
    """
    Write eigenvectors in ABACUS WFC_NAO_K*.txt plain-text format.

    Format per band:
        line 1:   <n_kpts> (index of k points)
        line 2:   kx ky kz
        line 3:   <n_bands> (number of bands)
        line 4:   <n_orb> (number of orbitals)
        line 5:   <iband> (band)
        line 6:   <energy_Ry> (Ry)
        line 7:   <occupation> (Occupations)
        line 8+:  <real> <imag> ... (5 pairs per line)
    """
    n_orb, n_bands = wfc_coeffs.shape

    with open(output_path, 'w') as f:
        f.write("1 (index of k points)\n")
        f.write("0 0 0\n")
        f.write(f"{n_bands} (number of bands)\n")
        f.write(f"{n_orb} (number of orbitals)\n")

        for b in range(n_bands):
            f.write(f"{b + 1} (band)\n")
            f.write(f"{eigvals_ry[b]:.25e} (Ry)\n")
            f.write("0.0000000000000000000000000e+00 (Occupations)\n")

            coeffs = wfc_coeffs[:, b]
            for i in range(0, n_orb, 5):
                chunk = coeffs[i:i + 5]
                parts = []
                for c in chunk:
                    parts.append(f"{c.real:.25e} 0.0000000000000000000000000e+00")
                f.write(" ".join(parts) + "\n")

    print(f"[done] {n_bands} band(s) written to {output_path}")


@click.command()
@click.argument("data_dir", type=click.Path(exists=True))
@click.argument("output_path", type=click.Path())
@click.option("--all-bands", is_flag=True, default=False,
              help="Export all bands (default: only first band)")
def main(data_dir, output_path, all_bands):
    """
    Read DeepH data from DATA_DIR, diagonalize at gamma, and
    write eigenvectors in ABACUS WFC_NAO_K*.txt format to OUTPUT_PATH.
    """
    data_dir = Path(data_dir)
    output_path = Path(output_path)

    # 1. Load Hamiltonian and diagonalize at gamma
    ham = HamiltonianObj(data_dir)
    eigvals, eigvecs = ham.diag(
        np.array([[0.0, 0.0, 0.0]]),
        n_jobs=-1,
        parallel_k=True,
        bands_only=False,
    )
    # eigvals: (n_bands, 1) --> gamma_eigvals: (n_bands,)
    # eigvecs: (n_orb, n_bands, 1) --> gamma_eigvecs: (n_orb, n_bands)
    gamma_eigvals = eigvals[:, 0].copy()
    gamma_eigvecs = eigvecs[:, :, 0].copy()

    # 2. Build global permutation from elements_orbital_map
    perm, sign = build_permutation(ham.elements_orbital_map, ham.elements)

    # 3. Apply transformation
    wfc_coeffs = apply_transform(gamma_eigvecs, perm, sign)

    # 4. Optionally keep only the first band
    if not all_bands:
        wfc_coeffs = wfc_coeffs[:, :1]
        gamma_eigvals = gamma_eigvals[:1]

    # 5. Write to file
    ev_to_ry = 1.0 / 13.605693122994
    write_wfc(output_path, gamma_eigvals * ev_to_ry, wfc_coeffs)


if __name__ == "__main__":
    main()
