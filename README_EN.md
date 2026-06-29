[**English**](README_EN.md) | [**中文**](README.md)

# DeepH → ABACUS WFC Converter

Convert eigenvectors from DeepH-format `hamiltonian.h5` to ABACUS `WFC_NAO_K*.txt`, then use ABACUS `get_wf` to obtain real-space `.cube` wavefunction files.

> Example: 64-atom Si with DZP basis (832 orbitals) at Gamma point.

---

## Dependencies

| Tool | Purpose |
|------|---------|
| DeepH-dock (deepx-dock) | Diagonalize `hamiltonian.h5` to obtain eigenvalues / eigenvectors |
| ABACUS = 3.10-LTS | `calculation = get_wf` to read WFC and produce `.cube` files |

---

## Data Preparation

The DeepH-format input directory (`001/` in this repository):

```
001/
├── hamiltonian.h5    # Hamiltonian matrix (real space)
├── overlap.h5        # Overlap matrix
├── info.json         # Metadata: #atoms, #orbitals, basis set, Fermi level
└── POSCAR            # Lattice + atomic coordinates
```

The `info.json` shipped in this repo:

```json
{"atoms_quantity": 64, "orbits_quantity": 832, "orthogonal_basis": false,
 "spinful": false, "fermi_energy_eV": 6.4157182229,
 "elements_orbital_map": {"Si": [0, 0, 1, 1, 2]}}
```

This means Si has 5 zeta groups (2s + 2p + 1d), and each Si atom contributes 2×1 + 2×3 + 1×5 = 13 orbitals.

---

## Usage

```bash
python deeph_to_abacus_wfc.py ./001 ./WFC_NAO_K1.txt              # first band only (default)
python deeph_to_abacus_wfc.py ./001 ./WFC_NAO_K1.txt --all-bands  # all bands
```

### Arguments

| Argument | Description |
|----------|-------------|
| `data_dir` | DeepH data directory (containing `hamiltonian.h5`, etc.) |
| `output_path` | Output WFC file path |
| `--all-bands` | Export all bands (default: first band only) |

### Output Format

Plain-text WFC_NAO_K*.txt, one block per band:

```
1 (index of k points)
0 0 0
154 (number of bands)
832 (number of orbitals)
1 (band)
-4.2796723504902273171879301e-01 (Ry)
0.0000000000000000000000000e+00 (Occupations)
<real_1> <imag_1> <real_2> <imag_2> ...   # 5 pairs per line, n_orb pairs total
```

Orbital coefficient order: **atom-by-atom**; within each atom:  
**s₁ s₂ p₁(m=-1,0,+1) p₂(m=-1,0,+1) d₁(m=-2,+1,+2,-1,0)**.

---

## Orbital Convention Mapping

DeepH (OpenMX convention) and ABACUS WFC differ in basis-function ordering; a **permutation + sign** transform is required.

### Per-atom ordering comparison

| l | DeepH (OpenMX) | m order | ABACUS WFC | m order |
|---|----------------|---------|------------|---------|
| 0 | `[s]` | `[0]` | `[-s]` | `[0]` |
| 1 | `[p_x, p_y, p_z]` | `[+1, -1, 0]` | `[-p_y, +p_z, +p_x]` | `[-1, 0, +1]` |
| 2 | `[d_z², d_x²-y², d_xy, d_xz, d_yz]` | `[0, +2, -2, +1, -1]` | `[-d_xy, +d_xz, +d_x²-y², -d_yz, -d_z²]` | `[-2, +1, +2, -1, 0]` |

### Per-l-shell mapping table

DeepH → ABACUS WFC: `wfc[abacus_pos] = sign[j] × deeph[perm_inv[j]]`

| l | `perm_inv` (deeph idx → abacus pos j) | `sign` |
|---|----------------------------------------|--------|
| 0 | `[0]` | `[-1]` |
| 1 | `[1, 2, 0]` | `[-1, +1, +1]` |
| 2 | `[2, 3, 1, 4, 0]` | `[-1, +1, +1, -1, -1]` |

> This mapping was validated on the Si 64-atom DZP system via majority voting: **s: 64/64, p: 60+/64, d: 64/64** atoms agree.

---

## Workflow: DeepH → Real-space Wavefunction

```mermaid
flowchart LR
    A[DeepH: hamiltonian.h5] --> B[deeph_to_abacus_wfc.py]
    B --> C[WFC_NAO_K1.txt]
    C --> D[ABACUS get_wf]
    D --> E[.cube real-space wavefunction]
```

### Step-by-step

1. **Generate WFC file**

   ```bash
   python deeph_to_abacus_wfc.py ./001 ./WFC_NAO_K1.txt --all-bands
   ```

2. **Copy into ABACUS output directory**

   ```bash
   cp WFC_NAO_K1.txt /path/to/abacus/OUT.{}/
   ```

3. **Configure ABACUS INPUT**

   ```
   calculation     get_wf
   out_wfc_norm    1          # band to output (see ABACUS manual)
   ```

   The `KPT` file should contain only the Gamma point `0 0 0`.

4. **Run ABACUS**

   ```bash
   mpirun -np 4 abacus
   ```

   Produces `BAND1_k_1_s_1_ENV.cube` and similar Cube files.

---

## Validation

> The `hamiltonian.h5` in `001/` was converted **directly from an ABACUS DFT result**, not from a DeepH neural-network prediction. Therefore this validates the **correctness of the conversion code** itself, not the accuracy of DeepH.

### Compared pipelines

| Pipeline | Description |
|----------|-------------|
| **A (Baseline)** | ABACUS DFT → `calculation get_wf` → `.cube` |
| **B (This tool)** | ABACUS DFT → DeepH-dock `ham.diag()` → `deeph_to_abacus_wfc.py` → ABACUS `get_wf` → `.cube` |

### Results (isosurface comparison)

| ABACUS DFT direct (Baseline) | After conversion by this tool |
|:---:|:---:|
| ![](./figs/abacus.png) | ![](./figs/deeph.png) |

Both figures show the isosurface of band 1 (VBM) at the Gamma point. The agreement confirms that the orbital permutation + sign mapping is correct.

---

## Limitations

- **Gamma point only**: only outputs `(0, 0, 0)` k-point
- **Non-spin-polarized**: `nspin = 2` is not supported
- **Mapping validated scope**: this mapping was validated on Si DZP; other elements / basis sets may need re-validation
  - `l = 0, 1, 2` are covered; `l >= 3` is not handled
- **Occupation**: written as 0 (placeholder), not read from `info.json` or `wavefunction_ao.h5`

-----
ABACUS get `WFC_NAO_K1.txt`
```
calculation     scf
out_wfc_lcao    1
```