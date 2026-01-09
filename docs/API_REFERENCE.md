# HECTOR API Reference

Complete API documentation for the HECTOR source code modules.

---

## Table of Contents

1. [hector_mapper.py](#hector_mapperpy)
2. [maps_analysis_pairs_vs_all.py](#maps_analysis_pairs_vs_allpy)
3. [aln_fltr_4_dots.py](#aln_fltr_4_dotspy)
4. [sim.pyx](#simpyx)
5. [Data Structures](#data-structures)

---

## hector_mapper.py

Core fingerprinting/mapping engine for the HECTOR algorithm.

**Location:** `code/hector_mapper.py`

### Functions

---

#### `parse_ply(ply_fn, vrtx_cols=6, fc_cols=7)`

Parse PLY (Polygon File Format) files generated from EDTSurf or MeshLab.

**Parameters:**
| Name | Type | Default | Description |
|------|------|---------|-------------|
| `ply_fn` | str | - | Input PLY file name |
| `vrtx_cols` | int | 6 | Number of vertex information columns |
| `fc_cols` | int | 7 | Number of face information columns |

**Returns:**
| Name | Type | Description |
|------|------|-------------|
| `vrtx_arr` | numpy.ndarray | Array of floats (N_vertices x vrtx_cols) |
| `fc_arr` | numpy.ndarray | Array of ints (N_faces x fc_cols) |

**Example:**
```python
vrtx_arr, fc_arr = parse_ply("protein.ply")
coordinates = vrtx_arr[:, :3]  # x, y, z
normals = vrtx_arr[:, 3:6]     # nx, ny, nz
```

---

#### `normalise(in_arr)`

Normalize vectors to unit length.

**Parameters:**
| Name | Type | Description |
|------|------|-------------|
| `in_arr` | numpy.ndarray | Matrix of surface normal vectors (N x 3) |

**Returns:**
| Name | Type | Description |
|------|------|-------------|
| `norm_arr` | numpy.ndarray | Matrix of unit surface normal vectors (N x 3) |

**Example:**
```python
unit_normals = normalise(raw_normals)
# Each row now has magnitude 1.0
```

---

#### `spn_map_srfc(srfc_vrtx_arr, srfc_nrml_arr, spprt_dstnc_x, spprt_dstnc_y, bin_sz, dot_skp_frq, dot_map_frq)`

Core fingerprinting function. Generates spin image maps from surface vertices.

**Parameters:**
| Name | Type | Description |
|------|------|-------------|
| `srfc_vrtx_arr` | numpy.ndarray | Matrix of input vertex position vectors (N x 3) |
| `srfc_nrml_arr` | numpy.ndarray | Matrix of input vertex surface normal vectors (N x 3) |
| `spprt_dstnc_x` | float | Support distance alpha; cylindrical basis radius (Angstroms) |
| `spprt_dstnc_y` | float | Support distance beta; cylindrical basis height (Angstroms) |
| `bin_sz` | float | Bin size; map resolution (Angstroms) |
| `dot_skp_frq` | int | Dot-skipping frequency for reference vertices |
| `dot_map_frq` | int | Dot-skipping frequency for mapped vertices |

**Returns:**
| Name | Type | Description |
|------|------|-------------|
| `vrtx_spn_maps` | list[numpy.ndarray] | List of 2D fingerprint matrices |

**Algorithm:**
```
S_O: R^3 -> R^2
S_O(x) -> (alpha, beta)

where:
  alpha = sqrt(||x-p||^2 - (n.(x-p))^2)
  beta = n.(x-p)
```

---

#### `generate_circle_matrix_3d(input_3d_matrix, radius_fade, slope_fade)`

Apply sigmoid radial fading to maps.

**Parameters:**
| Name | Type | Description |
|------|------|-------------|
| `input_3d_matrix` | numpy.ndarray | Input 3D matrix of maps (N x H x W) |
| `radius_fade` | float | Fading radius in bins |
| `slope_fade` | float | Sigmoid slope parameter |

**Returns:**
| Name | Type | Description |
|------|------|-------------|
| `result` | numpy.ndarray | Maps with radial fading applied |

---

#### `logistic(x, k=10)`

Logistic (sigmoid) function for fading calculations.

**Parameters:**
| Name | Type | Default | Description |
|------|------|---------|-------------|
| `x` | float/array | - | Input value(s) |
| `k` | float | 10 | Slope parameter |

**Returns:**
```
1 - 1 / (1 + exp(-k * x))
```

---

#### `minmax_norm(in_arr_maps)`

Apply min-max normalization to maps.

**Parameters:**
| Name | Type | Description |
|------|------|-------------|
| `in_arr_maps` | numpy.ndarray | Input maps array (N x H x W) |

**Returns:**
| Name | Type | Description |
|------|------|-------------|
| `maps_norm_arr` | numpy.ndarray | Normalized maps with values in [0, 1] |

---

#### `hector_mapper_vctrsd(ply_fn, spprt_dstnc_x, spprt_dstnc_y, bin_sz, dot_skp_frq, dot_map_frq, sign_flg, infltn, radius_fade, slope_fade)`

Main wrapper function for surface mapping and PLY preprocessing.

**Parameters:**
| Name | Type | Description |
|------|------|-------------|
| `ply_fn` | str | Input PLY file name |
| `spprt_dstnc_x` | float | Support distance alpha (Angstroms) |
| `spprt_dstnc_y` | float | Support distance beta (Angstroms) |
| `bin_sz` | float | Bin size/resolution (Angstroms) |
| `dot_skp_frq` | int | Reference vertex-skipping frequency |
| `dot_map_frq` | int | Mapped vertex-skipping frequency |
| `sign_flg` | str | Mapping polarity: "lgnd" (inverse) or "rcpt" (forward) |
| `infltn` | float | Surface inflation (Angstroms) |
| `radius_fade` | float | Fading radius (bins) |
| `slope_fade` | float | Fading slope |

**Output:**
Saves compressed NPZ file to `/results/` with:
- `coords`: Patch center coordinates
- `nrmls`: Surface normal vectors
- `maps`: Fingerprint matrices
- `comments`: Metadata string

---

### Command Line Interface

```bash
python hector_mapper.py <ply> <dist_x> <dist_y> <bin> <dot_skp> <map_skp> <infltn> <radius_fade> <slope_fade> <sign>
```

**Arguments:**
| Argument | Type | Description |
|----------|------|-------------|
| `ply` | str | Input PLY file path |
| `dist_x` | float | Support distance alpha (default: 12.0) |
| `dist_y` | float | Support distance beta (default: 6.0) |
| `bin` | float | Bin size (default: 0.4) |
| `dot_skp` | int | Reference vertex-skip frequency (default: 40) |
| `map_skp` | int | Mapped vertex-skip frequency (default: 8) |
| `infltn` | float | Surface inflation (default: 0.0) |
| `radius_fade` | float | Radial fading radius (default: 40.0) |
| `slope_fade` | float | Fading slope (default: 0.3) |
| `sign` | str | "lgnd" or "rcpt" |

---

## maps_analysis_pairs_vs_all.py

Map comparison and search for complementary surface pairs.

**Location:** `code/maps_analysis_pairs_vs_all.py`

### Functions

---

#### `inv_maps(lgnd_maps)`

Invert ligand maps for complementarity comparison.

**Parameters:**
| Name | Type | Description |
|------|------|-------------|
| `lgnd_maps` | numpy.ndarray | Ligand/query maps array (N x H x W) |

**Returns:**
| Name | Type | Description |
|------|------|-------------|
| `lgnd_inv_maps` | numpy.ndarray | Inverted maps (flipped along axis 0) |

**Implementation:**
```python
inverted_matrix = matrix[::-1, :]  # Flip rows
```

---

#### `find_map_indcs(atm_pair, coords_arr)`

Identify map indices closest to selected surface dots.

**Parameters:**
| Name | Type | Description |
|------|------|-------------|
| `atm_pair` | tuple | Pair of 3D coordinates (atm1, atm2) |
| `coords_arr` | numpy.ndarray | Array of map center coordinates (N x 3) |

**Returns:**
| Name | Type | Description |
|------|------|-------------|
| `map_indx_lst` | list[int] | List of indices of nearest maps |

---

#### `srch_diff_v4(qry_maps_mat, sbjct_maps_mat)`

Compare query maps against subject maps using SSIM.

**Parameters:**
| Name | Type | Description |
|------|------|-------------|
| `qry_maps_mat` | numpy.ndarray | Query maps matrix (N_q x H x W) |
| `sbjct_maps_mat` | numpy.ndarray | Subject maps matrix (N_s x H x W) |

**Returns:**
| Name | Type | Description |
|------|------|-------------|
| `r_map_ij_arr` | numpy.ndarray | Array of [R-factor, query_idx, subject_idx] |

**Notes:**
- Uses compiled `sim.ssim()` function for fast comparison
- R-factor is negated SSIM (-1 * SSIM)

---

#### `map_analysis(sbjct_path, qry_path, rf_cutoff, atm_pair)`

Main analysis function implementing two-vs-all search strategy.

**Parameters:**
| Name | Type | Description |
|------|------|-------------|
| `sbjct_path` | str | Path to subject NPZ file |
| `qry_path` | str | Path to query NPZ file |
| `rf_cutoff` | float | R-factor cutoff threshold |
| `atm_pair` | tuple | Pair of query atom coordinates |

**Returns:**
| Name | Type | Description |
|------|------|-------------|
| `filt_map_pairs_rf_coords_pdb` | numpy.ndarray | Filtered results with coordinates and metadata |

**Algorithm:**
1. Load query and subject maps
2. Find query maps nearest to specified coordinates
3. Invert query maps
4. Compare against all subject maps (R-factor)
5. Select top 5% hits for each query
6. Filter by inter-patch distance
7. Filter by average R-factor
8. Return filtered pairs with coordinates

---

### Global Variables

```python
n_cores = 10        # Number of parallel jobs
top_pct = 5         # Percentage of top hits to retain
dist_tol = 0.01     # Distance tolerance (Angstroms)
```

### Command Line Interface

```bash
python maps_analysis_pairs_vs_all.py <npz_dir> <qry_npz> <rf_cutoff>
```

**Arguments:**
| Argument | Type | Description |
|----------|------|-------------|
| `npz_dir` | str | Directory with subject NPZ files |
| `qry_npz` | str | Query NPZ filename |
| `rf_cutoff` | float | R-factor cutoff (e.g., -0.82) |

---

## aln_fltr_4_dots.py

Docking, alignment, and filtering of identified scaffolds.

**Location:** `code/aln_fltr_4_dots.py`

### Classes

---

#### `class atom`

Parse individual atoms from PDB format.

**Constructor:**
```python
atom(pdbatom_str)
```

**Attributes:**
| Attribute | Type | Description |
|-----------|------|-------------|
| `idx` | int | Atom index |
| `name` | str | Atom name (e.g., "CA", "CB") |
| `resn` | str | Residue name (e.g., "ALA", "GLY") |
| `chain_id` | str | Chain identifier |
| `resid` | int | Residue number |
| `coords` | numpy.ndarray | 3D coordinates [x, y, z] |

---

#### `class pdb`

Manage PDB structure with transform operations.

**Constructor:**
```python
pdb(atms_lst)
```

**Attributes:**
| Attribute | Type | Description |
|-----------|------|-------------|
| `atms_lst` | list[atom] | List of atom objects |
| `coords` | numpy.ndarray | Coordinates array (N x 3) |

**Methods:**

##### `update()`
Update atom coordinates from coords array.

##### `transform(rot_mat, trans_vec, centre=True)`
Apply rotation and translation to structure.

**Parameters:**
| Name | Type | Default | Description |
|------|------|---------|-------------|
| `rot_mat` | numpy.ndarray | Identity | 3x3 rotation matrix |
| `trans_vec` | numpy.ndarray | [0,0,0] | Translation vector |
| `centre` | bool | True | Center structure before transform |

---

### Functions

---

#### `read_pdb(fn)`

Read PDB file and return list of atom objects.

**Parameters:**
| Name | Type | Description |
|------|------|-------------|
| `fn` | str | PDB file path |

**Returns:**
| Name | Type | Description |
|------|------|-------------|
| `atms_lst` | list[atom] | List of atom objects |

---

#### `gen_pdb_rcrd(atm, idx=None)`

Generate PDB-format ATOM record string.

**Parameters:**
| Name | Type | Description |
|------|------|-------------|
| `atm` | atom | Atom object |
| `idx` | int/None | Override atom index |

**Returns:**
| Name | Type | Description |
|------|------|-------------|
| `atm_str` | str | PDB ATOM record (80 characters) |

---

#### `write_pdb(o_fn, atms_lst, rmrk=None)`

Write atom list to PDB file.

**Parameters:**
| Name | Type | Description |
|------|------|-------------|
| `o_fn` | str | Output filename |
| `atms_lst` | list[atom] | List of atom objects |
| `rmrk` | str/None | REMARK text |

**Returns:** 0 on success

---

#### `kabsch_v2(P, Q)`

Kabsch algorithm for optimal 3D alignment.

**Parameters:**
| Name | Type | Description |
|------|------|-------------|
| `P` | numpy.ndarray | First point set (N x 3) |
| `Q` | numpy.ndarray | Second point set (N x 3) |

**Returns:**
| Name | Type | Description |
|------|------|-------------|
| `P_alnd_trnsltd` | numpy.ndarray | Aligned and translated P |
| `rotation_matrix` | numpy.ndarray | 3x3 rotation matrix |
| `tvec_ref` | numpy.ndarray | Translation to reference |
| `tvec_mbl` | numpy.ndarray | Translation from mobile |
| `rmsd` | float | Root mean square deviation |

**Algorithm:**
1. Center both point sets
2. Compute covariance matrix H = P^T * Q
3. SVD: H = U * S * V^T
4. Rotation: R = U * V^T
5. Handle reflection if det(R) < 0
6. Calculate RMSD

---

#### `anchor_4_dots(s_mbl, q_ref, in_mbl_fn, dmp_alnd_pdb=False)`

Align candidate structures using 4 anchor points.

**Parameters:**
| Name | Type | Description |
|------|------|-------------|
| `s_mbl` | numpy.ndarray | Mobile (subject) anchor points (4 x 3) |
| `q_ref` | numpy.ndarray | Reference (query) anchor points (4 x 3) |
| `in_mbl_fn` | str | Input mobile PDB file path |
| `dmp_alnd_pdb` | bool | Write aligned PDB file |

**Returns:**
| Name | Type | Description |
|------|------|-------------|
| `tuple` | (pdb, float) | (Aligned structure, RMSD) or None |

---

#### `calc_ovrlp(qry_atms_lst, sbj_atms_lst, slv_rsltn, slv_n_vxls, slv_dict, spprt_dstnc)`

Quantify surface overlap and interface residues.

**Parameters:**
| Name | Type | Description |
|------|------|-------------|
| `qry_atms_lst` | list[atom] | Query structure atoms |
| `sbj_atms_lst` | list[atom] | Subject structure atoms |
| `slv_rsltn` | float | Solvation kernel resolution |
| `slv_n_vxls` | int | Number of voxels |
| `slv_dict` | dict | Solvation kernel dictionary |
| `spprt_dstnc` | float | Support distance |

**Returns:**
| Name | Type | Description |
|------|------|-------------|
| `ovrlp` | int | Number of overlapping voxels |
| `n_interface_res` | int | Number of interface residues |

---

#### `dock_hits(hit, pdbs_dir, in_qry_pdb_fn, ovrlp_co, intrfc_resids_co, rmsd_co, slv_rsltn, slv_n_vxls, slv_dict, spprt_dstnc)`

Main docking function with filtering.

**Parameters:**
| Name | Type | Description |
|------|------|-------------|
| `hit` | numpy.ndarray | Hit data from search results |
| `pdbs_dir` | str | Directory with PDB files |
| `in_qry_pdb_fn` | str | Query PDB file path |
| `ovrlp_co` | int | Maximum overlap cutoff |
| `intrfc_resids_co` | int | Minimum interface residues |
| `rmsd_co` | float | Maximum RMSD cutoff |
| `slv_*` | various | Solvation kernel parameters |

**Returns:**
| Name | Type | Description |
|------|------|-------------|
| `result` | list/None | [hit, rmsd, overlap, interface_res, query_pdb] |

**Filtering criteria:**
- `ovrlp <= ovrlp_co` (default: 125,000)
- `intrfc_resids >= intrfc_resids_co` (default: 25)
- `rmsd <= rmsd_co` (default: 0.5)

---

### Command Line Interface

```bash
python aln_fltr_4_dots.py <hits_fn> <pdbs_dir>
```

**Arguments:**
| Argument | Type | Description |
|----------|------|-------------|
| `hits_fn` | str | Path to search results (srch_rslts.npy) |
| `pdbs_dir` | str | Directory with scaffold PDB files |

---

## sim.pyx

Cython implementation of SSIM (Structural Similarity Index) calculation.

**Location:** `code/sim.pyx`

### Functions

---

#### `ssim(im1, im2, data_range=1, full=False, K1=0.01)`

Calculate element-wise structural similarity between 2D maps.

**Parameters:**
| Name | Type | Default | Description |
|------|------|---------|-------------|
| `im1` | numpy.ndarray[float32] | - | First image (H x W) |
| `im2` | numpy.ndarray[float32] | - | Second image (H x W) |
| `data_range` | float | 1 | Data range for SSIM |
| `full` | bool | False | Return full SSIM map |
| `K1` | float | 0.01 | SSIM constant |

**Returns:**
| Name | Type | Description |
|------|------|-------------|
| `mssim` | float | Mean SSIM value |
| `S` | numpy.ndarray | Full SSIM map (if full=True) |

**Formula:**
```
S[i,j] = (2 * im1[i,j] * im2[i,j] + C1) / (im1[i,j]^2 + im2[i,j]^2 + C1)

where C1 = (K1 * data_range)^2
```

**Performance:**
- Bounds checking disabled (`@cython.boundscheck(False)`)
- Wraparound disabled (`@cython.wraparound(False)`)

---

## Data Structures

### NPZ File Format (Fingerprints)

```python
# Loading
data = np.load('protein_rcpt.npz', allow_pickle=True)

# Contents
data['coords']    # numpy.ndarray, shape (N, 3) - Map center coordinates
data['nrmls']     # numpy.ndarray, shape (N, 3) - Surface normal vectors
data['maps']      # numpy.ndarray, shape (N, H, W) - Fingerprint matrices
data['comments']  # str - Metadata
```

### Search Results Array (srch_rslts.npy)

| Column | Content | Type |
|--------|---------|------|
| 0 | Distance between subject maps | float |
| 1 | Subject map 1 index | int |
| 2 | Subject map 2 index | int |
| 3 | Average R-factor | float |
| 4-6 | Subject map 1 coordinates | float |
| 7-9 | Subject map 1 normals | float |
| 10-12 | Subject map 2 coordinates | float |
| 13-15 | Subject map 2 normals | float |
| 16-18 | Query map 1 coordinates | float |
| 19-21 | Query map 1 normals | float |
| 22-24 | Query map 2 coordinates | float |
| 25-27 | Query map 2 normals | float |
| 28 | NPZ filename | str |

### Docking Results Array (dock_rslts.npy)

```python
# Each element is a list:
[
    hit,           # Original hit array from search
    rmsd,          # float - Alignment RMSD
    overlap,       # int - Voxel overlap count
    interface_res, # int - Number of interface residues
    query_pdb      # str - Query PDB filename
]
```

### Solvation Kernels (solv_krnls_*.npz)

```python
# Loading
slv_rsltn, slv_n_vxls, slv_dict = load_solv_krnls(solv_krnls_fn)

# Contents
slv_rsltn   # float - Voxel resolution (e.g., 0.5 Angstroms)
slv_n_vxls  # int - Number of voxels per dimension (e.g., 24)
slv_dict    # dict - Atom type to kernel mapping
```

---

## Error Codes

| Code | Description |
|------|-------------|
| 0 | Success |
| 2 | Invalid arguments or runtime error |

---

## Dependencies

### Python Packages

```python
import numpy as np          # Array operations
import scipy                # Scientific computing
from joblib import Parallel # Parallel processing
```

### Cython Modules

```python
import sim                  # SSIM calculation
from tensors import *       # Solvation kernel operations
```
