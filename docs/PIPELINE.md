# HECTOR Pipeline Usage Guide

This guide provides step-by-step instructions for running the HECTOR pipeline to identify complementary scaffolds for protein binder design.

---

## Table of Contents

1. [Pipeline Overview](#pipeline-overview)
2. [Prerequisites](#prerequisites)
3. [Step 1: Surface Generation](#step-1-surface-generation)
4. [Step 2: Surface Mapping](#step-2-surface-mapping)
5. [Step 3: Database Search](#step-3-database-search)
6. [Step 4: Docking and Filtering](#step-4-docking-and-filtering)
7. [Output Analysis](#output-analysis)
8. [Complete Example](#complete-example)

---

## Pipeline Overview

```
┌──────────────────┐     ┌──────────────────┐     ┌──────────────────┐
│  PDB Structure   │────>│   PLY Surface    │────>│  NPZ Fingerprint │
│  (query.pdb)     │     │  (query.ply)     │     │  (query_rcpt.npz)│
└──────────────────┘     └──────────────────┘     └──────────────────┘
        │                        EDTSurf              hector_mapper.py
        │
        v
┌──────────────────────────────────────────────────────────────────────┐
│                      Scaffold Database (pre-computed)                │
│  ┌────────────┐  ┌────────────┐  ┌────────────┐  ┌────────────┐     │
│  │ 5djl_rcpt  │  │ 5nlc_rcpt  │  │ 7z64_rcpt  │  │ 8brb_rcpt  │ ... │
│  └────────────┘  └────────────┘  └────────────┘  └────────────┘     │
└──────────────────────────────────────────────────────────────────────┘
        │
        │ maps_analysis_pairs_vs_all.py
        v
┌──────────────────┐     ┌──────────────────┐     ┌──────────────────┐
│  Search Results  │────>│  Docked Hits     │────>│  Aligned PDBs    │
│  (srch_rslts.npy)│     │ (dock_rslts.npy) │     │  (*_alnd.pdb)    │
└──────────────────┘     └──────────────────┘     └──────────────────┘
                             aln_fltr_4_dots.py
```

---

## Prerequisites

### Required Files

1. **Query protein structure** (PDB format)
2. **EDTSurf executable** for surface generation
3. **Scaffold database** with pre-computed fingerprints
4. **Solvation kernels** (`solv_krnls_0.50A_24vxl.npz`)

### Directory Structure

```bash
HECTOR/
├── code/
│   ├── hector_mapper.py
│   ├── maps_analysis_pairs_vs_all.py
│   ├── aln_fltr_4_dots.py
│   ├── sim.cpython-39-x86_64-linux-gnu.so
│   ├── tensors.cpython-39-x86_64-linux-gnu.so
│   └── solv_krnls_0.50A_24vxl.npz
├── data/
│   ├── query.pdb          # Your target protein
│   └── scaffolds_db/      # Pre-computed database
│       ├── scaffold1.pdb
│       ├── scaffold1_rcpt.npz
│       └── ...
└── results/               # Output directory
```

---

## Step 1: Surface Generation

Generate a molecular surface mesh from the PDB structure using EDTSurf.

### Command

```bash
EDTSurf -i /data/query.pdb -o /data/query.ply -s 3
```

### Parameters

| Option | Value | Description |
|--------|-------|-------------|
| `-i` | input.pdb | Input PDB file path |
| `-o` | output.ply | Output PLY file path |
| `-s` | 3 | Surface type (3 = molecular surface) |

### Output

**PLY file format:**
```
ply
format ascii 1.0
element vertex 12345
property float x
property float y
property float z
property float nx
property float ny
property float nz
element face 24680
property list uchar int vertex_indices
end_header
[vertex data...]
[face data...]
```

### Troubleshooting

- **Missing atoms**: Add missing atoms with PyMOL or pdbfixer
- **Multiple chains**: Process each chain separately
- **Large proteins**: May require increased memory

---

## Step 2: Surface Mapping

Generate surface fingerprints (maps) from the PLY surface mesh.

### Command

```bash
cd /path/to/HECTOR/code

python hector_mapper.py \
    /data/query.ply \    # Input PLY file
    10                \   # Support distance X (radius, Angstroms)
    20                \   # Support distance Y (height, Angstroms)
    0.2               \   # Bin size (resolution, Angstroms)
    5                 \   # Dot skip frequency (mapping every Nth dot)
    1                 \   # Map skip frequency (use every dot)
    0.5               \   # Surface inflation (Angstroms)
    40                \   # Radial fade radius (bins)
    0.3               \   # Fade slope
    rcpt                  # Mapping direction (rcpt = receptor/forward)
```

### Parameters Explained

| Parameter | Typical Value | Description |
|-----------|---------------|-------------|
| `ply` | - | Path to input PLY file |
| `dist_x` | 10.0 | Radial support distance (cylindrical basis radius) |
| `dist_y` | 20.0 | Axial support distance (cylindrical basis height) |
| `bin` | 0.2 | Map resolution (smaller = higher resolution) |
| `dot_skp` | 5 | Generate map every N dots (5 = ~4 maps/A^2) |
| `map_skp` | 1 | Use every N-th dot in projection (1 = all dots) |
| `infltn` | 0.5 | Surface inflation to account for uncertainty |
| `radius_fade` | 40 | Sigmoid fading radius in bins |
| `slope_fade` | 0.3 | Sigmoid fading slope |
| `sign` | rcpt/lgnd | Mapping direction (rcpt=forward, lgnd=inverse) |

### Output

**File:** `query_rcpt.npz` (NumPy compressed archive)

**Contents:**
- `coords`: Map center coordinates (N x 3 array)
- `nrmls`: Surface normal vectors (N x 3 array)
- `maps`: Fingerprint matrices (N x 50 x 100 array)
- `comments`: Metadata string

### Inspect Output

```python
import numpy as np

data = np.load('/results/query_rcpt.npz', allow_pickle=True)
print(f"Number of maps: {data['maps'].shape[0]}")
print(f"Map dimensions: {data['maps'].shape[1:]}")
print(f"Coordinates shape: {data['coords'].shape}")
print(f"Comments: {data['comments']}")
```

---

## Step 3: Database Search

Search the scaffold database for surface patches complementary to the query epitope.

### Select Query Points

First, identify 2-3 points at your target epitope. These should be:
- Located at the intended binding interface
- Separated by at least 10 Angstroms
- Away from glycosylation sites

**Example for IL-7Ralpha (from `maps_analysis_pairs_vs_all.py`):**
```python
# Coordinates of epitope points
atm0 = np.array([12.619, 18.25, -15.497])
atm1 = np.array([9.304, 19.759, -4.274])
atm2 = np.array([-2.767, 20.462, -7.706])

# All pairwise combinations
atm_prs = [(atm0, atm1), (atm1, atm2), (atm2, atm0)]
```

### Command

```bash
python maps_analysis_pairs_vs_all.py \
    /data/scaffolds_db/    \  # Directory with scaffold NPZ files
    query_rcpt.npz         \  # Query fingerprints
    -0.82                     # R-factor cutoff
```

### Parameters

| Parameter | Description |
|-----------|-------------|
| `npz_dir` | Directory containing subject scaffold NPZ files |
| `qry_npz` | Query protein fingerprint file (must be in npz_dir) |
| `rf_cutoff` | R-factor cutoff (recommended: -0.82) |

### Internal Parameters

```python
# In maps_analysis_pairs_vs_all.py
n_cores = 10           # Number of parallel processes
top_pct = 5            # Select top 5% of hits per query
dist_tol = 0.01        # Distance tolerance (Angstroms)
```

### Output

**File:** `srch_rslts.npy`

**Array columns:**
| Column | Content |
|--------|---------|
| 0 | Distance between two subject maps |
| 1 | Index of subject map 1 |
| 2 | Index of subject map 2 |
| 3 | Average R-factor |
| 4-6 | Coordinates of subject map 1 |
| 7-9 | Normal of subject map 1 |
| 10-12 | Coordinates of subject map 2 |
| 13-15 | Normal of subject map 2 |
| 16-18 | Coordinates of query map 1 |
| 19-21 | Normal of query map 1 |
| 22-24 | Coordinates of query map 2 |
| 25-27 | Normal of query map 2 |
| 28 | NPZ filename |

### Inspect Results

```python
import numpy as np

results = np.load('/results/srch_rslts.npy', allow_pickle=True)
print(f"Number of hits: {len(results)}")
print(f"R-factor range: {results[:, 3].min():.3f} to {results[:, 3].max():.3f}")

# Top 10 hits by R-factor
sorted_idx = np.argsort(results[:, 3])
for i in sorted_idx[:10]:
    print(f"  {results[i, 28]}: R={results[i, 3]:.3f}")
```

---

## Step 4: Docking and Filtering

Dock identified scaffolds against the target and filter by quality metrics.

### Command

```bash
python aln_fltr_4_dots.py \
    /results/srch_rslts.npy    \  # Search results from Step 3
    /data/scaffolds_db            # Directory with scaffold PDB files
```

### Filtering Parameters

```python
# In aln_fltr_4_dots.py
ovrlp_co = 125000         # Maximum surface overlap (voxels)
intrfc_resids_co = 25     # Minimum interface residues
rmsd_co = 0.5             # Maximum RMSD (Angstroms)
```

### Output Files

**1. Docking results:** `dock_rslts.npy`

Contains: [hit_data, rmsd, overlap, interface_residues, query_pdb]

**2. Aligned structures:** `*_alnd.pdb`

PDB files of scaffolds aligned to the target position.

### Inspect Docking Results

```python
import numpy as np

dock = np.load('/results/dock_rslts.npy', allow_pickle=True)
print(f"Validated hits: {len(dock)}")

for hit in dock:
    print(f"  Scaffold: {hit[0][-1][:6]}")
    print(f"    R-factor: {hit[0][3]:.3f}")
    print(f"    RMSD: {hit[1]:.3f} A")
    print(f"    Overlap: {hit[2]} voxels")
    print(f"    Interface residues: {hit[3]}")
```

---

## Output Analysis

### Visualizing Results

Use PyMOL, Chimera, or other molecular visualization software:

```bash
# PyMOL example
pymol /data/query.pdb /results/5nlc_A_080_alnd.pdb
```

### Ranking Candidates

Prioritize scaffolds by:

1. **R-factor**: Lower is better (more complementary)
2. **Interface residues**: Higher is better (larger interface)
3. **RMSD**: Lower is better (better structural fit)
4. **Overlap**: Lower is better (fewer clashes)

### Suggested Scoring

```python
def score_hit(hit):
    """Combined score for ranking hits (lower is better)."""
    r_factor = hit[0][3]
    rmsd = hit[1]
    overlap = hit[2] / 100000  # Normalize
    interface = -hit[3] / 50   # Higher is better (negative)

    return 3 * r_factor + rmsd + overlap + interface
```

---

## Complete Example

### IL-7Ralpha Binder Design

```bash
#!/bin/bash
# invoke_hector.sh - Complete HECTOR pipeline

# Step 1: Generate surface (if not already done)
# EDTSurf -i /data/il7ra.pdb -o /data/il7ra.ply -s 3

# Step 2: Generate fingerprints
python hector_mapper.py \
    /data/il7ra.ply \
    10 20 0.2 5 1 0.5 40 0.3 rcpt

# Step 3: Search database
python maps_analysis_pairs_vs_all.py \
    /data/scaffolds_db/ \
    il7ra_rcpt.npz \
    -0.82

# Step 4: Dock and filter
python aln_fltr_4_dots.py \
    /results/srch_rslts.npy \
    /data/scaffolds_db
```

### Expected Output

```
HECTOR Mapper (CPU) v13  - Copyright (C) 2018 Mohammad ElGamacy
generating local maps ...
1000
2000
...
mapped 45678 patches of il7ra in 423.5 sec

Processing scaffold database...
5djl_A_rcpt.npz
5nlc_A_rcpt.npz
...

hit: 5nlc_A  rf: -0.856  rmsd: 0.342  ovrlp: 98432  intrfc: 31
hit: 7z64_B  rf: -0.823  rmsd: 0.498  ovrlp: 115678  intrfc: 28
...
```

---

## Troubleshooting

### No Hits Found

1. **Relax R-factor cutoff**: Try -0.75 or -0.70
2. **Expand query points**: Use different epitope locations
3. **Check surface quality**: Inspect PLY file for artifacts

### Too Many Hits

1. **Tighten R-factor cutoff**: Use -0.85 or stricter
2. **Reduce distance tolerance**: Use 0.005 instead of 0.01
3. **Increase interface requirements**: Require more residues

### Memory Issues

1. **Reduce parallel cores**: Set `n_cores = 4` or lower
2. **Process databases in batches**: Split large databases
3. **Use Docker with memory limits**: `docker run --memory=16g`

### Slow Performance

1. **Use pre-compiled modules**: Ensure .so files match Python version
2. **Enable GPU acceleration**: Install ArrayFire + CUDA
3. **Reduce mapping density**: Use `dot_skp=10` for initial screening

---

## Next Steps

After identifying candidate scaffolds:

1. **Interface design**: Use Rosetta or Damietta for sequence optimization
2. **MD validation**: Run molecular dynamics to assess stability
3. **Experimental testing**: Express and characterize top candidates

See the original publication for detailed design protocols.
