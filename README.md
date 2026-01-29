# HECTOR

**High-throughput Epitope Complementarity-based Target-Oriented Ranking**

An ultra-fast computational algorithm for de novo protein binder design through surface complementarity evaluation.

[![DOI](https://img.shields.io/badge/DOI-10.1002%2Fadvs.202502015-blue)](https://doi.org/10.1002/advs.202502015)
[![License](https://img.shields.io/badge/License-Confidential-red)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.9-green)](https://python.org)

---

## Overview

HECTOR enables the identification of complementary protein scaffolds from structure databases by evaluating surface complementarity through vectorized, lower-dimensional surface representations called "fingerprints". The algorithm was developed as part of a complementarity-first approach to protein binder design.

### Key Features

- **Ultra-fast docking**: Sub-microsecond map-to-map comparisons on GPU
- **Training-free**: No reliance on machine learning or training data
- **Invertible fingerprints**: Query surfaces trivially inverted to find complementary matches
- **High success rate**: >25% nanomolar binders identified in experimental validation

### Publication

> **"A Complementarity-Based Approach to De Novo Binder Design"**
> Maksymenko K, Hatskovska V, Coles M, et al.
> *Advanced Science*, 2025, 12, e02015
> DOI: [10.1002/advs.202502015](https://doi.org/10.1002/advs.202502015)

---

## Quick Start

### Prerequisites

- Python 3.9+
- NumPy, SciPy, Cython, Joblib
- [EDTSurf](https://zhanggroup.org/EDTSurf/) for surface mesh generation

### Installation

```bash
# Using Docker (recommended)
cd environment
docker build -t hector:latest .
docker run -it -v $(pwd):/workspace hector:latest

# Or using conda
conda create -n hector python=3.9
conda activate hector
conda install cython=3.0.11 numpy=1.23.4 scipy=1.12.0 joblib=1.4.2
```

### Basic Usage

```bash
# Step 1: Generate surface mesh (requires EDTSurf)
EDTSurf -i protein.pdb -o protein.ply -s 3

# Step 2: Generate surface fingerprints
# Standard version (single-threaded)
python code/hector_mapper.py protein.ply 10 20 0.2 5 1 0.5 40 0.3 rcpt

# OR: Parallel version for faster processing (recommended)
python code/hector_mapper_parallel.py protein.ply 10 20 0.2 40 8 0.5 40 0.3 rcpt

# Step 3: Search scaffold database
python code/maps_analysis_pairs_vs_all.py /data/scaffolds_db/ protein_rcpt.npz -0.82

# Step 4: Dock and filter candidates
python code/aln_fltr_4_dots.py /results/srch_rslts.npy /data/scaffolds_db
```

---

## Project Structure

```
HECTOR/
├── code/                           # Source code
│   ├── hector_mapper.py            # Surface fingerprinting engine
│   ├── hector_mapper_parallel.py   # Parallelized fingerprinting (NEW)
│   ├── maps_analysis_pairs_vs_all.py   # Database search module
│   ├── aln_fltr_4_dots.py          # Docking and filtering
│   ├── sim.pyx                     # Cython SSIM implementation
│   ├── invoke_hector.sh            # Example execution script
│   └── solv_krnls_0.50A_24vxl.npz  # Solvation kernels
├── data/                           # Input data (not tracked)
│   └── scaffolds_db/               # Pre-computed scaffold database
├── results/                        # Output results (not tracked)
├── docs/                           # Documentation
│   ├── README.md                   # Documentation index
│   ├── INSTALLATION.md             # Setup instructions
│   ├── ALGORITHM.md                # Algorithm details
│   ├── PIPELINE.md                 # Usage guide
│   ├── API_REFERENCE.md            # Code documentation
│   └── SCIENTIFIC_BACKGROUND.md    # Research context
├── environment/                    # Docker configuration
│   └── Dockerfile
├── reference/                      # Research article
└── metadata/                       # Project metadata
```

---

## Pipeline Overview

```
┌─────────────┐     ┌─────────────┐     ┌─────────────┐     ┌─────────────┐
│  PDB File   │────>│ PLY Surface │────>│ Fingerprint │────>│   Search    │
│             │     │   (EDTSurf) │     │    (NPZ)    │     │   Results   │
└─────────────┘     └─────────────┘     └─────────────┘     └─────────────┘
                                                                   │
                    ┌─────────────┐     ┌─────────────┐            │
                    │  Aligned    │<────│   Docking   │<───────────┘
                    │    PDBs     │     │  & Filter   │
                    └─────────────┘     └─────────────┘
```

### Stage 1: Surface Mapping
- Input: Protein structure (PDB) → Surface mesh (PLY)
- Process: Generate spin image fingerprints from surface geometry
- Output: Compressed NPZ files with coordinates, normals, and 2D maps

### Stage 2: Database Search
- Input: Query maps + scaffold database maps
- Process: Compare query maps against all subjects using SSIM
- Filter: Inter-patch distance tolerance and R-factor threshold
- Output: Candidate hit pairs sorted by complementarity

### Stage 3: Docking & Validation
- Input: Hit pairs + original structures
- Process: Kabsch alignment, RMSD calculation, overlap quantification
- Filter: RMSD < 0.5 Å, overlap < 125,000 voxels, interface residues > 25
- Output: Aligned PDB structures of validated candidates

---

## Algorithm Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `dist_x` | 10.0 | Radial support distance (Å) |
| `dist_y` | 20.0 | Axial support distance (Å) |
| `bin` | 0.2 | Map resolution (Å) |
| `dot_skp` | 5 | Reference vertex-skip frequency |
| `map_skp` | 1 | Mapped vertex-skip frequency |
| `infltn` | 0.5 | Surface inflation (Å) |
| `rf_cutoff` | -0.82 | R-factor threshold |

---

## Performance Optimization

### Parallel Processing (NEW)

HECTOR now includes a parallelized version of the fingerprinting engine for multi-core systems:

```bash
# Use all CPU cores for faster processing
python code/hector_mapper_parallel.py protein.ply 10 20 0.2 40 8 0.5 40 0.3 rcpt

# Specify number of cores
python code/hector_mapper_parallel.py protein.ply 10 20 0.2 40 8 0.5 40 0.3 rcpt --n_jobs 8
```

**Expected speedup:**
- 8-core system: ~8× faster
- 16-core system: ~16× faster

### Optimization Strategies

1. **Use recommended skipping frequencies**: `dot_skp=40, map_skp=8` (64× fewer vertices than default)
2. **Enable parallel processing**: `hector_mapper_parallel.py` with `--n_jobs -1`
3. **Combine both**: 512× speedup for high-throughput screening

**Trade-offs:**
- Higher skipping frequencies reduce resolution but maintain biological relevance
- Parallel processing increases memory usage (~N× for N cores)
- For large surfaces (>500K vertices), adjust parameters accordingly

---

## Experimental Results

### VEGF Binders

| Design | Scaffold | Kd (nM) | Tm (°C) |
|--------|----------|---------|---------|
| Sam0.7 | Ketosteroid isomerase | 190 | 63 |
| Sima3.2 | Nitrophorin | 14 | 65 |

### IL-7Rα Binders

| Design | Scaffold | Kd (nM) |
|--------|----------|---------|
| des01 | TIM barrel | 26 |
| des03 | Four-helix bundle | 1.4 |
| des07 | Rossmann fold | 20 |

---

## Documentation

Detailed documentation is available in the [docs/](docs/) folder:

- [Installation Guide](docs/INSTALLATION.md) - Setup and dependencies
- [Algorithm Description](docs/ALGORITHM.md) - Technical details
- [Pipeline Usage](docs/PIPELINE.md) - Step-by-step workflow
- [API Reference](docs/API_REFERENCE.md) - Code documentation
- [Scientific Background](docs/SCIENTIFIC_BACKGROUND.md) - Research context

---

## Citation

If you use HECTOR in your research, please cite:

```bibtex
@article{maksymenko2025complementarity,
  title={A Complementarity-Based Approach to De Novo Binder Design},
  author={Maksymenko, Kateryna and Hatskovska, Valeriia and Coles, Murray and others},
  journal={Advanced Science},
  volume={12},
  pages={e02015},
  year={2025},
  publisher={Wiley-VCH GmbH},
  doi={10.1002/advs.202502015}
}
```

---

## Data Availability

- **HECTOR Software**: [CodeOcean](https://doi.org/10.24433/CO.9243108.v1)
- **Crystal Structures**: PDB [8BL5](https://www.rcsb.org/structure/8BL5), [8BL9](https://www.rcsb.org/structure/8BL9)
- **MD Trajectories**: [Zenodo](https://doi.org/10.5281/zenodo.14028991)

---

## License

- **Source Code**: Confidential (for peer-review only)
- **Data**: [CC-BY-NC-ND 4.0 International](https://creativecommons.org/licenses/by-nc-nd/4.0/)

---

## Authors

- **Mohammad ElGamacy** - University Hospital Tübingen
- **Kateryna Maksymenko** - Max Planck Institute for Biology
- **Valeriia Hatskovska** - University Hospital Tübingen

For questions, contact: Mohammad.Elgamacy@med.uni-tuebingen.de

---

## Acknowledgments

This work was supported by:
- Horizon Europe European Research Council (Grant No.: 863952)
- Deutsche Forschungsgemeinschaft (Grant No.: 500215849)
- Max-Planck-Gesellschaft
