# HECTOR Documentation

## High-throughput Epitope Complementarity-based Target-Oriented Ranking

**Version:** 0.2
**Author:** Mohammad ElGamacy (University Hospital Tubingen)
**Copyright:** 2018 Mohammad ElGamacy / Max Planck Society
**License:** Confidential (source code); CC-BY-NC-ND 4.0 (data)

---

## Overview

HECTOR (High-throughput Epitope Complementarity-based Target-Oriented Ranking) is an ultra-fast computational algorithm for de novo protein binder design. It enables the identification of complementary scaffolds from protein structure databases by evaluating surface complementarity through a vectorized, lower-dimensional surface representation called "fingerprints" or "maps."

The algorithm was developed as part of a complementarity-first approach to protein binder design, as described in the research article:

> **"A Complementarity-Based Approach to De Novo Binder Design"**
> Maksymenko K, Hatskovska V, Coles M, et al.
> *Advanced Science*, 2025, 12, e02015
> DOI: 10.1002/advs.202502015

---

## Key Features

- **Ultra-fast docking evaluation**: Sub-microsecond map-to-map comparisons on GPU
- **Training-free approach**: No reliance on machine learning or training data
- **Invertible fingerprints**: Query surfaces can be trivially inverted to find complementary matches
- **Dimensionality reduction**: 3D surface patches compressed to 2D maps for efficient comparison
- **Rotational invariance**: Fingerprints are rotation-invariant through cylindrical basis projection
- **High success rate**: Demonstrated nanomolar binders against VEGF and IL-7Ralpha targets

---

## Documentation Contents

| Document | Description |
|----------|-------------|
| [Installation Guide](INSTALLATION.md) | Setup instructions and dependencies |
| [Algorithm Description](ALGORITHM.md) | Detailed explanation of the HECTOR algorithm |
| [Pipeline Usage](PIPELINE.md) | Step-by-step workflow guide |
| [API Reference](API_REFERENCE.md) | Source code documentation |
| [Scientific Background](SCIENTIFIC_BACKGROUND.md) | Research context and methodology |

---

## Quick Start

### 1. Prerequisites

- Python 3.9+
- Docker (recommended) or conda environment
- EDTSurf for surface mesh generation

### 2. Basic Workflow

```bash
# Step 1: Generate surface fingerprints for query protein
python hector_mapper.py /data/query.ply 10 20 0.2 5 1 0.5 40 0.3 rcpt

# Step 2: Search scaffold database for complementary surfaces
python maps_analysis_pairs_vs_all.py /data/scaffolds_db/ query_rcpt.npz -0.82

# Step 3: Dock and filter candidate scaffolds
python aln_fltr_4_dots.py /results/srch_rslts.npy /data/scaffolds_db
```

### 3. Output

- `*_rcpt.npz`: Surface fingerprints with coordinates and normals
- `srch_rslts.npy`: Search results with hit pairs
- `dock_rslts.npy`: Filtered docking results
- `*_alnd.pdb`: Aligned scaffold structures

---

## Project Structure

```
HECTOR/
├── code/                          # Core implementation
│   ├── hector_mapper.py           # Surface fingerprinting engine
│   ├── maps_analysis_pairs_vs_all.py  # Database search module
│   ├── aln_fltr_4_dots.py         # Docking and filtering module
│   ├── sim.pyx                    # Cython SSIM implementation
│   ├── tensors.cpython-*.so       # Compiled tensor operations
│   └── invoke_hector.sh           # Example execution script
├── data/                          # Input data
│   ├── il7ra.pdb                  # Example query structure
│   └── scaffolds_db/              # Pre-computed scaffold database
├── results/                       # Output results
├── environment/                   # Docker configuration
├── reference/                     # Research article
└── docs/                          # Documentation (this folder)
```

---

## Biological Applications

HECTOR has been successfully applied to design:

1. **Anti-VEGF binders** (Sam and Sima scaffolds)
   - Nanomolar affinity (Kd: 14-190 nM)
   - Demonstrated anti-angiogenic activity in vitro
   - Tumor-inhibiting activity in vivo (zebrafish xenografts)

2. **Anti-IL-7Ralpha binders** (des01-des08)
   - Multiple scaffold folds (TIM barrel, four-helix bundle, Rossmann fold)
   - Sub-nanomolar to nanomolar affinities (Kd: 1.4-26 nM)
   - Inhibition of IL-7 signaling pathway

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

- **HECTOR software**: [CodeOcean Repository](https://doi.org/10.24433/CO.9243108.v1)
- **Crystal structures**: PDB accession codes 8BL5 (Sam0.2) and 8BL9 (Sam0.7)
- **MD trajectories**: [Zenodo Repository](https://doi.org/10.5281/zenodo.14028991)

---

## Support

For questions or issues, please refer to the original publication or contact the authors at the University Hospital Tubingen.
