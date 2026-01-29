# Installation Guide

This guide covers the installation and setup of the HECTOR software for protein binder design.

---

## System Requirements

### Hardware
- **CPU**: Multi-core processor (10+ cores recommended for parallel processing)
- **RAM**: 16 GB minimum (64 GB recommended for large databases)
- **GPU**: NVIDIA GPU with CUDA support (optional, for accelerated R-factor calculations)
- **Storage**: 10+ GB for scaffold databases

### Operating System
- Linux (Ubuntu 20.04 recommended)
- macOS (with Homebrew)
- Windows (via WSL2)

---

## Installation Methods

### Method 1: Docker (Recommended)

Docker provides the easiest and most reproducible setup.

#### 1. Build the Docker Image

```bash
cd /path/to/HECTOR
docker build -f docker/Dockerfile -t hector:latest .
```

#### 2. Run HECTOR in Docker

```bash
# Mount the HECTOR directory to /workspace in the container
docker run -it -v $(pwd):/workspace hector:latest

# Inside the container, you'll be in /workspace (the project root)
# All paths are relative to this directory
```

#### Docker Image Contents

The Dockerfile installs:
- **Base**: Ubuntu 20.04 with Miniconda3, Python 3.9
- **Compiler**: GCC 9.3.0
- **Python packages**:
  - Cython 3.0.11
  - NumPy 1.23.4
  - SciPy 1.12.0
  - Matplotlib 3.6.2
  - Jupyter 1.0.0
  - Joblib 1.4.2

---

### Method 2: Conda Environment

#### 1. Create Environment

```bash
conda create -n hector python=3.9
conda activate hector
```

#### 2. Install Dependencies

```bash
conda install -c conda-forge \
    cython=3.0.11 \
    numpy=1.23.4 \
    scipy=1.12.0 \
    matplotlib=3.6.2 \
    joblib=1.4.2 \
    jupyter=1.0.0
```

#### 3. Install GCC Compiler

```bash
# Ubuntu/Debian
sudo apt-get install gcc

# macOS
brew install gcc
```

#### 4. Compile Cython Extensions

```bash
cd /path/to/HECTOR/code

# Compile sim.pyx
cython sim.pyx
gcc -shared -pthread -fPIC -fwrapv -O2 -Wall \
    -fno-strict-aliasing -I/path/to/python/include \
    -L/path/to/python/lib \
    -o sim.cpython-39-x86_64-linux-gnu.so sim.c \
    -lpython3.9
```

---

### Method 3: Manual Installation

#### 1. Install Python 3.9

```bash
# Ubuntu
sudo apt-get install python3.9 python3.9-dev python3.9-venv

# Create virtual environment
python3.9 -m venv hector_env
source hector_env/bin/activate
```

#### 2. Install pip packages

```bash
pip install numpy==1.23.4 scipy==1.12.0 cython==3.0.11 \
    matplotlib==3.6.2 joblib==1.4.2 jupyter==1.0.0
```

---

## External Dependencies

### EDTSurf (Surface Mesh Generation)

EDTSurf is required to generate PLY surface meshes from PDB files.

#### Installation

```bash
# Download EDTSurf
wget https://zhanggroup.org/EDTSurf/EDTSurf.zip
unzip EDTSurf.zip
cd EDTSurf

# Compile (Linux)
g++ -O3 -o EDTSurf EDTSurf.cpp

# Add to PATH
export PATH=$PATH:/path/to/EDTSurf
```

#### Usage

```bash
# Generate molecular surface
EDTSurf -i input_structure.pdb -o surface_file.ply -s 3
```

**Options:**
- `-i`: Input PDB file
- `-o`: Output PLY file
- `-s 3`: Surface type (3 = molecular surface)

---

## Pre-compiled Cython Modules

The repository includes pre-compiled Cython modules for Python 3.9 on Linux x86_64:

| Module | File | Description |
|--------|------|-------------|
| SSIM | `sim.cpython-39-x86_64-linux-gnu.so` | Structural similarity calculation |
| Tensors | `tensors.cpython-39-x86_64-linux-gnu.so` | Solvation kernel projection |

If using a different Python version or platform, recompile from source.

---

## Solvation Kernels

The docking module requires pre-computed solvation kernels:

```
solv_krnls_0.50A_24vxl.npz
```

This file contains voxelized solvation representations for overlap calculations.

---

## Verifying Installation

### Test Import

```python
python -c "
import numpy as np
import scipy
from joblib import Parallel, delayed
import sim
print('All imports successful!')
print(f'NumPy: {np.__version__}')
print(f'SciPy: {scipy.__version__}')
"
```

### Test Mapper

```bash
cd /path/to/HECTOR/code
python hector_mapper.py --help
```

Expected output:
```
usage: hector_mapper.py [-h] ply dist_x dist_y bin dot_skp map_skp infltn radius_fade slope_fade sign
```

---

## Directory Setup

Ensure the following directory structure:

```bash
# Create necessary directories (if not already present)
mkdir -p data/scaffolds_db
mkdir -p results
mkdir -p input
mkdir -p output
```

---

## Scaffold Database Setup

**Important:** Scaffold database files (*.pdb, *.npz) are not included in the git repository due to their size. You need to download and generate them.

### Option 1: Automated Setup (Recommended)

Use the provided setup scripts to automatically download and process example scaffolds:

#### Using Bash (Linux/Mac/Docker):
```bash
./scripts/setup_scaffold_database.sh
```

#### Using Python (Cross-platform):
```bash
python scripts/setup_scaffold_database.py
```

**What the scripts do:**
1. Download 8 example PDB structures from RCSB PDB
2. Extract specific chains (1CGI_A, 1CGI_B, 1ky2_A, 1uzi_A, 5djl_A, 5nlc_A, 7z64_B, 8brb_B)
3. Generate surface meshes using EDTSurf
4. Generate fingerprints using hector_mapper_parallel.py

**Processing time:** ~30-40 minutes for all scaffolds

**Expected output:**
```
data/scaffolds_db/
├── 1CGI_A.pdb          # PDB structure
├── 1CGI_A_rcpt.npz     # Pre-computed fingerprints
├── 1CGI_B.pdb
├── 1CGI_B_rcpt.npz
└── ... (16 files total)
```

### Option 2: Manual Setup for Custom Scaffolds

To add your own scaffolds to the database:

```bash
# 1. Download PDB structure
wget https://files.rcsb.org/download/XXXX.pdb -O data/scaffolds_db/XXXX_A.pdb

# 2. Generate surface mesh
EDTSurf -i data/scaffolds_db/XXXX_A.pdb -o output/XXXX_A.ply -s 3

# 3. Generate fingerprints (use parallel version for speed)
python code/hector_mapper_parallel.py output/XXXX_A.ply 10 20 0.2 40 8 0.5 40 0.3 rcpt

# 4. Move fingerprints to scaffold database
mv results/XXXX_A_rcpt.npz data/scaffolds_db/
```

---

## Troubleshooting

### Common Issues

#### 1. ImportError: No module named 'sim'

**Cause**: Cython module not compiled or wrong Python version.

**Solution**: Recompile with correct Python version:
```bash
cython sim.pyx
python setup.py build_ext --inplace
```

#### 2. Memory Error

**Cause**: Insufficient RAM for large databases.

**Solution**:
- Reduce `n_cores` in parallel processing
- Process scaffolds in batches
- Increase swap space

#### 3. EDTSurf Segmentation Fault

**Cause**: Invalid PDB format or missing atoms.

**Solution**:
- Clean PDB file (remove HETATM, fix numbering)
- Use a PDB preparation tool (e.g., PyMOL, pdbfixer)

#### 4. GPU Not Detected

**Cause**: Missing CUDA libraries or ArrayFire.

**Solution**:
- Install NVIDIA CUDA toolkit
- Install ArrayFire with CUDA backend
- Set `CUDA_VISIBLE_DEVICES` environment variable

---

## Performance Optimization

### CPU Optimization

```python
# Adjust number of parallel cores
n_cores = 10  # Set based on available CPU cores
```

### GPU Optimization (Optional)

For GPU-accelerated R-factor calculations:

1. Install ArrayFire with CUDA backend
2. Use GPU implementation from published code

---

## Next Steps

After installation, proceed to:
1. [Algorithm Description](ALGORITHM.md) - Understand the HECTOR methodology
2. [Pipeline Usage](PIPELINE.md) - Run your first analysis
