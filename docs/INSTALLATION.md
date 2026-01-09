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
cd /path/to/HECTOR/environment
docker build -t hector:latest .
```

#### 2. Run HECTOR in Docker

```bash
docker run -it -v /path/to/HECTOR:/workspace hector:latest
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
# Create necessary directories
mkdir -p /path/to/HECTOR/data
mkdir -p /path/to/HECTOR/results
mkdir -p /path/to/HECTOR/data/scaffolds_db
```

---

## Scaffold Database Setup

### Option 1: Use Provided Database

The repository includes a sample scaffold database with pre-computed fingerprints:
- `1CGI_A_rcpt.npz`, `1CGI_B_rcpt.npz`
- `5djl_A_rcpt.npz`, `5nlc_A_rcpt.npz`
- `7z64_B_rcpt.npz`, `8brb_B_rcpt.npz`

### Option 2: Build Custom Database

To create fingerprints for new scaffolds:

```bash
# 1. Download PDB structures
wget https://files.rcsb.org/download/XXXX.pdb

# 2. Generate surface mesh
EDTSurf -i XXXX.pdb -o XXXX.ply -s 3

# 3. Generate fingerprints
python hector_mapper.py XXXX.ply 10 20 0.2 5 1 0.5 40 0.3 rcpt
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
