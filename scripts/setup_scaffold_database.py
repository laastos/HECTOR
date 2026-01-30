#!/usr/bin/env python
"""
Setup script for HECTOR scaffold database
Downloads PDB files and generates fingerprints for example scaffolds
"""

import os
import sys
import subprocess
import urllib.request
from pathlib import Path

# Example scaffolds to download
SCAFFOLDS = [
    ("1CGI", "E"),  # 1CGI has chains E and I, not A and B
    ("1CGI", "I"),
    ("1ky2", "A"),
    ("1uzi", "A"),
    ("5djl", "A"),
    ("5nlc", "A"),
    ("7z64", "B"),
    ("8brb", "B"),
]

# ANSI color codes
GREEN = '\033[0;32m'
YELLOW = '\033[1;33m'
RED = '\033[0;31m'
NC = '\033[0m'  # No Color


def print_color(message, color=NC):
    """Print colored message"""
    print(f"{color}{message}{NC}")


def check_dependencies():
    """Check if required tools are available"""
    # Check EDTSurf
    try:
        subprocess.run(['EDTSurf'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except FileNotFoundError:
        print_color("Error: EDTSurf not found. Please install EDTSurf first.", RED)
        sys.exit(1)

    # Check Python modules
    try:
        import numpy
    except ImportError:
        print_color("Error: NumPy not found. Please install required Python packages.", RED)
        sys.exit(1)


def download_pdb(pdb_id, chain, output_dir):
    """Download PDB file and extract specific chain"""
    output_name = f"{pdb_id}_{chain}"
    output_path = output_dir / f"{output_name}.pdb"

    if output_path.exists():
        print(f"  PDB file already exists, skipping download")
        return output_path

    print(f"  Downloading {pdb_id}...")
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

    try:
        with urllib.request.urlopen(url) as response:
            pdb_content = response.read().decode('utf-8')
    except Exception as e:
        print_color(f"  Failed to download {pdb_id}: {e}", RED)
        return None

    # Extract specific chain
    print(f"  Extracting chain {chain}...")
    chain_lines = []
    header_lines = []

    # First pass: collect header records
    for line in pdb_content.split('\n'):
        if line.startswith(('HEADER', 'TITLE', 'COMPND', 'SOURCE', 'KEYWDS',
                           'EXPDTA', 'AUTHOR', 'REVDAT', 'JRNL', 'REMARK',
                           'DBREF', 'SEQRES', 'CRYST1', 'ORIGX', 'SCALE', 'MTRIX')):
            header_lines.append(line)

    # Second pass: collect ATOM/HETATM records for specified chain
    for line in pdb_content.split('\n'):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            # Chain ID is at column 21 (0-indexed) in standard PDB format
            if len(line) > 21 and line[21] == chain:
                # Truncate to column 66 to remove non-standard extra fields
                # Standard PDB format is columns 1-66 (ATOM to B-factor)
                clean_line = line[:66].rstrip()
                chain_lines.append(clean_line)

    # Write PDB file with proper structure
    with open(output_path, 'w') as f:
        # Write headers
        if header_lines:
            f.write('\n'.join(header_lines))
            f.write('\n')

        # Write atoms
        if chain_lines:
            f.write('\n'.join(chain_lines))
            f.write('\n')

        # Write END record
        f.write('TER\n')
        f.write('END\n')

    return output_path


def generate_surface(pdb_path, output_dir):
    """Generate surface mesh using EDTSurf"""
    # EDTSurf automatically adds .ply extension, so don't include it in the path
    ply_base = output_dir / pdb_path.stem
    ply_path = output_dir / f"{pdb_path.stem}.ply"

    if ply_path.exists():
        print("  Surface mesh already exists, skipping")
        return ply_path

    print("  Generating surface mesh...")
    # Note: EDTSurf returns exit code 1 even on success, so we can't use check=True
    subprocess.run([
        'EDTSurf',
        '-i', str(pdb_path),
        '-o', str(ply_base),  # Don't include .ply - EDTSurf adds it automatically
        '-s', '3'
    ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Verify the output file was actually created
    if not ply_path.exists():
        print_color(f"  Failed to generate surface: output file not created", RED)
        return None

    return ply_path


def generate_fingerprints(ply_path, output_dir, use_parallel=True):
    """Generate fingerprints using hector_mapper"""
    npz_path = output_dir / f"{ply_path.stem}_rcpt.npz"

    if npz_path.exists():
        print("  Fingerprints already exist, skipping")
        return npz_path

    print("  Generating fingerprints (this may take a few minutes)...")

    # Choose mapper script
    if use_parallel and Path("code/hector_mapper_parallel.py").exists():
        mapper_script = "code/hector_mapper_parallel.py"
        extra_args = ["--n_jobs", "-1"]
    else:
        mapper_script = "code/hector_mapper.py"
        extra_args = []

    try:
        subprocess.run([
            sys.executable,
            mapper_script,
            str(ply_path),
            "10", "20", "0.2", "40", "8", "0.5", "40", "0.3", "rcpt"
        ] + extra_args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

        # Move from results/ to scaffolds_db/
        results_npz = Path("results") / f"{ply_path.stem}_rcpt.npz"
        if results_npz.exists():
            results_npz.rename(npz_path)
        else:
            print_color(f"  Warning: Expected output not found: {results_npz}", YELLOW)
            return None

    except subprocess.CalledProcessError as e:
        print_color(f"  Failed to generate fingerprints: {e}", RED)
        return None

    return npz_path


def main():
    """Main setup function"""
    print_color("=" * 40, GREEN)
    print_color("HECTOR Scaffold Database Setup", GREEN)
    print_color("=" * 40, GREEN)
    print()

    # Check dependencies
    check_dependencies()

    # Create directories
    print_color("Creating directories...", YELLOW)
    scaffolds_db = Path("data/scaffolds_db")
    scaffolds_ply = Path("output/scaffolds")
    results_dir = Path("results")

    scaffolds_db.mkdir(parents=True, exist_ok=True)
    scaffolds_ply.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    print_color("Downloading and processing scaffolds...", GREEN)
    print()

    success_count = 0
    total_count = len(SCAFFOLDS)

    for pdb_id, chain in SCAFFOLDS:
        output_name = f"{pdb_id}_{chain}"
        print_color(f"Processing {output_name}...", YELLOW)

        # Download PDB
        pdb_path = download_pdb(pdb_id, chain, scaffolds_db)
        if not pdb_path:
            print()
            continue

        # Generate surface
        ply_path = generate_surface(pdb_path, scaffolds_ply)
        if not ply_path:
            print()
            continue

        # Generate fingerprints
        npz_path = generate_fingerprints(ply_path, scaffolds_db)
        if not npz_path:
            print()
            continue

        print_color(f"  ✓ {output_name} complete", GREEN)
        print()
        success_count += 1

    print_color("=" * 40, GREEN)
    print_color("Scaffold database setup complete!", GREEN)
    print_color("=" * 40, GREEN)
    print()
    print(f"Summary:")
    print(f"  Successfully processed: {success_count}/{total_count} scaffolds")
    print(f"  Location: {scaffolds_db}")
    print()
    print("You can now run the HECTOR pipeline with:")
    print("  python code/maps_analysis_pairs_vs_all.py data/scaffolds_db/ results/query_rcpt.npz -0.82")


if __name__ == "__main__":
    main()
