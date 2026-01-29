#!/bin/bash
# Setup script for HECTOR scaffold database
# Downloads PDB files and generates fingerprints for example scaffolds

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}HECTOR Scaffold Database Setup${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""

# Create directories
echo -e "${YELLOW}Creating directories...${NC}"
mkdir -p data/scaffolds_db
mkdir -p output/scaffolds
mkdir -p results

# Example scaffolds to download
# Format: "PDB_ID:CHAIN"
SCAFFOLDS=(
    "1CGI:A"
    "1CGI:B"
    "1ky2:A"
    "1uzi:A"
    "5djl:A"
    "5nlc:A"
    "7z64:B"
    "8brb:B"
)

# Check if EDTSurf is available
if ! command -v EDTSurf &> /dev/null; then
    echo -e "${RED}Error: EDTSurf not found. Please install EDTSurf first.${NC}"
    exit 1
fi

# Check if Python is available
if ! command -v python &> /dev/null; then
    echo -e "${RED}Error: Python not found.${NC}"
    exit 1
fi

echo -e "${GREEN}Downloading and processing scaffolds...${NC}"
echo ""

for scaffold in "${SCAFFOLDS[@]}"; do
    IFS=':' read -r pdb_id chain <<< "$scaffold"
    pdb_lower=$(echo "$pdb_id" | tr '[:upper:]' '[:lower:]')
    output_name="${pdb_id}_${chain}"

    echo -e "${YELLOW}Processing ${output_name}...${NC}"

    # Download PDB file
    if [ ! -f "data/scaffolds_db/${output_name}.pdb" ]; then
        echo "  Downloading ${pdb_id}..."
        wget -q "https://files.rcsb.org/download/${pdb_id}.pdb" -O "/tmp/${pdb_id}.pdb"

        if [ ! -f "/tmp/${pdb_id}.pdb" ]; then
            echo -e "  ${RED}Failed to download ${pdb_id}${NC}"
            continue
        fi

        # Extract specific chain
        echo "  Extracting chain ${chain}..."
        grep "^ATOM.*\ ${chain}\ " "/tmp/${pdb_id}.pdb" > "data/scaffolds_db/${output_name}.pdb" || true
        grep "^HETATM.*\ ${chain}\ " "/tmp/${pdb_id}.pdb" >> "data/scaffolds_db/${output_name}.pdb" || true
        echo "END" >> "data/scaffolds_db/${output_name}.pdb"

        rm "/tmp/${pdb_id}.pdb"
    else
        echo "  PDB file already exists, skipping download"
    fi

    # Generate surface with EDTSurf
    if [ ! -f "output/scaffolds/${output_name}.ply" ]; then
        echo "  Generating surface mesh..."
        EDTSurf -i "data/scaffolds_db/${output_name}.pdb" \
                -o "output/scaffolds/${output_name}.ply" \
                -s 3 > /dev/null 2>&1
    else
        echo "  Surface mesh already exists, skipping"
    fi

    # Generate fingerprints with parallel version
    if [ ! -f "data/scaffolds_db/${output_name}_rcpt.npz" ]; then
        echo "  Generating fingerprints (this may take a few minutes)..."

        # Use parallel version if available, otherwise use standard version
        if [ -f "code/hector_mapper_parallel.py" ]; then
            python code/hector_mapper_parallel.py \
                "output/scaffolds/${output_name}.ply" \
                10 20 0.2 40 8 0.5 40 0.3 rcpt \
                --n_jobs -1 > /dev/null 2>&1
        else
            python code/hector_mapper.py \
                "output/scaffolds/${output_name}.ply" \
                10 20 0.2 40 8 0.5 40 0.3 rcpt > /dev/null 2>&1
        fi

        # Move fingerprints to scaffolds_db
        mv "results/${output_name}_rcpt.npz" "data/scaffolds_db/${output_name}_rcpt.npz"
    else
        echo "  Fingerprints already exist, skipping"
    fi

    echo -e "  ${GREEN}✓ ${output_name} complete${NC}"
    echo ""
done

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}Scaffold database setup complete!${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo "Summary:"
echo "  Location: data/scaffolds_db/"
ls -lh data/scaffolds_db/ | grep -E "\\.pdb|\\.npz" | wc -l | xargs echo "  Files created:"
echo ""
echo "You can now run the HECTOR pipeline with:"
echo "  python code/maps_analysis_pairs_vs_all.py data/scaffolds_db/ results/query_rcpt.npz -0.82"
