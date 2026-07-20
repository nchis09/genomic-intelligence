#!/bin/bash
# Download Kraken2 viral database (~2.5GB)
# This script downloads a viral-specific database for Kraken2 classification

set -e

DB_NAME="kraken2_viral"
DB_DIR="databases"

echo "=========================================="
echo "Kraken2 Viral Database Download"
echo "=========================================="
echo ""

# Create databases directory if it doesn't exist
mkdir -p $DB_DIR

# Check if kraken2-build is available
if ! command -v kraken2-build &> /dev/null; then
    echo "ERROR: kraken2-build not found"
    echo ""
    echo "To install Kraken2:"
    echo "  conda install -c bioconda kraken2"
    echo "  OR"
    echo "  Download from: https://github.com/DerrickWood/kraken2"
    echo ""
    exit 1
fi

echo "Step 1: Downloading viral library from NCBI (~2.5GB)..."
kraken2-build --download-library viral --db $DB_DIR/$DB_NAME

echo ""
echo "Step 2: Building the database..."
kraken2-build --build --db $DB_DIR/$DB_NAME

echo ""
echo "=========================================="
echo "Database download complete!"
echo "=========================================="
echo "Database location: $DB_DIR/$DB_NAME"
echo "Database size: $(du -sh $DB_DIR/$DB_NAME | cut -f1)"
echo ""
echo "To use with Stage 1:"
echo "  python stage1_classification.py --fasta input_FASTA.fasta \\"
echo "    --taxonomy reference_library/ebola/taxonomy/taxonomy.yaml \\"
echo "    --kraken-db $DB_DIR/$DB_NAME \\"
echo "    --output output_dir"
echo ""
