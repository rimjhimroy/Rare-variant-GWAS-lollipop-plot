#!/bin/bash

# Example usage script for Rare Variant GWAS Lollipop Plot Generator
# This script demonstrates various ways to use the tool

# Make sure you have activated the conda environment first:
# conda activate lollipop

echo "=== Rare Variant GWAS Lollipop Plot Generator - Examples ==="
echo ""

# Example 1: Basic usage with gene symbol
echo "Example 1: Basic usage with gene symbol GCDH"
echo "Command: Rscript lollipop_maker.R --symbol GCDH"
echo ""
Rscript lollipop_maker.R --symbol GCDH
echo ""

# Example 2: Using a different q-value threshold
echo "Example 2: Using a stricter q-value threshold (1e-6)"
echo "Command: Rscript lollipop_maker.R --symbol GCDH --qvalue 1e-6 --output plot/strict/"
echo ""
# Rscript lollipop_maker.R --symbol GCDH --qvalue 1e-6 --output plot/strict/
echo ""

# Example 3: Custom plot dimensions for publication
echo "Example 3: Custom plot dimensions (publication-ready)"
echo "Command: Rscript lollipop_maker.R --symbol GCDH --width 25 --height 10 --dpi 600 --output plot/publication/"
echo ""
# Rscript lollipop_maker.R --symbol GCDH --width 25 --height 10 --dpi 600 --output plot/publication/
echo ""

# Example 4: Using an Ensembl transcript ID
echo "Example 4: Using Ensembl transcript ID (if you know it)"
echo "Command: Rscript lollipop_maker.R --enst ENST00000123456"
echo "(Note: Replace with actual ENST ID from your data)"
echo ""

# Example 5: Custom input file paths
echo "Example 5: Custom input file paths"
echo "Command: Rscript lollipop_maker.R --symbol GCDH \\"
echo "  --variants /path/to/variants.rds \\"
echo "  --annotations /path/to/annotations.rds \\"
echo "  --transcripts /path/to/transcripts.tsv.gz \\"
echo "  --exons /path/to/exons.txt.gz \\"
echo "  --output results/"
echo ""

echo "=== Examples complete ==="
echo "Check the plot/ directory for output files"
