# Quick Reference Guide

## Quick Start

1. **Install dependencies:**
   ```bash
   conda env create -f environment.yml
   conda activate lollipop
   ```

2. **Run with default settings:**
   ```bash
   Rscript lollipop_maker.R --symbol GCDH
   ```

3. **Check output:**
   ```bash
   ls plot/GCDH_lolli.*
   ```

## Common Commands

### Generate plot for a gene
```bash
Rscript lollipop_maker.R --symbol GCDH
```

### Use custom q-value threshold
```bash
Rscript lollipop_maker.R --symbol GCDH --qvalue 1e-6
```

### Specify custom output directory
```bash
Rscript lollipop_maker.R --symbol GCDH --output results/
```

### High-resolution plot for publication
```bash
Rscript lollipop_maker.R --symbol GCDH --width 25 --height 10 --dpi 600
```

### Use custom input files
```bash
Rscript lollipop_maker.R \
  --symbol GCDH \
  --variants my_data/variants.rds \
  --annotations my_data/annotations.rds
```

## File Locations

- **Script:** `lollipop_maker.R`
- **Helper scripts:** `scripts/make_gene_model.R`, `scripts/plot_lollipop.R`, `scripts/pfam_dom.R`
- **Input data:** `data/` directory
- **Output:** `plot/` directory (default)

## Troubleshooting

### "Gene not found in variants data"
- Check that the gene symbol or ENST ID exists in your variants.rds file
- The gene name is case-sensitive

### "File not found" errors
- Verify all input files exist in the specified paths
- Use absolute paths if relative paths don't work

### "Multiple transcripts found"
- The tool will use the first transcript found
- For precise control, use `--enst` with a specific transcript ID

### Missing Pfam domains
- Ensure you have internet access (required for biomaRt queries)
- The tool connects to Ensembl to fetch protein domain information

## Input File Format Details

### variants.rds
R named list where names are gene symbols/ENST IDs:
```r
variants <- list(
  "GCDH" = data.frame(
    CHROM = "19",
    GENPOS = 12896252,
    ID = "19:12896252:G>A",
    # ... other columns
  ),
  # ... more genes
)
```

### variant_annotations.rds
R named list matching variants.rds structure:
```r
annotations <- list(
  "GCDH" = data.frame(
    variant = "19:12896252:G>A",
    consequence = "missense_variant"
  ),
  # ... more genes
)
```

## Dependencies

All dependencies are managed through the conda environment (`environment.yml`). Key packages:
- R 4.1.3
- biomaRt (for Ensembl queries)
- ggplot2 (for plotting)
- data.table (for data manipulation)
- ensembldb + EnsDb.Hsapiens.v86 (for annotations)

## Getting Help

Display all command-line options:
```bash
Rscript lollipop_maker.R --help
```
