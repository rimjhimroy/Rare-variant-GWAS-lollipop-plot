# Rare Variant GWAS Lollipop Plot Generator

A visualization tool for creating publication-ready lollipop plots of rare variant associations per gene. The plots display genetic variants with their association statistics, protein domains (Pfam), exon structure, allele frequencies, and effect sizes (Beta values).

## Features

- **Comprehensive Visualization**: Shows variants with their genomic positions, association q-values, effect directions (positive/negative Beta), and allele frequencies
- **Protein Domain Annotation**: Integrates Pfam domain information from Ensembl
- **Gene Structure**: Displays exon/intron structure with UTRs and coding sequences (CDS)
- **Customizable**: Configurable q-value thresholds, output formats, and plot dimensions
- **Flexible Input**: Works with Regenie association studies and can be adapted to other GWAS tools

## Installation

### Using Conda (Recommended)

Create and activate the conda environment with all required dependencies:

```bash
conda env create -f environment.yml
conda activate lollipop
```

### Manual Installation

If you prefer to install dependencies manually, you'll need R (version 4.1+) with the following packages:

**Bioconductor packages:**
- biomaRt
- ensembldb
- EnsDb.Hsapiens.v86
- (see environment.yml for complete list)

**CRAN packages:**
- data.table
- ggplot2
- tidyverse
- ggrepel
- lemon
- stringr
- bedr
- optparse

## Input Data Requirements

The tool requires the following input files:

### 1. Variant Association Results (RDS format)

A named list of data frames, where each element corresponds to a gene. Each data frame must contain:

- `CHROM`: Chromosome
- `GENPOS`: Genomic position
- `ID`: Variant identifier
- `ALLELE0`: Reference allele
- `ALLELE1`: Alternate allele
- `A1FREQ`: Alternate allele frequency
- `N`: Sample size
- `TEST`: Test name
- `BETA`: Effect size (Beta coefficient)
- `SE`: Standard error
- `trait`: Trait name
- `qval`: Q-value (multiple-testing corrected p-value)
- `variant`: Variant in format "chr:pos:ref>alt" (e.g., "19:12896252:G>A")

### 2. Variant Annotations (RDS format)

A named list of data frames with variant functional consequences:

- `variant`: Variant identifier (matching the variants file)
- `consequence`: Variant consequence (e.g., "missense_variant", "synonymous_variant", "stop_gained")

### 3. Transcript Information (TSV.GZ format)

Tab-delimited file with transcript annotations containing at minimum:
- `ENST`: Ensembl transcript ID
- `SYMBOL`: Gene symbol
- (Additional columns as needed)

### 4. Exon Models (TXT.GZ format)

Tab-delimited file with exon structure information:
- Transcript ID
- CDS start
- CDS end
- Number of exons
- Exon starts (comma-separated)
- Exon ends (comma-separated)
- Transcript type

## Usage

### Basic Usage

Generate a lollipop plot for a specific gene:

```bash
Rscript lollipop_maker.R --symbol GCDH
```

Or using an Ensembl transcript ID:

```bash
Rscript lollipop_maker.R --enst ENST00000123456
```

### Advanced Usage

Customize input files, thresholds, and output:

```bash
Rscript lollipop_maker.R \
  --symbol GCDH \
  --variants data/my_variants.rds \
  --annotations data/my_annotations.rds \
  --transcripts data/my_transcripts.tsv.gz \
  --exons data/my_exons.txt.gz \
  --qvalue 1e-5 \
  --output results/ \
  --width 25 \
  --height 10 \
  --dpi 600
```

### Command-Line Options

| Option | Short | Description | Default |
|--------|-------|-------------|---------|
| `--symbol` | `-s` | Gene symbol (e.g., GCDH) | - |
| `--enst` | `-e` | Ensembl transcript ID | - |
| `--variants` | `-v` | Path to variants RDS file | `data/variants.rds` |
| `--annotations` | `-a` | Path to variant annotations RDS file | `data/variant_annotations.rds` |
| `--transcripts` | `-t` | Path to transcripts file | `data/transcripts.tsv.gz` |
| `--exons` | `-x` | Path to exon models file | `data/exon_models.txt.gz` |
| `--qvalue` | `-q` | Q-value threshold for labeling | `5.0e-8` |
| `--output` | `-o` | Output directory | `plot` |
| `--width` | `-w` | Plot width in inches | `20` |
| `--height` | `-H` | Plot height in inches | `8` |
| `--dpi` | `-d` | Plot resolution (DPI) | `450` |

**Note:** Either `--symbol` or `--enst` must be provided.

### Getting Help

```bash
Rscript lollipop_maker.R --help
```

## Output Files

The tool generates three output files per gene:

1. **`{SYMBOL}_lolli.png`**: High-resolution PNG plot
2. **`{SYMBOL}_lolli.tiff`**: TIFF format plot (for publications)
3. **`{SYMBOL}_lolli.tsv`**: Tab-separated table of variants with annotations and statistics

### Plot Interpretation

- **Lollipop stems**: Represent individual variants at their genomic positions
- **Lollipop heads (circles)**: Size indicates allele frequency
- **Y-axis position**: Height shows -log10(q-value); positive/negative indicates Beta direction
  - Upper half: Positive Beta (risk-increasing)
  - Lower half: Negative Beta (protective)
- **Red dashed line**: Q-value significance threshold
- **Blue rectangles**: Coding sequence (CDS) exons
- **Thin rectangles**: Untranslated regions (UTRs)
- **Colored rectangles at bottom**: Pfam protein domains
- **Labels**: Significant variants (below q-value threshold) are labeled with ID, frequency, trait, and consequence

## Example

Using the provided example data:

```bash
# Activate conda environment
conda activate lollipop

# Generate plot for GCDH gene
Rscript lollipop_maker.R --symbol GCDH

# Output will be created in plot/ directory:
# - plot/GCDH_lolli.png
# - plot/GCDH_lolli.tiff
# - plot/GCDH_lolli.tsv
```

## Adapting for Other GWAS Tools

While designed for Regenie output, the tool can be adapted for other GWAS tools by:

1. Converting your association results to the required format (see Input Data Requirements)
2. Ensuring variant IDs match between association results and annotations
3. Computing q-values from p-values if not provided (e.g., using the `qvalue` R package)

## Citation

If you use this tool in your research, please cite:

R.R. Choudhury (2025). Rare variant GWAS lollipop plots. https://github.com/rimjhimroy/Rare-variant-GWAS-lollipop-plot

Suggested BibTeX for this wrapper:

```bibtex
@misc{rchoudhury_lollipop_plot_2025,
	author = {Choudhury, R. R.},
	title = {Rare variant GWAS lollipop plots},
	year = {2025},
	howpublished = {Repository / workflow in project},
	note = {URL: https://github.com/rimjhimroy/Rare-variant-GWAS-lollipop-plot}
}
```

## License

MIT License - see LICENSE file for details

## Support

For issues, questions, or contributions, please visit:
https://github.com/rimjhimroy/Rare-variant-GWAS-lollipop-plot

## Acknowledgments

This tool integrates data from:
- Ensembl (biomaRt)
- Pfam protein domain database
- EnsDb.Hsapiens.v86 annotation package
