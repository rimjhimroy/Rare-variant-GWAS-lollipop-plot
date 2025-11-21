#!/usr/bin/env Rscript

## Rare variant GWAS Lollipop Plot Generator
## Creates lollipop plots per gene with exon, pfam domain, qvalue, allele frequency and Beta information
## 
## Requirements:
## - Gene symbol or ENST transcript ID
## - Association test results with columns: "CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "TEST", "BETA", "SE", "trait", "qval", "variant"
## - Variant annotation file with variant id and consequence (e.g., missense, synonymous)
## - Transcript file (transcripts.tsv.gz)
## - Exon models file (exon_models.txt.gz)
##
## Usage: Rscript lollipop_maker.R --symbol GCDH --variants data/variants.rds --annotations data/variant_annotations.rds --output plot/

## Parse command line arguments
suppressMessages(suppressWarnings(library("optparse", quietly = TRUE, warn.conflicts = FALSE)))

option_list <- list(
  make_option(c("-s", "--symbol"), type = "character", default = NULL,
              help = "Gene symbol (e.g., GCDH)", metavar = "character"),
  make_option(c("-e", "--enst"), type = "character", default = NULL,
              help = "Ensembl transcript ID (e.g., ENST00000123456)", metavar = "character"),
  make_option(c("-v", "--variants"), type = "character", default = "data/variants.rds",
              help = "Path to variants RDS file [default: %default]", metavar = "character"),
  make_option(c("-a", "--annotations"), type = "character", default = "data/variant_annotations.rds",
              help = "Path to variant annotations RDS file [default: %default]", metavar = "character"),
  make_option(c("-t", "--transcripts"), type = "character", default = "data/transcripts.tsv.gz",
              help = "Path to transcripts file [default: %default]", metavar = "character"),
  make_option(c("-x", "--exons"), type = "character", default = "data/exon_models.txt.gz",
              help = "Path to exon models file [default: %default]", metavar = "character"),
  make_option(c("-q", "--qvalue"), type = "double", default = 5.0e-8,
              help = "Q-value threshold for plotting [default: %default]", metavar = "number"),
  make_option(c("-o", "--output"), type = "character", default = "plot",
              help = "Output directory for plots and tables [default: %default]", metavar = "character"),
  make_option(c("-w", "--width"), type = "integer", default = 20,
              help = "Plot width in inches [default: %default]", metavar = "number"),
  make_option(c("-H", "--height"), type = "integer", default = 8,
              help = "Plot height in inches [default: %default]", metavar = "number"),
  make_option(c("-d", "--dpi"), type = "integer", default = 450,
              help = "Plot resolution (DPI) [default: %default]", metavar = "number")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

## Validate required arguments
if (is.null(opt$symbol) && is.null(opt$enst)) {
  print_help(opt_parser)
  stop("Either --symbol or --enst must be provided", call. = FALSE)
}

## Quietly Load Libraries
load.package <- function(name) {
    suppressMessages(suppressWarnings(library(name, quietly = T, warn.conflicts = F, character.only = T)))
}

cat("Loading required libraries...\n")
load.package("biomaRt")
load.package("data.table")
load.package("ggplot2")
load.package("bedr")
load.package("lemon")
load.package("stringr")
load.package("ggrepel")
load.package("tidyverse")
load.package("ensembldb")
load.package("EnsDb.Hsapiens.v86")

## Get script directory for sourcing helper scripts
script_dir <- dirname(sys.frame(1)$ofile)
if (length(script_dir) == 0 || script_dir == "") {
  script_dir <- "."
}

source(file.path(script_dir, "scripts/make_gene_model.R"))
source(file.path(script_dir, "scripts/plot_lollipop.R"))
source(file.path(script_dir, "scripts/pfam_dom.R"))

## Set parameters from command line arguments
enst <- opt$enst
symbol <- opt$symbol
qvalue <- opt$qvalue

## Create output directory if it doesn't exist
if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
  cat("Created output directory:", opt$output, "\n")
}

## Load input files
cat("Loading variant data from:", opt$variants, "\n")
if (!file.exists(opt$variants)) {
  stop("Variants file not found: ", opt$variants, call. = FALSE)
}
res.df <- readRDS(opt$variants)

cat("Loading variant annotations from:", opt$annotations, "\n")
if (!file.exists(opt$annotations)) {
  stop("Variant annotations file not found: ", opt$annotations, call. = FALSE)
}
variant_annotations <- readRDS(opt$annotations)

## Determine gene symbol for data extraction
gene_name <- if (!is.null(symbol)) symbol else enst

## Extract variant data for the gene
if (!gene_name %in% names(res.df)) {
  stop("Gene '", gene_name, "' not found in variants data. Available genes: ", 
       paste(head(names(res.df), 10), collapse = ", "), "...", call. = FALSE)
}
dat <- res.df[[gene_name]]

if (!gene_name %in% names(variant_annotations)) {
  stop("Gene '", gene_name, "' not found in variant annotations. Available genes: ", 
       paste(head(names(variant_annotations), 10), collapse = ", "), "...", call. = FALSE)
}
vannot <- variant_annotations[[gene_name]]

## Merge variant data with annotations
dat1 <- dplyr::left_join(dat, vannot[, c("variant", "consequence")], by = c("ID" = "variant"))
variants <- data.table(dat1)

cat("Processing", nrow(variants), "variants for gene:", gene_name, "\n")

## Load transcripts and get info for gene of interest
cat("Loading transcript data from:", opt$transcripts, "\n")
if (!file.exists(opt$transcripts)) {
  stop("Transcripts file not found: ", opt$transcripts, call. = FALSE)
}
transcripts <- fread(opt$transcripts)

## Get gene info from transcripts
if (!is.null(enst)) {
    cat("Looking up transcript:", enst, "\n")
    gene_info <- transcripts[ENST == enst]
    if (nrow(gene_info) == 0) {
      stop("Transcript '", enst, "' not found in transcripts file", call. = FALSE)
    }
    # Get symbol from ENST if not provided
    if (is.null(symbol)) {
      symbol <- gene_info[1, SYMBOL]
      cat("Using gene symbol from transcript:", symbol, "\n")
    }
} else if (!is.null(symbol)) {
    cat("Looking up gene symbol:", symbol, "\n")
    gene_info <- transcripts[SYMBOL == symbol]
    if (nrow(gene_info) == 0) {
      stop("Gene symbol '", symbol, "' not found in transcripts file", call. = FALSE)
    }
    if (nrow(gene_info) > 1) {
        warning(paste0("Found ", nrow(gene_info), " transcripts for gene symbol ", symbol, ". Using first one."))
        gene_info <- gene_info[1, ]
    }
}

## Alternatively, get gene info and canonical transcript ID from biomart
# ensembl <- useMart("ensembl")
# ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
# gene_info <- getBM(
#     attributes = c("chromosome_name", "start_position", "end_position", "strand", "ensembl_transcript_id", "ensembl_gene_id", "transcript_mane_select", "transcript_length", "hgnc_symbol", "transcript_is_canonical", "transcript_biotype"),
#     filters = c("hgnc_symbol"),
#     values = list(genename),
#     mart = ensembl
# )
# gene_info <- data.table(gene_info[which(gene_info$transcript_is_canonical == 1), ])
# names(gene_info) <- c("#chrom", "start", "end", "strand", "ENST", "ENSG", "MANE", "transcript_length", "SYMBOL", "CANONICAL", "BIOTYPE")

cat("Building gene model...\n")
gene_model <- make_gene_model(gene_info, variants, opt$exons)

cat("Fetching protein domains...\n")
pfam_dom <- make_pfam_dom(gene_model, symbol)

cat("Creating lollipop plot...\n")
lollipop_plot <- plot_lollipop(gene_model, pfam_dom, qvalue, symbol)

## Save outputs
file_prefix <- paste0(symbol, "_lolli")
output_base <- file.path(opt$output, file_prefix)

cat("Saving plot to:", output_base, "\n")

# Save as TIFF
tiff(filename = paste0(output_base, ".tiff"), units = "in", width = opt$width, height = opt$height, res = opt$dpi)
print(lollipop_plot)
dev.off()

# Save as PNG
ggsave(filename = paste0(output_base, ".png"), plot = lollipop_plot, width = opt$width, height = opt$height, dpi = opt$dpi)

# Save variant table
out_table <- gene_model$variants[order(-log.q), ]
fwrite(out_table, file = paste0(output_base, ".tsv"), col.names = TRUE, row.names = FALSE, quote = FALSE, na = "NA", sep = "\t")

cat("\nDone! Output files:\n")
cat("  - Plot (PNG):", paste0(output_base, ".png"), "\n")
cat("  - Plot (TIFF):", paste0(output_base, ".tiff"), "\n")
cat("  - Variant table:", paste0(output_base, ".tsv"), "\n")
