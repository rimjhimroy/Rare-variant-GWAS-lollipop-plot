## create and load conda environment with environment.yml
## requires Gene symbol in "symbol" and results of association tests with the following column: "CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "TEST", "BETA", "SE", "trait", "qval", "variant" -> (eg. 19:12896252:G>A)
## also needs variant annotation file with variant id and consequence, i.e. missense, synonymous, etc.
## Quietly Load Libraries
load.package <- function(name) {
    suppressMessages(suppressWarnings(library(name, quietly = T, warn.conflicts = F, character.only = T)))
}
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

source("scripts/make_gene_model.R")
source("scripts/plot_lollipop.R")
source("scripts/pfam_dom.R")

enst <- NA
symbol <- "GCDH"

res.df <- readRDS("data/variants.rds")
variant_annotations <- readRDS("data/variant_annotations.rds")
dat <- res.df[[symbol]]
vannot <- variant_annotations[[symbol]]
dat1 <- dplyr::left_join(dat, vannot[, c("variant", "consequence")], by = c("ID" = "variant"))
variants <- data.table(dat1)
qvalue <- 5.0e-8 # qvalue threshold to plot (genomewide pvalue threshold used)

## Load transcripts and get info for gene of interest
transcripts <- fread("data/transcripts.tsv.gz")

if (!is.na(enst)) {
    gene_info <- transcripts[ENST == enst]
} else if (!is.na(symbol)) {
    gene_info <- transcripts[SYMBOL == symbol]
    if (nrow(gene_info) > 1) {
        stop(paste0("Found more than one gene with symbol ", symbol, ". Please try again, possibly with ENST."))
    }
} else {
    stop("Neither gene ENST (-e/--enst) nor Symbol (-s/--symbol) were provided. Please try again!")
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

gene_model <- make_gene_model(gene_info, variants)
pfam_dom <- make_pfam_dom(gene_model, symbol)
lollipop_plot <- plot_lollipop(gene_model, pfam_dom, qvalue)

file_prefix <- paste0(symbol, "_lolli")
tiff(filename = paste0("plot/", file_prefix, ".tiff"), , units = "in", width = 15, height = 8, res = 450)
lollipop_plot
dev.off()
ggsave(filename = paste0("plot/", file_prefix, ".png"), plot = lollipop_plot, width = 20, height = 8, dpi = 450)
out_table=gene_model$variants[order(-log.q),]
fwrite(out_table, file = paste0("plot/",file_prefix,".tsv"), col.names = T, row.names = F, quote = F, na = "NA", sep = "\t")
