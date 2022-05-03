## This function is to setup the appropriate gene model with "tiny" introns so that the reader can focus on coding sequence
## Then this function modifies the position of the variants according to the new gene model and returns the modified gene model and variants file.
## Input is the gene_info created in the main lollipop_maker.R script and a data.table of variants with qvalue "qval", assiciation betas "BETA", and position of the variants in the genome "GENPOS".

make_gene_model <- function(gene_info, variants) {
    # Load the gene's exon model:
    exon_models <- fread("data/exon_models.txt.gz")
    setnames(exon_models, names(exon_models), c("transcript", "cdsStart", "cdsEnd", "numExons", "exonStarts", "exonEnds", "transcriptType"))
    exon_models[, transcript := str_split(transcript, "\\.", simplify = T)[1], by = 1:nrow(exon_models)]
    ensembl_annotation_gene <- exon_models[transcript == gene_info[, ENST]]

    # Here we plot the actual lolliplot

    # Create initial exon breakpoints
    coding.start <- ensembl_annotation_gene[, cdsStart]
    coding.end <- ensembl_annotation_gene[, cdsEnd]

    exon.starts <- as.integer(unlist(str_split(ensembl_annotation_gene[, exonStarts], ",")))
    exon.starts <- exon.starts[1:length(exon.starts) - 1]
    exon.ends <- as.integer(unlist(str_split(ensembl_annotation_gene[, exonEnds], ",")))
    exon.ends <- exon.ends[1:length(exon.ends) - 1]

    pos.map <- data.table()
    gene.map <- data.table()
    last.pos <- 1
    for (i in 1:length(exon.starts)) {
        if (exon.starts[i] < coding.start && exon.ends[i] < coding.start) {
            # all UTR
            gene.map <- rbind(
                gene.map,
                data.table(
                    start = exon.starts[i],
                    end = exon.ends[i],
                    ymin = -0.25,
                    ymax = 0.25,
                    annotation = "utr"
                )
            )
        } else if (exon.starts[i] > coding.end && exon.ends[i] > coding.end) {
            gene.map <- rbind(
                gene.map,
                data.table(
                    start = exon.starts[i],
                    end = exon.ends[i],
                    ymin = -0.25,
                    ymax = 0.25,
                    annotation = "utr"
                )
            )
        } else if (exon.starts[i] < coding.start && exon.ends[i] > coding.end) {
            # Single exon gene <-  make three exons:
            gene.map <- rbind(
                gene.map,
                data.table(
                    start = c(exon.starts[i], coding.start, coding.end),
                    end = c(coding.start, coding.end, exon.ends[i]),
                    ymin = c(-0.25, -0.5, -0.25),
                    ymax = c(0.25, 0.5, 0.25),
                    annotation = c("utr", "cds", "utr")
                )
            )
        } else if (exon.starts[i] < coding.start && exon.ends[i] > coding.start) {
            # CDS + UTR <- make two exons
            gene.map <- rbind(
                gene.map,
                data.table(
                    start = c(exon.starts[i], coding.start),
                    end = c(coding.start, exon.ends[i]),
                    ymin = c(-0.25, -0.5),
                    ymax = c(0.25, 0.5),
                    annotation = c("utr", "cds")
                )
            )
        } else if (exon.starts[i] < coding.end && exon.ends[i] > coding.end) {
            # CDS + UTR <- make two exons
            # CDS + UTR <- make two exons
            gene.map <- rbind(
                gene.map,
                data.table(
                    start = c(exon.starts[i], coding.end),
                    end = c(coding.end, exon.ends[i]),
                    ymin = c(-0.5, -0.25),
                    ymax = c(0.5, 0.25),
                    annotation = c("cds", "utr")
                )
            )
        } else {
            # All other exons <- make one exon
            gene.map <- rbind(
                gene.map,
                data.table(
                    start = exon.starts[i],
                    end = exon.ends[i],
                    ymin = -0.5,
                    ymax = 0.5,
                    annotation = "cds"
                )
            )
        }

        len.exon <- (exon.ends[i] - exon.starts[i]) + 20
        last.fake <- last.pos + len.exon
        pos.map <- rbind(
            pos.map,
            data.table(
                pos = (exon.starts[i] - 10):(exon.ends[i] + 10),
                fake.pos = last.pos:last.fake
            )
        )
        last.pos <- last.fake + 1
    }

    # Then we merge that "fake" map on top of the actual gene map
    gene.map <- merge(gene.map, pos.map, by.x = "start", by.y = "pos")
    setnames(gene.map, "fake.pos", "fake.start")
    gene.map <- merge(gene.map, pos.map, by.x = "end", by.y = "pos")
    setnames(gene.map, "fake.pos", "fake.end")


    # Process variants
    # variants[,log.p:=-log10(pval)]
    variants[, log.q := -log10(qval)]
    variants <- merge(variants, pos.map, by.x = "GENPOS", by.y = "pos", all.x = T)
    variants <- variants[!is.na(fake.pos)]
    variants[, mod.p := if_else(BETA < 0, T, F)]

    if (max(variants[, log.q], na.rm = T) > 10) {
        warning(paste0("A variant in ", gene_info[, SYMBOL], " has a log10 q. value greater than 10..."))
    }

    fake.coding.start <- pos.map[pos == coding.start, fake.pos]
    fake.coding.end <- pos.map[pos == coding.end, fake.pos]
    return(list(gene.map = gene.map, variants = variants, fake.coding.start = fake.coding.start, fake.coding.end = fake.coding.end,pos.map=pos.map))
}