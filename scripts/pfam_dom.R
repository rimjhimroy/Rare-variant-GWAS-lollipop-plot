
make_pfam_dom <- function(gene_model, genename) {
    ensembl <- useMart("ensembl")
    ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
    searchAttributes(ensembl, "mane")


    gene_loc <- getBM(
        attributes = c("chromosome_name", "start_position", "end_position", "strand", "ensembl_transcript_id", "transcript_is_canonical", "refseq_mrna", "transcript_mane_select", "hgnc_symbol", "pfam", "pfam_start", "pfam_end", "namespace_1003"),
        filters = c("hgnc_symbol"),
        values = list(genename),
        mart = ensembl
    )
    gene_loc <- gene_loc[which(gene_loc$transcript_is_canonical == 1), ]
    gene_length <- gene_loc$end_position - gene_loc$start_position
    canonical_trans <- unique(gene_loc$ensembl_transcript_id)
    gene_Start <- unique(gene_loc$start_position)
    gene_End <- unique(gene_loc$end_position)
    refseqid <- gsub("\\..*", "", unique(gene_loc$transcript_mane_select))
    txid <- unique(gene_loc$ensembl_transcript_id)



    gene_annot <- data.table(getBM(
        attributes = c("ensembl_transcript_id", "ensembl_transcript_id_version", "chromosome_name", "start_position", "end_position", "strand", "cdna_coding_start", "cdna_coding_end", "genomic_coding_start", "genomic_coding_end", "exon_chrom_start", "exon_chrom_end", "cds_start", "cds_end"),
        filters = c("hgnc_symbol"),
        values = list(genename),
        mart = ensembl
    ))
    gene_annot <- gene_annot[which(gene_annot$ensembl_transcript_id == canonical_trans), ]

    # refseqids = c("NM_005359","NM_000546")


    edb <- EnsDb.Hsapiens.v86
    pdoms <- proteins(edb,
        filter = ~ tx_id %in% txid &
            protein_domain_source == "pfam",
        columns = c(
            "protein_domain_id", "prot_dom_start",
            "prot_dom_end", "protein_domain_source", "interpro_accession"
        )
    )
    pdoms <- as.data.table(pdoms)

    ipro <- getBM(
        attributes = c("refseq_mrna", "interpro", "interpro_description"),
        filters = "refseq_mrna",
        values = refseqid,
        mart = ensembl
    )
    pdoms <- merge(pdoms, ipro, by.x = "interpro_accession", by.y = "interpro")
    pdoms_rng <- IRanges(
        start = pdoms$prot_dom_start, end = pdoms$prot_dom_end,
        names = pdoms$protein_id # nolint
    )

    pdoms_gnm <- proteinToGenome(pdoms_rng, edb)
    pdoms_gnm_grng <- unlist(GRangesList(pdoms_gnm))
    pdoms_gnm_grng$id <- rep(pdoms$protein_domain_id, lengths(pdoms_gnm))
    pdoms_gnm_grng$grp <- rep(1:nrow(pdoms), lengths(pdoms_gnm))
    pdom_data <- as.data.table(pdoms_gnm_grng)
    pdom_data <- merge(pdom_data, pdoms, by.x = "id", by.y = "protein_domain_id")
    pdom_dd <- data.table(pdom_data %>%
        group_by(id) %>%
        mutate(Start.Value = min(start)) %>%
        mutate(End.Value = max(end)))
    pdom_dd1 <- pdom_dd[, c("id", "seqnames"), with = F]
    pdom_dd1 <- unique(pdom_dd[, c("id", "seqnames", "Start.Value", "End.Value", "protein_id.x", "tx_id.x", "cds_ok", "prot_dom_start", "prot_dom_end", "interpro_accession", "protein_domain_source", "refseq_mrna", "interpro_description"), with = FALSE])
    names(pdom_dd1) <- c("id", "Chrom", "start", "end", "protein_id", "tx_id", "cds_ok", "prot_dom_start", "prot_dom_end", "interpro_accession", "protein_domain_source", "refseq_mrna", "Domains")

    pdom.map <- merge(pdom_dd1, gene_model$pos.map, by.x = "start", by.y = "pos")
    setnames(pdom.map, "fake.pos", "fake.start")
    pdom.map <- merge(pdom.map, gene_model$pos.map, by.x = "end", by.y = "pos")
    setnames(pdom.map, "fake.pos", "fake.end")

    # txs <- getGeneRegionTrackForGviz(edb, filter = ~ genename == genename & tx_id %in% txid)
    # library(Gviz)

    # ## Define the individual tracks:
    # ## - Ideogram
    # ## ideo_track <- IdeogramTrack(genome = "hg38", chromosome = "chr21")
    # ## - Genome axis
    # gaxis_track <- GenomeAxisTrack()
    # ## - Transcripts
    # gene_track <- GeneRegionTrack(txs,
    #     showId = TRUE, just.group = "right",
    #     name = "", geneSymbol = TRUE, size = 0.5
    # )
    # ## - Protein domains
    # pdom_track <- AnnotationTrack(pdoms_gnm_grng,
    #     group = pdoms_gnm_grng$grp,
    #     id = pdoms_gnm_grng$id, groupAnnotation = "id",
    #     just.group = "right", shape = "box",
    #     name = "Protein domains", size = 0.5
    # )

    ## Generate the plot
    # plotTracks(list(gaxis_track, gene_track, pdom_track))
    # dev.off()
    return(pdom.map)
}