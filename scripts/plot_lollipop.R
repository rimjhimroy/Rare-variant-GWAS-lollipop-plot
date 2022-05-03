plot_lollipop <- function(gene_model, pfam_dom, qvalue) {
    # Default Theme:
    theme <- theme(
        panel.background = element_rect(fill = "white"),
        line = element_line(size = 1, colour = "black", lineend = "round"),
        axis.line = element_line(size = 1),
        text = element_text(size = 16, face = "bold", colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(size = 1, colour = "black"),
        axis.ticks.length = unit(.1, "cm"),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.text.x = element_blank(),
        legend.position = "right",
        # axis.line.x = element_blank(), # remove x axis line
        panel.grid.major = element_blank(),
        legend.key = element_blank()
     
    )

    max(gene_model$variants[, log.q])
    ylimits <- as.integer((max(gene_model$variants[, log.q]) / 2)) * 2 + 2
    ylabs <- c(rev(seq(0, ylimits, by = 2)), seq(0, ylimits, by = 2))
    ybreaks <- c(seq(-ylimits, 0, by = 2) - 1.5, seq(0, ylimits, by = 2) + 1.5)
    ylim_min <- -13
    ylim_max <- 30

    pfam_dom$ymin <- ylim_min
    pfam_dom$ymax <- ylim_min + 1
    x_axis_name <- symbol
    gene.plot <- ggplot() +
        geom_segment(aes(x = gene_model$fake.coding.start, xend = gene_model$fake.coding.end, y = 0, yend = 0), size = 1) +
        geom_hline(yintercept = c(1.5 + -log10(qvalue), -1.5 - -log10(qvalue)), colour = "red", linetype = 2) +
        geom_segment(data = gene_model$variants, aes(x = fake.pos, xend = fake.pos, y = 0, yend = if_else(mod.p == T, -1.5 - log.q, 1.5 + log.q))) +
        geom_segment(aes(x = gene_model$fake.coding.start - 60, xend = gene_model$fake.coding.start - 60, y = 2.5, yend = 8.5), arrow = arrow(length = unit(0.02, "npc")), size = 1, lineend = "round", linejoin = "round", colour = "darkgrey") +
        geom_segment(aes(x = gene_model$fake.coding.start - 60, xend = gene_model$fake.coding.start - 60, y = -2.5, yend = -8.5), arrow = arrow(length = unit(0.02, "npc")), size = 1, lineend = "round", linejoin = "round", colour = "darkgrey") +
        annotate(geom = "text", x = gene_model$fake.coding.start - 40, y = 5.5, hjust = 0, label = "Pos. β", size = 5) +
        annotate(geom = "text", x = gene_model$fake.coding.start - 40, y = -5.5, hjust = 0, label = "Neg. β", size = 5) +
        geom_rect(data = gene_model$gene.map[annotation != "utr"], aes(xmin = fake.start, xmax = fake.end, ymin = ymin, ymax = ymax), fill = "lightblue", colour = "black") +
        geom_point(data = gene_model$variants, aes(fake.pos, if_else(mod.p == T, -1.5 - log.q, 1.5 + log.q), size = A1FREQ)) +
        scale_x_continuous(name = x_axis_name, limits = c(gene_model$fake.coding.start - 75, gene_model$fake.coding.end + 50), breaks = c(gene_model$fake.coding.start, gene_model$fake.coding.end)) +
        scale_y_continuous(name = expression(bold(-log[10](italic(q)))), breaks = ybreaks, labels = ylabs, limits = c(ylim_min, ylim_max)) + # change limits to adjust labels
        scale_size_continuous(guide = guide_legend(title = "Allele FREQ")) +
        coord_flex_cart(bottom = capped_horizontal(capped = "both"), left = capped_vertical(capped = "both")) + # Comes from the "lemon" package
        theme 


    gene.plot1 <- gene.plot + geom_segment(data = pfam_dom,aes(x = gene_model$fake.coding.start, xend = gene_model$fake.coding.end, y =  (ymin+ ymax)/2, yend = (ymin+ ymax)/2), size = 1) +
        geom_rect(data = pfam_dom, aes(xmin = fake.start, xmax = fake.end, ymin = ymin, ymax = ymax, fill = Domains), colour = "black")
    
    if (nrow(gene_model$variants[log.q > -log10(qvalue)]) > 0) {
        gene.plot2 <- gene.plot1 + geom_text_repel(
            data = gene_model$variants[log.q > -log10(qvalue)], aes(fake.pos, if_else(mod.p == T, -1.5 - log.q, 1.5 + log.q), label = paste0("ID:", ID, "\nA1FQ: ", A1FREQ, "\nTrait: ", gsub("_", " ", trait), "\nCSQ: ", consequence)),
            nudge_x = 1,
            box.padding = 0.5,
            nudge_y = 4,
            segment.linetype = 6,
            segment.curvature = -0.1,
            segment.ncp = 3,
            point.padding = 0,
            segment.angle = 20,
            fontface = "bold",
            segment.color = "brown"
        ) 
    }
    return(gene.plot2)
}