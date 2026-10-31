# A small sashimi-style plot built from recount3 files only:
#   - coverage: per-sample BigWig files, scaled to "per million read pairs"
#     and averaged within each group;
#   - arcs: junction read counts (mean per sample in each group);
#   - gene model: GENCODE v26 exons of all transcripts, merged.
# It is a teaching tool; for publication figures see ggsashimi or IGV.

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(ggplot2)
  library(patchwork)
})

group_coverage <- function(region, runs, mapped_pairs) {
  m <- vapply(runs, function(r) {
    v <- as.numeric(import(BigWigFile(get_bigwig(r)), which = region,
                           as = "NumericList")[[1]])
    v / (mapped_pairs[r] / 1e6)
  }, numeric(width(region)))
  rowMeans(m)
}

sashimi_plot <- function(chrom, from, to, strand, groups, jxn_counts, jxn_ranges,
                         novel, txdb, mapped_pairs, min_mean_reads = 1,
                         highlight = NULL, title = NULL, group_cols = NULL) {
  region <- GRanges(chrom, IRanges(from, to))
  # Coverage, one track per group
  cov <- do.call(rbind, lapply(names(groups), function(g)
    data.frame(group = g, pos = from:to,
               cov = group_coverage(region, groups[[g]], mapped_pairs))))
  cov$group <- factor(cov$group, levels = names(groups))
  ymax <- max(cov$cov)

  # Junction arcs below each coverage track
  in_reg <- which(as.character(seqnames(jxn_ranges)) == chrom &
                    start(jxn_ranges) >= from & end(jxn_ranges) <= to &
                    as.character(strand(jxn_ranges)) == strand)
  arcs <- do.call(rbind, lapply(names(groups), function(g) {
    mu <- rowMeans(jxn_counts[in_reg, groups[[g]], drop = FALSE])
    keep <- in_reg[mu >= min_mean_reads]
    mu <- mu[mu >= min_mean_reads]
    if (!length(keep)) return(NULL)
    do.call(rbind, lapply(seq_along(keep), function(i) {
      j <- keep[i]
      t <- seq(0, 1, length.out = 60)
      h <- ymax * (0.15 + 0.45 * log10(1 + mu[i]) / log10(1 + max(mu)))
      data.frame(group = g, id = names(jxn_ranges)[j],
                 x = start(jxn_ranges)[j] - 1 + t * (end(jxn_ranges)[j] - start(jxn_ranges)[j] + 2),
                 y = -h * sin(pi * t), mean_reads = mu[i], novel = novel[j],
                 apex = t == 0.5 | abs(t - 0.5) == min(abs(t - 0.5)))
    }))
  }))
  arcs$group <- factor(arcs$group, levels = names(groups))
  labs_df <- arcs[arcs$apex, ]
  labs_df <- labs_df[!duplicated(paste(labs_df$group, labs_df$id)), ]

  fill_cols <- if (is.null(group_cols)) setNames(rep("grey50", length(groups)), names(groups)) else group_cols
  p1 <- ggplot() +
    geom_area(data = cov, aes(pos, cov, fill = group), alpha = 0.8, show.legend = FALSE) +
    geom_path(data = arcs, aes(x, y, group = id, colour = novel), linewidth = 0.5) +
    geom_label(data = labs_df, aes(x, y, label = round(mean_reads), colour = novel),
               size = 2.6, label.size = 0, label.padding = unit(0.1, "lines"),
               show.legend = FALSE) +
    facet_wrap(~ group, ncol = 1) +
    scale_fill_manual(values = fill_cols) +
    scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey35"),
                        labels = c(`TRUE` = "not in GENCODE v26", `FALSE` = "annotated"),
                        name = "Junction") +
    scale_x_continuous(limits = c(from, to), labels = scales::comma,
                       breaks = scales::breaks_pretty(n = 4)) +
    labs(x = NULL, y = "Coverage per million pairs\n(arcs: mean junction reads)", title = title) +
    theme_tutorial() +
    theme(legend.position = "top", axis.text.x = element_blank(), axis.ticks.x = element_blank())

  # Gene model: merged exons of all GENCODE v26 transcripts on this strand
  ex <- exons(txdb, filter = list(tx_strand = strand))
  ex <- reduce(subsetByOverlaps(ex, region))
  gm <- data.frame(start = pmax(start(ex), from), end = pmin(end(ex), to))
  p2 <- ggplot() +
    geom_segment(aes(x = from, xend = to, y = 0, yend = 0), colour = "grey60") +
    geom_rect(data = gm, aes(xmin = start, xmax = end, ymin = -0.4, ymax = 0.4), fill = "grey30")
  if (!is.null(highlight)) {
    p2 <- p2 + geom_rect(data = highlight, aes(xmin = start, xmax = end, ymin = -0.4, ymax = 0.4),
                         fill = "#D55E00")
  }
  p2 <- p2 +
    scale_x_continuous(limits = c(from, to), labels = scales::comma,
                       breaks = scales::breaks_pretty(n = 4)) +
    scale_y_continuous(limits = c(-0.6, 0.6)) +
    labs(x = paste0(chrom, " (", strand, " strand)"), y = NULL,
         caption = paste0("Grey boxes: GENCODE v26 exons (all transcripts merged).",
                          if (!is.null(highlight)) " Orange box: cryptic exon inferred from junctions." else "")) +
    theme_tutorial() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid = element_blank())
  p1 / p2 + plot_layout(heights = c(5, 1))
}
