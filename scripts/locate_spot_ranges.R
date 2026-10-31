# Estimate where the reads of each panel gene sit inside an SRA run.
#
# Why this works: the submitters uploaded coordinate-sorted BAM files, and SRA
# keeps the reads of each run in that order (chr1, chr2, ... with the
# alignments stripped). So the reads of one gene form one contiguous block of
# spot numbers, and the position of that block is roughly the number of
# fragments that map before the gene. recount3 gives us that number: the
# per-sample coverage BigWig files sum to (fragments x mapped length).
#
# Usage: Rscript scripts/locate_spot_ranges.R LOCI.bed RUN [RUN ...] > ranges.tsv
# Output columns: run, locus, chrom, start, end, first_spot, last_spot

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
})
source("R/helpers.R")

args <- commandArgs(trailingOnly = TRUE)
loci <- read.table(args[1], col.names = c("chrom", "start0", "end", "locus"))
runs <- args[-1]
samples <- read.csv("results/sample_sheet.csv")

# Assumed chromosome order of the submitted BAM files (checked in the chapter).
chroms <- paste0("chr", c(1:22, "X", "Y", "M"))

auc_before <- function(bw, chrom, pos, chrom_auc) {
  inside <- if (pos > 1) {
    summary(bw, GRanges(chrom, IRanges(1, pos)), size = 1L, type = "mean")[[1]]$score * pos
  } else 0
  sum(chrom_auc[seq_len(match(chrom, chroms) - 1)]) + inside
}

out <- list()
for (run in runs) {
  bw <- BigWigFile(get_bigwig(run))
  len <- seqlengths(seqinfo(bw))[chroms]
  chrom_auc <- vapply(chroms, function(ch)
    summary(bw, GRanges(ch, IRanges(1, len[ch])), size = 1L, type = "mean")[[1]]$score * len[ch],
    numeric(1))
  total <- sum(chrom_auc)
  mapped <- samples$mapped_pairs[samples$run == run]
  for (i in seq_len(nrow(loci))) {
    l <- loci[i, ]
    first <- auc_before(bw, l$chrom, l$start0 + 1, chrom_auc) / total * mapped
    last  <- auc_before(bw, l$chrom, l$end, chrom_auc) / total * mapped
    out[[length(out) + 1]] <- data.frame(run = run, locus = l$locus, chrom = l$chrom,
                                         start = l$start0 + 1, end = l$end,
                                         first_spot = round(first), last_spot = round(last))
  }
}
write.table(do.call(rbind, out), stdout(), sep = "\t", quote = FALSE, row.names = FALSE)
