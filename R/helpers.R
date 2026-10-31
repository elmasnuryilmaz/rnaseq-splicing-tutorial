# Shared helpers for the tutorial chapters.
#
# Every chapter sources this file, so each chapter can be rendered on its own.
# Downloads are cached in `data-cache/` (ignored by git). Small result tables
# that later chapters need are written to `results/` (tracked by git).

suppressPackageStartupMessages({
  library(SummarizedExperiment)
})

# ---- Fixed project settings ------------------------------------------------
# One study, one reference, one annotation for the whole book.
PROJECT_ID  <- "SRP166282"        # SRA study (GEO series GSE121569)
GENOME      <- "GRCh38"
ANNOTATION  <- "gencode_v26"      # recount3 code G026 = GENCODE release 26
GTF_NAME    <- "gencode.v26.chr_patch_hapl_scaff.annotation.gtf.gz"

# recount3 files are mirrored on the AWS Open Data registry. We point recount3
# at that mirror; the default host (duffel.rail.bio) serves identical files.
R3_URL <- "https://recount-opendata.s3.amazonaws.com/recount3/release"
options(recount3_url = R3_URL)

CACHE_DIR <- Sys.getenv("TUTORIAL_CACHE", unset = "data-cache")

# Default thresholds (fixed before looking at any result; never tuned).
DE_PADJ      <- 0.05
DE_LFC       <- 1
SPLICE_FDR   <- 0.05
SPLICE_DPSI  <- 0.1

cache_path <- function(...) {
  dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
  file.path(CACHE_DIR, ...)
}

# ---- recount3 downloads ----------------------------------------------------
recount3_project <- function() {
  proj <- recount3::available_projects(organism = "human")
  proj[proj$project == PROJECT_ID & proj$project_type == "data_sources", ]
}

load_gene_rse <- function() {
  f <- cache_path(paste0("rse_gene_", PROJECT_ID, "_G026.rds"))
  if (!file.exists(f)) {
    rse <- recount3::create_rse(recount3_project(), type = "gene",
                                annotation = ANNOTATION)
    saveRDS(rse, f)
  }
  readRDS(f)
}

load_jxn_rse <- function() {
  f <- cache_path(paste0("rse_jxn_", PROJECT_ID, "_ALL.rds"))
  if (!file.exists(f)) {
    rse <- recount3::create_rse(recount3_project(), type = "jxn",
                                jxn_format = "ALL")
    saveRDS(rse, f)
  }
  readRDS(f)
}

# Gene-level read counts (coverage sums divided by the mapped read length,
# see chapter 2), with the tidy sample sheet stored in colData.
load_counts_rse <- function() {
  f <- cache_path("rse_gene_counts.rds")
  if (!file.exists(f)) {
    rse <- load_gene_rse()
    assay(rse, "counts") <- recount3::compute_read_counts(rse)
    saveRDS(rse, f)
  }
  rse <- readRDS(f)
  assays(rse) <- assays(rse)["counts"]   # DESeq2 expects "counts" first
  st <- sample_table(rse)
  rse <- rse[, st$run]
  colData(rse) <- cbind(colData(rse), st[, c("gsm", "title", "sirna", "replicate", "batch")])
  rse
}

# GENCODE v26 GTF, the same file recount3 used for its gene counts.
gencode_gtf <- function() {
  f <- cache_path(GTF_NAME)
  if (!file.exists(f)) {
    download.file(file.path(R3_URL, "human/annotations", GTF_NAME), f,
                  mode = "wb", quiet = TRUE)
  }
  f
}

gencode_txdb <- function() {
  f <- cache_path("gencode_v26.sqlite")
  if (!file.exists(f)) {
    txdb <- txdbmaker::makeTxDbFromGFF(gencode_gtf(), format = "gtf",
                                       dataSource = "GENCODE v26",
                                       organism = "Homo sapiens")
    AnnotationDbi::saveDb(txdb, f)
  }
  AnnotationDbi::loadDb(f)
}

# Per-sample coverage BigWig files produced by recount3 (one per run).
bigwig_url <- function(run) {
  sprintf("%s/human/data_sources/sra/base_sums/%s/%s/%s/sra.base_sums.%s_%s.ALL.bw",
          R3_URL, substr(PROJECT_ID, nchar(PROJECT_ID) - 1, nchar(PROJECT_ID)),
          PROJECT_ID, substr(run, nchar(run) - 1, nchar(run)), PROJECT_ID, run)
}

get_bigwig <- function(run) {
  dir.create(cache_path("bigwig"), showWarnings = FALSE)
  f <- cache_path("bigwig", basename(bigwig_url(run)))
  if (!file.exists(f)) download.file(bigwig_url(run), f, mode = "wb", quiet = TRUE)
  f
}

# ---- Sample sheet ----------------------------------------------------------
# SRA stores sample attributes as "key;;value|key;;value". Pull one key out.
parse_attr <- function(x, key) {
  sub(paste0(".*", key, ";;([^|]*).*"), "\\1", x)
}

sample_table <- function(rse) {
  cd <- colData(rse)
  attr <- cd$sra.sample_attributes
  sirna <- parse_attr(attr, "sirna type")
  date  <- parse_attr(attr, "date")
  title <- sub("^GSM[0-9]+: ([^;]*);.*$", "\\1", cd$sra.experiment_title)
  st <- data.frame(
    run       = cd$external_id,
    gsm       = sub(":.*$", "", cd$sra.experiment_title),
    title     = title,
    sirna     = factor(sirna, levels = c("none", "fluorescent", "scramble", "TDP43")),
    replicate = as.integer(sub(".*rep", "", title)),
    date      = date,
    batch     = factor(ifelse(date == "2013-01-10", "B1", "B2")),
    layout    = cd$sra.library_layout,
    read_pairs   = as.numeric(cd$sra.num_spots),
    mapped_pairs = as.numeric(cd$`recount_qc.star.all_mapped_reads`),
    unique_pct   = as.numeric(cd$`recount_qc.star.uniquely_mapped_reads_%_both`),
    row.names = cd$external_id
  )
  st[order(st$sirna, st$batch, st$replicate), ]
}

# Shared colours (colour-blind safe, Okabe-Ito).
SIRNA_COLS <- c(none = "#999999", fluorescent = "#E69F00",
                scramble = "#0072B2", TDP43 = "#D55E00")

theme_tutorial <- function(base_size = 12) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(fill = "grey95"))
}
