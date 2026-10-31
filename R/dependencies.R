# This file is never run. It lists packages that the book uses indirectly
# (for example apeglm through DESeq2::lfcShrink(type = "apeglm")) so that
# renv::snapshot() records them in renv.lock.
library(apeglm)
library(txdbmaker)
library(AnnotationDbi)
library(BiocManager)
library(knitr)
library(rmarkdown)
library(msigdbr)
library(patchwork)
library(ggrepel)
library(scales)
library(tibble)
