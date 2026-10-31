# From counts to cryptic exons

A hands-on tutorial on **differential expression** and **alternative splicing**
analysis with public RNA-seq data, written as a Quarto book.

The running example is TDP-43 knockdown in human motor neurons
(Klim et al., *Nature Neuroscience* 2019,
[doi:10.1038/s41593-018-0300-4](https://doi.org/10.1038/s41593-018-0300-4);
SRA **SRP166282**, GEO **GSE121569**). Loss of TDP-43 lets "cryptic" exons into
mRNAs such as *STMN2*, which makes it a good teaching case: the expected
result is known, and the most interesting events are, by definition, missing
from the annotation. Gene and junction counts come from
[recount3](https://rna.recount.bio/), so the main analysis needs no alignment.

## What the tutorial covers

| Chapter | Topic |
|---|---|
| 1 | Setup and reproducibility: renv, Quarto freeze, session info |
| 2 | Experimental design and data inventory: groups, replicate unit, batches, library type, annotation version |
| 3 | Quality control: library sizes, knockdown check, PCA, sample correlation, batch, sex check |
| 4 | Strandedness: measuring it with RSeQC and featureCounts `-s 0/1/2`, and keeping featureCounts, HISAT2 and rMATS settings consistent |
| 5 | Differential expression with DESeq2 (apeglm shrinkage) and edgeR quasi-likelihood; svaseq sensitivity analysis |
| 6 | MA plot, volcano plot, heatmap and Hallmark GSEA |
| 7 | Cryptic exons from junction counts, with sashimi-style plots |
| 8 | Event-level splicing with rMATS-turbo: the command, the output, and filtering by column name |
| 9 | Reporting: Methods text, wording, "significant" versus "validated" |

Every chapter ends with a short *Check your understanding* box.

### Analysis standards taught throughout

- One genome (GRCh38) and one annotation (GENCODE v26) for everything.
- DESeq2 with an explicit reference level (`relevel`), a low-count filter
  (at least 10 counts in at least as many samples as the smallest group), and
  apeglm-shrunken fold changes reported next to the raw ones.
- Fixed thresholds: padj < 0.05 and |log2FC| > 1 for DE; FDR < 0.05 and
  |ΔPSI| > 0.1 for splicing.
- Volcano plots on padj, heatmaps on VST values with row z-scores, GSEA with a
  stated ranking metric.
- rMATS with the real read length, measured strandedness, `--novelSS` only for
  cryptic events, JC and JCEC kept separate, filters by column name.
- Cryptic exon claims need junction counts and a visual check; RNA-seq alone
  does not show protein or functional change.

## Building the book locally

Requirements: R 4.4, Quarto 1.4 or newer, and about 1.5 GB of free disk space
for downloads.

```bash
git clone https://github.com/elmasnuryilmaz/rnaseq-splicing-tutorial.git
cd rnaseq-splicing-tutorial
```

Restore the R packages (exact versions from `renv.lock`, Bioconductor 3.20):

```r
install.packages("renv")
renv::restore()
```

Render:

```bash
quarto render
```

Because the project uses `freeze: auto`, `quarto render` reuses the stored
results in `_freeze/` for chapters whose source has not changed. To re-run a
chapter, render it on its own, for example
`quarto render 05-differential-expression.qmd`; data are downloaded on first
use into `data-cache/`.

An alternative to renv is the conda environment in `envs/r-conda.yml`, which is
the environment the book was actually rendered in.

### Rebuilding the demo data (optional)

Chapters 4 and 8 use a small set of aligned reads (6 samples, 17 genes) whose
outputs are stored in `data/demo/`. To rebuild them you need the command-line
tools in `envs/ngs.yml`:

```bash
mamba env create -f envs/ngs.yml
mamba activate rnaseq-tutorial-ngs
RSCRIPT=$(which Rscript) bash scripts/make_demo_data.sh
```

The script streams only the reads of the panel genes from the public SRA copy
on AWS; it never downloads full FASTQ files or builds a full genome index.

## Publishing

`.github/workflows/publish.yml` renders the book from the frozen results and
deploys it to GitHub Pages on every push to `main`. It needs neither R nor
data. GitHub Pages must be enabled in the repository settings with
**Source: GitHub Actions**. Pages for private repositories requires a paid
GitHub plan; on a free plan the repository has to be public.

## Repository layout

```
_quarto.yml        book configuration
*.qmd              chapters
R/                 shared helpers and the sashimi plot function
scripts/           demo data pipeline (bash, Python, R)
data/demo/         stored outputs of the demo pipeline
results/           result tables written by the chapters
_freeze/           stored chapter outputs used by CI
renv.lock, envs/   software versions
```

## Citation and licence

When you use the analyses, please cite the original study (Klim et al. 2019)
and the methods papers listed in the book's reference section.

The text and figures are licensed under [CC BY 4.0](LICENSE-TEXT.md) and the
code under the [MIT licence](LICENSE).
