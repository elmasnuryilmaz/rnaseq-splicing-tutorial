#!/usr/bin/env bash
# Build the small demo data set used in the strandedness and rMATS chapters.
#
# What it does (and why it stays cheap):
#   1. Subsets the GENCODE v26 GTF to a short gene panel.
#   2. Builds a HISAT2 index of GRCh38 in which everything outside that panel
#      is masked to N (real coordinates, tiny index).
#   3. Estimates where each gene's reads sit inside each SRA run. The
#      submitters uploaded coordinate-sorted BAM files, so SRA stores every
#      run in genomic order and a gene's reads form one block of spot numbers
#      (see scripts/locate_spot_ranges.R).
#   4. Streams only those blocks from the AWS Open Data copy of SRA, aligns
#      them and keeps the alignments inside the panel (no full FASTQ download;
#      see scripts/fetch_locus_reads.py).
#   5. Measures strandedness (RSeQC and featureCounts -s 0/1/2).
#   6. Runs rMATS-turbo, knockdown (b1) versus scramble (b2).
#
# Tools: run inside the conda environment described in envs/ngs.yml. Step 3
# needs R with rtracklayer (the renv library); point RSCRIPT at that Rscript.
# Runtime on 4 cores: roughly 30 to 60 minutes, dominated by read streaming.
set -euo pipefail

THREADS=${THREADS:-4}
RSCRIPT=${RSCRIPT:-Rscript}
PAD=${PAD:-10000}                 # bp added on both sides of each gene
READ_LEN=76                       # 2 x 76 bp (SRA spot length 152)
WORK=${WORK:-work/demo}
OUT=${OUT:-data/demo}

GTF_URL=https://recount-opendata.s3.amazonaws.com/recount3/release/human/annotations/gencode.v26.chr_patch_hapl_scaff.annotation.gtf.gz
# GRCh38 sequence (1000 Genomes analysis set; read by region, never downloaded whole)
FASTA_URL=https://1000genomes.s3.amazonaws.com/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
SRA_URL=https://sra-pub-run-odp.s3.amazonaws.com/sra

# b1 = TDP-43 knockdown, b2 = scramble siRNA. All six are from batch B2.
KD_RUNS=(SRR8083874 SRR8083875 SRR8083880)     # TDP43 siRNA rep3, rep4, rep5
CT_RUNS=(SRR8083872 SRR8083873 SRR8083878)     # scramble siRNA rep3, rep4, rep5

mkdir -p "$WORK"/{ref,bam} "$OUT"/{strandedness,rmats,align}

# ---- 1. Gene panel and annotation subset ------------------------------------
[ -s "$WORK/ref/gencode.v26.gtf.gz" ] || curl -sS -o "$WORK/ref/gencode.v26.gtf.gz" "$GTF_URL"
python3 scripts/subset_gtf.py "$WORK/ref/gencode.v26.gtf.gz" scripts/demo_genes.txt \
  "$PAD" "$WORK/ref/panel"
cp "$WORK/ref/panel.loci.bed" "$OUT/panel_loci.bed"

# ---- 2. Masked GRCh38 and HISAT2 index --------------------------------------
if [ ! -s "$WORK/ref/panel_masked.fa" ]; then
  curl -sS -o "$WORK/ref/grch38.fa.fai" "$FASTA_URL.fai"
  regions=$(awk '{printf "%s:%d-%d ", $1, $2 + 1, $3}' "$WORK/ref/panel.loci.bed")
  # shellcheck disable=SC2086
  samtools faidx "$FASTA_URL" $regions > "$WORK/ref/panel_loci.fa"
  rm -f ./*.fai   # htslib may leave a copy of the remote index here
  python3 scripts/mask_genome.py "$WORK/ref/panel.loci.bed" "$WORK/ref/grch38.fa.fai" \
    "$WORK/ref/panel_loci.fa" "$WORK/ref/panel_masked.fa"
fi
if [ ! -s "$WORK/ref/idx.1.ht2" ]; then
  hisat2_extract_splice_sites.py "$WORK/ref/panel.gtf" > "$WORK/ref/panel.ss"
  hisat2_extract_exons.py "$WORK/ref/panel.gtf" > "$WORK/ref/panel.exon"
  hisat2-build -p "$THREADS" --ss "$WORK/ref/panel.ss" --exon "$WORK/ref/panel.exon" \
    "$WORK/ref/panel_masked.fa" "$WORK/ref/idx" > "$WORK/ref/hisat2-build.log" 2>&1
fi

# ---- 3. Where are the panel genes inside each run? --------------------------
"$RSCRIPT" scripts/locate_spot_ranges.R "$WORK/ref/panel.loci.bed" \
  "${KD_RUNS[@]}" "${CT_RUNS[@]}" > "$WORK/spot_ranges.tsv"

# ---- 4. Stream those blocks, align, filter ----------------------------------
rm -f "$OUT/align/spot_blocks.tsv"
for run in "${KD_RUNS[@]}" "${CT_RUNS[@]}"; do
  python3 scripts/fetch_locus_reads.py "$WORK/spot_ranges.tsv" "$run" "$WORK/ref/idx" \
    "$WORK/bam/$run.bam" "$OUT/align/spot_blocks.tsv" "$THREADS"
done

# ---- 5. Strandedness --------------------------------------------------------
for run in "${KD_RUNS[@]}" "${CT_RUNS[@]}"; do
  infer_experiment.py -r "$WORK/ref/panel.bed12" -i "$WORK/bam/$run.bam" -s 200000 \
    > "$OUT/strandedness/$run.infer_experiment.txt"
done
bams=()
for run in "${KD_RUNS[@]}" "${CT_RUNS[@]}"; do bams+=("$WORK/bam/$run.bam"); done
for s in 0 1 2; do
  featureCounts -T "$THREADS" -p --countReadPairs -t exon -g gene_id -s "$s" \
    -a "$WORK/ref/panel.gtf" -o "$WORK/featureCounts_s$s.txt" "${bams[@]}" \
    > /dev/null 2>&1
  cp "$WORK/featureCounts_s$s.txt.summary" "$OUT/strandedness/featureCounts_s$s.summary"
done

# ---- 6. rMATS-turbo ---------------------------------------------------------
kd_bams=(); ct_bams=()
for run in "${KD_RUNS[@]}"; do kd_bams+=("$(realpath "$WORK/bam/$run.bam")"); done
for run in "${CT_RUNS[@]}"; do ct_bams+=("$(realpath "$WORK/bam/$run.bam")"); done
(IFS=,; echo "${kd_bams[*]}") > "$WORK/b1.txt"     # b1: TDP-43 knockdown
(IFS=,; echo "${ct_bams[*]}") > "$WORK/b2.txt"     # b2: scramble control
rm -rf "$WORK/rmats_tmp"
rmats.py --b1 "$WORK/b1.txt" --b2 "$WORK/b2.txt" --gtf "$WORK/ref/panel.gtf" \
  -t paired --readLength "$READ_LEN" --libType fr-unstranded --novelSS \
  --nthread "$THREADS" --od "$WORK/rmats" --tmp "$WORK/rmats_tmp" \
  > "$WORK/rmats.log" 2>&1
cp "$WORK"/rmats/*.MATS.JC.txt "$WORK"/rmats/*.MATS.JCEC.txt "$WORK"/rmats/summary.txt "$OUT/rmats/"
# Keep the group files, with paths made relative to the project root.
for b in b1 b2; do sed "s#$(pwd)/##g" "$WORK/$b.txt" > "$OUT/rmats/$b.txt"; done

# Small STMN2 slices for viewing in IGV (one sample per group).
mkdir -p "$OUT/igv"
stmn2=$(awk '$4 == "STMN2" {printf "%s:%d-%d", $1, $2 + 1, $3}' "$WORK/ref/panel.loci.bed")
samtools view -b -o "$OUT/igv/STMN2_KD_${KD_RUNS[0]}.bam" "$WORK/bam/${KD_RUNS[0]}.bam" "$stmn2"
samtools view -b -o "$OUT/igv/STMN2_scramble_${CT_RUNS[0]}.bam" "$WORK/bam/${CT_RUNS[0]}.bam" "$stmn2"
for b in "$OUT"/igv/*.bam; do samtools index "$b"; done

# Record tool versions next to the outputs.
{
  echo "date: $(date -u +%Y-%m-%d)"
  hisat2 --version | head -1
  samtools --version | head -1
  featureCounts -v 2>&1 | grep -i featurecounts
  echo "RSeQC infer_experiment.py $(infer_experiment.py --version 2>&1 | tail -1)"
  echo "rMATS $(rmats.py --version)"
  echo "sra-tools $(fastq-dump --version | grep -o '[0-9.]*$')"
} > "$OUT/tool_versions.txt"
echo "Demo data written to $OUT"
