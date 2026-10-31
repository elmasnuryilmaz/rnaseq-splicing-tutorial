#!/usr/bin/env python3
"""Stream only the reads of selected genes from a coordinate-ordered SRA run.

For every locus, a window of spots around the estimated block (from
locate_spot_ranges.R) is streamed with fastq-dump, aligned with HISAT2 to the
masked GRCh38 index, and the alignments inside the locus are kept. If the
locus reads reach the edge of the window, the window is widened and the step
repeated, so a poor estimate costs time but never truncates a gene.

Usage: fetch_locus_reads.py RANGES.tsv RUN INDEX OUT.bam LOG.tsv [THREADS]
"""
import csv
import os
import subprocess
import sys
import tempfile

ranges_tsv, run, index, out_bam, log_tsv = sys.argv[1:6]
threads = sys.argv[6] if len(sys.argv) > 6 else "4"
SRA_URL = f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{run}/{run}"
MARGIN = 100_000   # spots added on each side of the estimated block
GUARD = 5_000      # locus reads closer than this to a window edge -> widen

loci = [r for r in csv.DictReader(open(ranges_tsv), delimiter="\t") if r["run"] == run]
tmp = tempfile.mkdtemp(prefix=f"{run}_", dir=os.path.dirname(out_bam) or ".")
locus_bams, log_rows = [], []


def run_cmd(cmd):
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


for locus in loci:
    region = f"{locus['chrom']}:{locus['start']}-{locus['end']}"
    lo = max(1, int(locus["first_spot"]) - MARGIN)
    hi = int(locus["last_spot"]) + MARGIN
    for attempt in range(8):
        prefix = os.path.join(tmp, locus["locus"])
        run_cmd(f"fastq-dump -N {lo} -X {hi} --split-files -O {tmp} {SRA_URL} > /dev/null 2>&1")
        # Unstranded library: no --rna-strandness. --no-softclip keeps reads at 76 bp.
        run_cmd(f"hisat2 -p {threads} --no-softclip -x {index} "
                f"-1 {tmp}/{run}_1.fastq -2 {tmp}/{run}_2.fastq 2> /dev/null | "
                f"samtools view -u -f 2 -q 60 -e '[NM] <= 3' - | "
                f"samtools sort -o {prefix}.all.bam - && samtools index {prefix}.all.bam")
        names = subprocess.run(f"samtools view {prefix}.all.bam {region} | cut -f1",
                               shell=True, capture_output=True, text=True, check=True).stdout.split()
        spots = sorted(int(n.split(".")[-1]) for n in names)
        if not spots:
            first = last = None
            widen_lo = widen_hi = True
        else:
            first, last = spots[0], spots[-1]
            widen_lo = first - lo < GUARD and lo > 1
            widen_hi = hi - last < GUARD
        if not (widen_lo or widen_hi):
            break
        lo = max(1, lo - MARGIN) if widen_lo else lo
        hi = hi + MARGIN if widen_hi else hi
    run_cmd(f"samtools view -b -o {prefix}.bam {prefix}.all.bam {region}")
    locus_bams.append(f"{prefix}.bam")
    log_rows.append([run, locus["locus"], locus["first_spot"], locus["last_spot"],
                     first, last, len(spots), lo, hi, attempt + 1])
    for f in (f"{tmp}/{run}_1.fastq", f"{tmp}/{run}_2.fastq", f"{prefix}.all.bam", f"{prefix}.all.bam.bai"):
        os.remove(f)

run_cmd(f"samtools merge -f -o {tmp}/merged.bam {' '.join(locus_bams)} && "
        f"samtools sort -o {out_bam} {tmp}/merged.bam && samtools index {out_bam}")
run_cmd(f"rm -rf {tmp}")

write_header = not os.path.exists(log_tsv)
with open(log_tsv, "a") as fh:
    w = csv.writer(fh, delimiter="\t")
    if write_header:
        w.writerow(["run", "locus", "est_first_spot", "est_last_spot", "obs_first_spot",
                    "obs_last_spot", "reads_in_locus", "window_lo", "window_hi", "attempts"])
    w.writerows(log_rows)
print(f"{run}: {len(loci)} loci written to {out_bam}")
