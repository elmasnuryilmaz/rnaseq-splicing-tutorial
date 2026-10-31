#!/usr/bin/env python3
"""Subset a GENCODE GTF to a list of genes and write three files.

1. <prefix>.gtf     all GTF lines for the requested genes (for featureCounts,
                    HISAT2 splice sites and rMATS)
2. <prefix>.bed12   one BED12 line per transcript (for RSeQC)
3. <prefix>.loci.bed  one padded region per gene (used to mask the genome)

Usage: subset_gtf.py GENCODE.gtf.gz GENE_LIST.txt PAD_BP OUT_PREFIX
Only the Python standard library is used.
"""
import gzip
import re
import sys
from collections import defaultdict

gtf_path, gene_list_path, pad, prefix = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
wanted = {g.strip() for g in open(gene_list_path) if g.strip() and not g.startswith("#")}

attr_re = re.compile(r'(\S+) "([^"]*)"')
primary_re = re.compile(r'^chr([0-9]+|X|Y|M)$')
genes = {}                      # gene_name -> (chrom, start, end, strand)
tx_exons = defaultdict(list)    # transcript_id -> [(start, end)]
tx_info = {}                    # transcript_id -> (chrom, strand, gene_name)
found = set()

with gzip.open(gtf_path, "rt") as gtf, open(prefix + ".gtf", "w") as out_gtf:
    for line in gtf:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        attrs = dict(attr_re.findall(f[8]))
        name = attrs.get("gene_name")
        if name not in wanted:
            continue
        # Keep primary chromosomes only (skip patches, haplotypes, scaffolds).
        if not primary_re.match(f[0]):
            continue
        out_gtf.write(line)
        found.add(name)
        start, end = int(f[3]), int(f[4])
        if f[2] == "gene":
            genes[name] = (f[0], start, end, f[6])
        elif f[2] == "exon":
            tid = attrs["transcript_id"]
            tx_exons[tid].append((start, end))
            tx_info[tid] = (f[0], f[6], name)

missing = wanted - found
if missing:
    sys.exit("Genes not found in GTF: " + ", ".join(sorted(missing)))

with open(prefix + ".bed12", "w") as bed:
    for tid, exons in tx_exons.items():
        chrom, strand, name = tx_info[tid]
        exons.sort()
        s0 = exons[0][0] - 1                       # BED is 0-based, half-open
        e0 = exons[-1][1]
        sizes = ",".join(str(e - s + 1) for s, e in exons)
        starts = ",".join(str(s - 1 - s0) for s, e in exons)
        bed.write(f"{chrom}\t{s0}\t{e0}\t{tid}|{name}\t0\t{strand}\t{s0}\t{e0}\t0\t"
                  f"{len(exons)}\t{sizes}\t{starts}\n")

with open(prefix + ".loci.bed", "w") as loci:
    for name, (chrom, start, end, strand) in sorted(genes.items(), key=lambda x: (x[1][0], x[1][1])):
        loci.write(f"{chrom}\t{max(0, start - 1 - pad)}\t{end + pad}\t{name}\n")

print(f"{len(genes)} genes, {len(tx_exons)} transcripts written with prefix {prefix}")
