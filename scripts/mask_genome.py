#!/usr/bin/env python3
"""Build a 'masked' GRCh38 FASTA that keeps real sequence only at target loci.

Every chromosome that contains a target locus is written at its full GRCh38
length, but all bases outside the loci are replaced by N. Alignments to this
reference therefore carry true GRCh38 coordinates, so the demo BAM files can
be compared directly with recount3 junctions and GENCODE v26 annotation.

Usage: mask_genome.py LOCI.bed CHROM_SIZES.tsv LOCI_SEQ.fa OUT.fa
  LOCI.bed        chrom, 0-based start, end, name
  CHROM_SIZES.tsv chrom and length (the first two columns of a .fai file)
  LOCI_SEQ.fa     output of `samtools faidx REF chrom:start-end ...`
"""
import sys

loci_bed, sizes_tsv, loci_fa, out_fa = sys.argv[1:5]

sizes = {}
for line in open(sizes_tsv):
    f = line.split("\t")
    sizes[f[0]] = int(f[1])

# Read the region sequences fetched by samtools (headers look like >chr8:100-200).
regions, name, chunks = {}, None, []
for line in open(loci_fa):
    line = line.strip()
    if line.startswith(">"):
        if name:
            regions[name] = "".join(chunks)
        name, chunks = line[1:], []
    else:
        chunks.append(line)
if name:
    regions[name] = "".join(chunks)

by_chrom = {}
for line in open(loci_bed):
    chrom, start, end, _ = line.rstrip("\n").split("\t")[:4]
    by_chrom.setdefault(chrom, []).append((int(start), int(end)))

with open(out_fa, "w") as out:
    for chrom in sorted(by_chrom):
        seq = bytearray(b"N" * sizes[chrom])
        for start, end in by_chrom[chrom]:
            key = f"{chrom}:{start + 1}-{end}"
            s = regions[key].encode()
            seq[start:start + len(s)] = s
        out.write(f">{chrom}\n")
        for i in range(0, len(seq), 60):
            out.write(seq[i:i + 60].decode() + "\n")
print(f"wrote {len(by_chrom)} masked chromosomes to {out_fa}")
