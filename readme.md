# Promoter Region Extraction for GO Enrichment and Motif Analysis

## Overview

This workflow constructs a set of 500 bp promoter intervals anchored to the Transcription Start Site (TSS) of human genes. Promoters are extracted in a strand-sensitive manner — meaning the upstream direction is determined by each gene's orientation — using the `bedtools` suite. The resulting intervals are intended for motif enrichment and Gene Ontology (GO) analysis.

---

## Generated Files

| File | Description |
|---|---|
| `genes_tss_final.bed` | Quality-filtered BED file of TSS positions |
| `promoters_500bp.bed` | 500 bp upstream promoter windows, strand-adjusted |

---

## Setting Up the Environment

Set up a clean Mamba environment running Python 3.12:

```bash
mamba create -n go_enrichment python=3.12
mamba activate go_enrichment
```

Then install all required tools from the Bioconda channel:

```bash
mamba install -c bioconda bedtools samtools emboss
```

---

## Required Input Files

### Reference Genome (hg38)

Obtain the human genome FASTA file from UCSC:

```
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/
```

### Gene Annotation

```
human_gene_annotation.tsv.gz
```

---

## Step 1 — Build a Chromosome Size Reference

Generate a FASTA index to allow fast chromosome-level lookups:

```bash
samtools faidx hg38.fa
```

Extract chromosome names and their lengths into a genome file:

```bash
cut -f1,2 hg38.fa.fai > hg38.genome
```

---

## Step 2 — Convert Annotation Data into TSS BED Format

Parse the gzipped annotation file and reformat it as a standard 6-column BED file. Each row in the output represents a single-nucleotide TSS position:

| Column | Value |
|---|---|
| 1 | Chromosome name (prefixed with `chr`) |
| 2 | TSS start position |
| 3 | TSS end position (start + 1) |
| 4 | Unique feature label (`chr@start-end\|gene`) |
| 5 | Score placeholder (`.`) |
| 6 | Strand (`+` or `-`) |

```bash
zcat human_gene_annotation.tsv.gz | \
awk 'BEGIN{OFS="\t"} NR>1{

    # Discard entries with undefined TSS or gene name
    if($8==-1 || $8=="" || $7=="")
        next

    # Build UCSC-style chromosome name; rename MT → M
    chrom = "chr"$5
    if(chrom == "chrMT")
        chrom = "chrM"

    tss = $8 + 0

    # Map numeric strand encoding to BED convention
    strand = "+"
    if($6 == -1)
        strand = "-"

    # Emit the BED record
    print chrom, tss, tss+1, chrom"@"tss"-"(tss+1)"|"$7, ".", strand

}' > genes_tss_clean.bed
```

Discard any records that reference chromosomes not present in the hg38 genome file:

```bash
grep -Fwf <(cut -f1 hg38.genome) genes_tss_clean.bed > genes_tss_final.bed
```

---

## Step 3 — Expand TSS Points to 500 bp Promoter Intervals

Apply `bedtools slop` to widen each single-nucleotide TSS into a 500 bp upstream window. The `-s` flag ensures the expansion respects each gene's strand:

```bash
bedtools slop \
    -i genes_tss_final.bed \
    -g hg38.genome \
    -l 500 \
    -r 0 \
    -s \
    > promoters_500bp.bed
```

### Strand-Aware Logic

Without strand awareness, "upstream" would always mean lower chromosomal coordinates. The `-s` flag corrects for this:

| Strand | Upstream direction | Expansion moves toward |
|---|---|---|
| `+` (sense) | 5′ end | Smaller coordinates |
| `-` (antisense) | 5′ end | Larger coordinates |

**Worked example — reverse-strand gene:**

```
chrM    4400    4901    chrM@4400-4401|MT-TQ    .    -
```

Because the gene sits on the minus strand, the 500 bp window extends rightward (toward higher coordinates), correctly capturing the upstream regulatory region.

---

## Full Pipeline Script

The complete workflow can be executed as a single shell script:

```bash
#!/usr/bin/env bash
set -euo pipefail

mamba activate go_enrichment

samtools faidx hg38.fa
cut -f1,2 hg38.fa.fai > hg38.genome

zcat human_gene_annotation.tsv.gz | \
awk 'BEGIN{OFS="\t"} NR>1{

    if($8==-1 || $8=="" || $7=="")
        next

    chrom = "chr"$5
    if(chrom == "chrMT")
        chrom = "chrM"

    tss    = $8 + 0
    strand = "+"
    if($6 == -1)
        strand = "-"

    print chrom, tss, tss+1, chrom"@"tss"-"(tss+1)"|"$7, ".", strand

}' > genes_tss_clean.bed

grep -Fwf <(cut -f1 hg38.genome) genes_tss_clean.bed > genes_tss_final.bed

bedtools slop \
    -i genes_tss_final.bed \
    -g hg38.genome \
    -l 500 \
    -r 0 \
    -s \
    > promoters_500bp.bed

echo "Pipeline complete. Output: promoters_500bp.bed"
```

---

## Final Output

The file `promoters_500bp.bed` provides a strand-corrected 500 bp promoter window for every gene in the annotation and is directly usable for motif scanning tools and GO enrichment pipelines.
