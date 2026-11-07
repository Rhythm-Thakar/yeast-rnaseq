# Yeast bulk RNA-seq: snf2 mutant vs wild type

Raw reads to differential expression, using 3 wild-type and 3 snf2 knock-out samples from the public study PRJEB5348 (ENA).

## Pipeline
FastQC/MultiQC -> HISAT2 -> samtools -> featureCounts -> DESeq2 (R)

Commands: `scripts/pipeline.sh` (WSL) and `scripts/analysis.R` (R, run from the repo root).
Versions: HISAT2 2.2.1, featureCounts 2.0.3, samtools 1.13, MultiQC 1.12, DESeq2 1.52.0, R 4.6.2.

## Design choices
- The study has 48 libraries per condition, each sequenced in 7 lanes. Lanes of one library are technical repeats, so I took one lane-1 run from each of 3 different libraries per condition to get independent replicates. Selection rules were fixed before looking at candidates: lane 1, at least 1M reads, first 3 libraries per condition in sorted order.
- Reads are 51 bp single-end. No trimming: per-base quality and adapter content passed FastQC in all six samples.
- Library is unstranded (`-s 0` assigned about 896k reads on one sample vs about 465k and 436k for `-s 1`/`-s 2`).
- Duplicates not removed (RNA-seq). Multi-mapping reads (18-26% of alignment records per sample) are excluded by featureCounts.
- Genes kept: at least 10 reads in at least 3 samples (5,403 of 6,477).

## Results
- 95.4-97.3% overall alignment rate across the six samples
- At adjusted p < 0.1: 249 genes up and 489 down in the mutant; 284 genes pass adjusted p < 0.05 and |log2 fold change| > 1
- SNF2 (YOR290C) is strongly down in the mutant (log2 fold change -7.66, adjusted p 2.4e-22), as expected for the knock-out
- PC1 (82% of variance) separates the two conditions

Files: `results/yeast_deseq2_results.csv`, `figures/pca.png`, `figures/volcano.png`, `data/`.

## Limits
Only 3 vs 3 samples, so less power than the full 48 vs 48 design. No functional enrichment yet.
