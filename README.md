# Yeast bulk RNA-seq: snf2 mutant vs wild type

Raw reads to differential expression, using 3 wild-type and 3 snf2 knock-out samples from the public study PRJEB5348 (ENA).

## Pipeline
FastQC/MultiQC -> HISAT2 -> samtools -> featureCounts -> DESeq2 (R) -> clusterProfiler (GO enrichment)

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
  Raw counts: 148, 220, 169 in wild type vs 2, 0, 2 in the mutant
- PC1 (82% of variance) separates the two conditions

Files: `results/yeast_deseq2_results.csv`, `figures/pca.png`, `figures/volcano.png`, `data/`.

## Limits
Only 3 vs 3 samples, so less power than the full 48 vs 48 design.

## Functional enrichment (GO Biological Process)
`scripts/enrichment.R` (clusterProfiler 4.20.3, org.Sc.sgd.db). Test set: the 284 strongly changed genes (adjusted p < 0.05, |log2FC| > 1; 211 down, 73 up). Background: all 5,403 tested genes, not the whole genome.
- 55 terms enriched at adjusted p < 0.05, 19 after merging redundant terms
- All of the signal comes from the down genes (73 terms); the 73 up genes give none, which may reflect low power for a short list
- Broad terms: oxoacid metabolic process (42/282 genes vs 368/5,397 in the background, adjusted p 4.4e-4) and transmembrane transport (44/282 vs 446/5,397)
- Sharper terms: fructose transmembrane transport (6 genes, all down: HXT2, HXT3, HXT4, HXT6, HXT7, HXK1) and IMP biosynthesis (6 genes)
- This is an association, not proof that a pathway changes, and with 3 vs 3 samples these are hypotheses to test

Files: `results/go_bp_enrichment*.csv`, `figures/go_bp_dotplot.png`.
