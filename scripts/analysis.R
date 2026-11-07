# DESeq2: snf2 mutant vs wild type (3 vs 3). Run from the repo root.
library(DESeq2)
library(ggplot2)

counts <- read.delim("data/counts_matrix.tsv", row.names = 1, check.names = FALSE)
info <- read.delim("data/sample_info.tsv", row.names = 1)
stopifnot(all(rownames(info) == colnames(counts)))
info$condition <- factor(info$condition, levels = c("wild_type", "mutant"))

dds <- DESeqDataSetFromMatrix(countData = counts, colData = info, design = ~ condition)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "mutant", "wild_type"))
summary(res)

gn <- read.delim("data/gene_names.tsv", header = FALSE, col.names = c("gene_id", "symbol"))
res$symbol <- gn$symbol[match(rownames(res), gn$gene_id)]

vsd <- vst(dds, blind = FALSE)
ggsave("figures/pca.png", plotPCA(vsd, intgroup = "condition"), width = 6, height = 5, dpi = 150)

df <- as.data.frame(res)
df <- df[!is.na(df$padj), ]
df$sig <- ifelse(df$padj < 0.05 & abs(df$log2FoldChange) > 1, "significant", "not significant")
p <- ggplot(df, aes(log2FoldChange, -log10(padj), colour = sig)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_colour_manual(values = c("grey60", "red")) +
  labs(title = "snf2 mutant vs wild type (yeast)", x = "log2 fold change", y = "-log10 adjusted p-value") +
  theme_minimal()
ggsave("figures/volcano.png", p, width = 7, height = 5, dpi = 150)

out <- as.data.frame(res)
out$gene_id <- rownames(out)
write.csv(out[order(out$padj), ], "results/yeast_deseq2_results.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
