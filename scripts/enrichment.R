# GO (Biological Process) enrichment of strongly changed genes. Run from the repo root.
library(clusterProfiler)
library(org.Sc.sgd.db)
library(ggplot2)

res_df <- read.csv("results/yeast_deseq2_results.csv")
gn <- read.delim("data/gene_names.tsv", header = FALSE, col.names = c("gene_id", "symbol"))
res_df$symbol <- gn$symbol[match(res_df$gene_id, gn$gene_id)]

universe <- res_df$gene_id  # background = all 5,403 tested genes
sig <- res_df$gene_id[!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1]
stopifnot(length(universe) == 5403, length(sig) == 284)

run_go <- function(g) enrichGO(g, universe = universe, OrgDb = org.Sc.sgd.db,
                               keyType = "ORF", ont = "BP", pAdjustMethod = "BH",
                               pvalueCutoff = 0.05, qvalueCutoff = 0.2)

ego <- run_go(sig)
ego2 <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
up <- sig[res_df$log2FoldChange[match(sig, res_df$gene_id)] > 0]
dn <- setdiff(sig, up)
ego_up <- run_go(up)
ego_dn <- run_go(dn)

ggsave("figures/go_bp_dotplot.png", dotplot(ego, showCategory = 15), width = 8, height = 6, dpi = 150)
write.csv(as.data.frame(ego), "results/go_bp_enrichment.csv", row.names = FALSE)
write.csv(as.data.frame(ego2), "results/go_bp_enrichment_simplified.csv", row.names = FALSE)
write.csv(as.data.frame(ego_dn), "results/go_bp_enrichment_down_genes.csv", row.names = FALSE)

print(c(sig = length(sig), up = length(up), down = length(dn),
        terms_all = nrow(as.data.frame(ego)), terms_simplified = nrow(as.data.frame(ego2)),
        terms_up = nrow(as.data.frame(ego_up)), terms_down = nrow(as.data.frame(ego_dn))))
