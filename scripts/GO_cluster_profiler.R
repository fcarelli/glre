require(clusterProfiler)
require(enrichplot)
require(org.Ce.eg.db)

him17_direct_targets = read.table("DE_analysis/genes_downregulated_him17_vs_N2.direct.txt")

pdf(file = "plots/GO_enrichment_direct_targets.pdf", width = 10, height = 8, useDingbats = F, pointsize = 5)
ck_targets_BP <- compareCluster(geneCluster = list("cluster1"=as.character(him17_direct_targets$V1)), fun = "enrichGO", ont = "BP", pvalueCutoff = 0.001, minGSSize = 5, pAdjustMethod = "BH", keyType = "WORMBASE", OrgDb='org.Ce.eg.db')
dotplot(ck_targets_BP, showCategory = 10)
dev.off()




