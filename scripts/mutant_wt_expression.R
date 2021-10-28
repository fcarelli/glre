# compare expression in him-17 mutants and wt of GL-specific genes regulated by unique promoter
# compare m1m2 vs non-m1m2 genes
source("scripts/common_R_functions.R")
library(Hmisc)

him17_DESeq = read.table("DE_analysis/genes_DESeq_him17_vs_N2.txt")

gl_m1m2_genes = read.table("gene_annotation/elegans.GL_specific_genes.m1m2_unique.genes")
gl_no_m1m2_genes = read.table("gene_annotation/elegans.GL_specific_genes.no_m1m2_unique.genes")

gl_m1m2_genes_him17_DESeq = him17_DESeq[row.names(him17_DESeq) %in% gl_m1m2_genes$V1,]
gl_no_m1m2_genes_him17_DESeq = him17_DESeq[row.names(him17_DESeq) %in% gl_no_m1m2_genes$V1,]

him17_ttest = wilcox.test(gl_m1m2_genes_him17_DESeq$log2FoldChange,
                     gl_no_m1m2_genes_him17_DESeq$log2FoldChange)

pdf(file = "plots/mut_vs_wt_expression.m1m2_promoters.pdf", width = 5, height = 5, useDingbats = F)

par(mfrow=c(1,1), oma = c(0,0,0,0), mar=c(5, 4, 3, 3))
k=boxplot(gl_m1m2_genes_him17_DESeq$log2FoldChange,
        gl_no_m1m2_genes_him17_DESeq$log2FoldChange,
        notch=T, col = c("darkorchid1", "mediumpurple4"), main = "him-17 vs wt expression", 
        las=1, ylab = "LCF mut vs wt", names=c("m1m2", "not m1m2"))
abline(h=0, col="grey", lwd=2, lty=2)
text(1.5, max(k$stats, k$out), labels = graphical_pvalue(him17_ttest$p.value))

dev.off()

# compare ubiquitously expressed genes
him17_bound = read.table("DE_analysis/elegans_HIM17.bound_genes.txt")

jacques_gene_annotation = read.table("data/external_data/gene_tissues_annotation_serizay.txt")
ubiq_genes = jacques_gene_annotation[jacques_gene_annotation$V2 == "Ubiq.",]

not_gl_m1m2_genes = read.table("promoter_annotation/promoters_not_gl_specific.elegans.m1m2.any_promoter.genes")
gl_m1m2_genes = read.table("promoter_annotation/promoters_gl_specific.elegans.m1m2.any_promoter.genes")

ubiq_m1m2 = not_gl_m1m2_genes[not_gl_m1m2_genes$V1 %nin% gl_m1m2_genes$V1 & not_gl_m1m2_genes$V1 %in% ubiq_genes$V1 & not_gl_m1m2_genes$V1 %in% him17_bound$V1,]
ubiq_not_m1m2 = ubiq_genes[ubiq_genes$V1 %nin% not_gl_m1m2_genes$V1 & ubiq_genes$V1 %nin% gl_m1m2_genes$V1 & ubiq_genes$V1 %nin% him17_bound$V1, 1]

ubiq_m1m2_genes_him17_DESeq = him17_DESeq[row.names(him17_DESeq) %in% ubiq_m1m2 & !is.na(him17_DESeq$log2FoldChange),]
ubiq_not_m1m2_genes_him17_DESeq = him17_DESeq[row.names(him17_DESeq) %in% ubiq_not_m1m2 & !is.na(him17_DESeq$log2FoldChange),]


him17_ttest = wilcox.test(ubiq_m1m2_genes_him17_DESeq$log2FoldChange,
                          ubiq_not_m1m2_genes_him17_DESeq$log2FoldChange)

pdf(file = "plots/mut_vs_wt_expression.m1m2_promoters.ubiq_genes.pdf", width = 5, height = 5, useDingbats = F)

par(mfrow=c(1,1), oma = c(0,0,0,0), mar=c(5, 4, 3, 3))
k=boxplot(ubiq_m1m2_genes_him17_DESeq$log2FoldChange,
          ubiq_not_m1m2_genes_him17_DESeq$log2FoldChange,
          notch=T, col = c("darkorchid1", "mediumpurple4"), main = "him-17 vs wt expression", 
          las=1, ylab = "LCF mut vs wt", names=c(paste("m1m2", dim(ubiq_m1m2_genes_him17_DESeq)[1], sep = "\n"),
                                                 paste("no m1m2", dim(ubiq_not_m1m2_genes_him17_DESeq)[1], sep = "\n")))
abline(h=0, col="grey", lwd=2, lty=2)
text(1.5, max(k$stats, k$out), labels = graphical_pvalue(him17_ttest$p.value))

dev.off()


# compare bound vs unbound genes (GL-specific/unique promoter)
GL_unique = c(as.character(gl_m1m2_genes$V1), as.character(gl_no_m1m2_genes$V1))

genes_HIM17_bound_DESeq = him17_DESeq[row.names(him17_DESeq) %in% GL_unique & row.names(him17_DESeq) %in% him17_bound$V1,]
genes_HIM17_not_bound_DESeq = him17_DESeq[row.names(him17_DESeq) %in% GL_unique & row.names(him17_DESeq) %nin% him17_bound$V1,]

him17_ttest = wilcox.test(genes_HIM17_bound_DESeq$log2FoldChange,
                     genes_HIM17_not_bound_DESeq$log2FoldChange)

pdf(file = "plots/mut_vs_wt_expression.bound_promoters.pdf", width = 5, height = 5, useDingbats = F)

par(mfrow=c(1,1))
k=boxplot(genes_HIM17_bound_DESeq$log2FoldChange,
          genes_HIM17_not_bound_DESeq$log2FoldChange,
          notch=T, col = c("darkorchid1", "mediumpurple4"), main = "him-17 vs wt expression", 
          las=1, ylab = "LCF mut vs wt", names=c("bound", "not bound"))
abline(h=0, col="grey", lwd=2, lty=2)
text(1.5, max(k$stats, k$out), labels = graphical_pvalue(him17_ttest$p.value))

dev.off()

