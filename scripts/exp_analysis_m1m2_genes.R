# comparison of (not )m1m2-associated (not )germline-specific gene expression using RNAseq data from Boeck et al. 
source("scripts/common_R_functions.R")
library("vioplot")
stages = c("embryo_series_1", "embryo_series_2", "embryo_series_3","embryo_series_4", "postembryonic")

pdf(file = "plots/expression_m1m2_genes_elegans_development.pdf", width = 7.20, height = 9.72)
par(mar=c(7,4,3,2), mfrow=c(3,2), oma=c(2,0,0,0))

for (stage in stages){
  exp_m1m2_GL = read.table(paste("GL_genes_exp/promoters_gl_specific.elegans.m1m2.unique_promoter.", stage, ".exp", sep=""))
  exp_no_m1m2_GL = read.table(paste("GL_genes_exp/promoters_gl_specific.elegans.no_m1m2.unique_promoter.", stage, ".exp", sep=""))
  exp_m1m2_noGL = read.table(paste("GL_genes_exp/promoters_not_gl_specific.elegans.m1m2.unique_promoter.", stage, ".exp", sep=""))
  exp_no_m1m2_noGL = read.table(paste("GL_genes_exp/promoters_not_gl_specific.elegans.no_m1m2.unique_promoter.", stage, ".exp", sep=""))
  
  exp_m1m2_GL_log = log2(exp_m1m2_GL + 1e-06)
  exp_no_m1m2_GL_log = log2(exp_no_m1m2_GL + 1e-06)
  exp_m1m2_noGL_log = log2(exp_m1m2_noGL + 1e-06)
  exp_no_m1m2_noGL_log = log2(exp_no_m1m2_noGL + 1e-06)
  
  exp_GL_wilcox = c()
  exp_noGL_wilcox = c()
  for (timepoint in c(1:dim(exp_m1m2_GL)[2])){
    exp_GL_wilcox = append(exp_GL_wilcox, wilcox.test(exp_m1m2_GL_log[,timepoint], exp_no_m1m2_GL_log[,timepoint])$p.value)
    exp_noGL_wilcox = append(exp_noGL_wilcox, wilcox.test(exp_m1m2_noGL_log[,timepoint], exp_no_m1m2_noGL_log[,timepoint])$p.value)
  }

  k = boxplot(c(exp_m1m2_GL_log, exp_no_m1m2_GL_log), plot=F)
  plot_lim = c(min(k$stats), max(k$stats) + (max(k$stats) - min(k$stats))/10)
  boxplot(c(exp_m1m2_GL_log, exp_no_m1m2_GL_log), outline=F, main = paste("expression", stage),
          at=c(seq(1,dim(exp_m1m2_GL_log)[2]*2-1, by=2), seq(2,dim(exp_m1m2_GL_log)[2]*2, by=2)), notch=T, ylim = plot_lim,
          col=rep(c("darkorchid2", "grey"), each=dim(exp_m1m2_GL_log)[2]), xaxt="n", las=1, ylab="log2 dcpm")
  if (stage == "postembryonic"){
    #axis(1, at=seq(1.5,dim(exp_m1m2_GL_log)[2]*2-0.5, by=2), labels = names(exp_m1m2_GL), las=2)}
    axis(1, at=seq(1.5,dim(exp_m1m2_GL_log)[2]*2-0.5, by=2), labels = sapply(strsplit(names(exp_m1m2_GL), split="[_]"), "[[", 1), las=2)}
else{
    axis(1, at=seq(1.5,dim(exp_m1m2_GL_log)[2]*2-0.5, by=2), labels = sapply(strsplit(sapply(strsplit(names(exp_m1m2_GL), split="[.]"), "[[", 2), split="[_]"), "[[", 1), las=2)}
  abline(v=seq(2.5, dim(exp_m1m2_GL_log)[2]*2-1.5, by=2), lty=2, lwd = 2, col="grey")
  text(seq(1.5,dim(exp_m1m2_GL_log)[2]*2-0.5, by=2), max(k$stats) + (max(k$stats) - min(k$stats))/20, labels = graphical_pvalue(exp_GL_wilcox))
}


exp_m1m2_GL = read.table("GL_genes_exp/promoters_gl_specific.elegans.m1m2.unique_promoter.embryo_germline.exp")
exp_no_m1m2_GL = read.table("GL_genes_exp/promoters_gl_specific.elegans.no_m1m2.unique_promoter.embryo_germline.exp")
exp_m1m2_noGL = read.table("GL_genes_exp/promoters_not_gl_specific.elegans.m1m2.unique_promoter.embryo_germline.exp")
exp_no_m1m2_noGL = read.table("GL_genes_exp/promoters_not_gl_specific.elegans.no_m1m2.unique_promoter.embryo_germline.exp")

exp_m1m2_GL_log = log2(exp_m1m2_GL + 1e-06)
exp_no_m1m2_GL_log = log2(exp_no_m1m2_GL + 1e-06)
exp_m1m2_noGL_log = log2(exp_m1m2_noGL + 1e-06)
exp_no_m1m2_noGL_log = log2(exp_no_m1m2_noGL + 1e-06)

exp_data = c(3, 6, 9)
exp_GL_wilcox = c()
exp_noGL_wilcox = c()
for (timepoint in exp_data){
  exp_GL_wilcox = append(exp_GL_wilcox, wilcox.test(exp_m1m2_GL_log[,timepoint], exp_no_m1m2_GL_log[,timepoint])$p.value)
  exp_noGL_wilcox = append(exp_noGL_wilcox, wilcox.test(exp_m1m2_noGL_log[,timepoint], exp_no_m1m2_noGL_log[,timepoint])$p.value)
}

k = boxplot(c(exp_m1m2_GL_log[,c(3, 6, 9)], exp_no_m1m2_GL_log[,c(3, 6, 9)]), plot=F)
plot_lim = c(min(k$stats), max(k$stats) + (max(k$stats) - min(k$stats))/10)
boxplot(c(exp_m1m2_GL_log[,c(3, 6, 9)], exp_no_m1m2_GL_log[,c(3, 6, 9)]), outline=F, main = "expression embryonic germline",
        at=c(seq(1,5, by=2), seq(2,6, by=2)), notch=T, ylim = plot_lim,
        col=rep(c("darkorchid2", "grey"), each=3), xaxt="n", las=1, ylab="log2 bootstrap.TPM")
axis(1, at=c(1.5,3.5,5.5), labels = c("pseudotime_1", "pseudotime_2", "pseudotime_3"), las=2)
abline(v=c(2.5, 4.5), lty=2, lwd = 2, col="grey")
text(c(1.5,3.5,5.5), max(k$stats) + (max(k$stats) - min(k$stats))/20, labels = graphical_pvalue(exp_GL_wilcox))

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(0, -0.9, legend = c("m1m2", "no m1m2"), fill = c("darkorchid3", "grey"), bty="n", ncol = 2, xjust = 0.5)

dev.off()

# comparison of (not )m1m2-associated (not )germline-specific gene expression using RNAseq data from PGCs (Lee et al. 2017)
pgc_exp = read.table("data/external_data/rnaseq/elegans/gene_expression.pgc.txt")
pgc_exp_avg = data.frame(cbind(rowMeans(pgc_exp[,c(1:4)]), rowMeans(pgc_exp[,c(5:7)])))
pgc_exp_avg = pgc_exp_avg+1
names(pgc_exp_avg) = c("starved", "fed")

postdev_exp = read.table("data/external_data/rnaseq/elegans/gene_expression.stages.txt")
postdev_exp_avg = data.frame(cbind(rowMeans(postdev_exp[,c(1:3)]), postdev_exp[,c(4:8)]))
postdev_exp_avg = postdev_exp_avg+1
names(postdev_exp_avg) = c("L1", "L2", "L3", "L4", "lL4", "YA")


unique_m1m2_gl = read.table("promoter_annotation/promoters_gl_specific.elegans.m1m2.unique_promoter.genes")
unique_no_m1m2_gl = read.table("promoter_annotation/promoters_gl_specific.elegans.no_m1m2.unique_promoter.genes")

pdf(paste("plots/expression_unique_GL_promoter_genes.pdf", sep="."), width=14, height=6, useDingbats = F)

par(mfrow=c(1,1), mar=c(4, 4, 3, 2), oma=c(2,0,0,0))

k=vioplot(log2(pgc_exp_avg$starved[row.names(pgc_exp_avg) %in% unique_m1m2_gl$V1]),
          log2(pgc_exp_avg$starved[row.names(pgc_exp_avg) %in% unique_no_m1m2_gl$V1]),
          log2(postdev_exp_avg$lL4[row.names(postdev_exp_avg) %in% unique_m1m2_gl$V1]),
          log2(postdev_exp_avg$lL4[row.names(postdev_exp_avg) %in% unique_no_m1m2_gl$V1]),
          log2(postdev_exp_avg$YA[row.names(postdev_exp_avg) %in% unique_m1m2_gl$V1]),
          log2(postdev_exp_avg$YA[row.names(postdev_exp_avg) %in% unique_no_m1m2_gl$V1]), las = 1,
          main = "expression unique GL promoter genes", ylab = "log2 TPM+1", xaxt = "n",
          col = rep(c(t_col("darkorchid3"), t_col("grey")), 7))

all_tests_p = c()
all_tests_p = append(all_tests_p, wilcox.test(pgc_exp_avg$starved[row.names(pgc_exp_avg) %in% unique_m1m2_gl$V1],
                                              pgc_exp_avg$starved[row.names(pgc_exp_avg) %in% unique_no_m1m2_gl$V1])$p.value)
for (stage in c(5,6)){
  all_tests_p = append(all_tests_p, wilcox.test(postdev_exp_avg[row.names(postdev_exp_avg) %in% unique_m1m2_gl$V1, stage],
                                                postdev_exp_avg[row.names(postdev_exp_avg) %in% unique_no_m1m2_gl$V1, stage])$p.value)
}

axis(1, at = c(1.5, 3.5, 5.5), labels = c("PGC", "late L4", "YA"))
text(c(1.5, 3.5, 5.5), max(k$upper), labels = graphical_pvalue(p.adjust(all_tests_p, method="fdr")), pos = 1)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(0, -0.9, legend = c("m1m2", "no m1m2"), fill = c("darkorchid3", "grey"), bty="n", ncol = 2, xjust = 0.5)

dev.off()



