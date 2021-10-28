# ATAC-seq coverage over GL-specific REs
source("scripts/common_R_functions.R")
library("vioplot")

RE_m1m2_atac_wt = read.table("atac_cov/ele_m1m2_RE_atac_adult.txt")
RE_no_m1m2_atac_wt = read.table("atac_cov/ele_no_m1m2_RE_atac_adult.txt")

prom_m1m2_atac_wt = RE_m1m2_atac_wt[grep("coding.promoter", RE_m1m2_atac_wt$V4),]
prom_no_m1m2_atac_wt = RE_no_m1m2_atac_wt[grep("coding.promoter", RE_no_m1m2_atac_wt$V4),]

RE_m1m2_atac_glp1 = read.table("atac_cov/ele_m1m2_RE_atac_glp1.txt")
RE_no_m1m2_atac_glp1 = read.table("atac_cov/ele_no_m1m2_RE_atac_glp1.txt")

prom_m1m2_atac_glp1 = RE_m1m2_atac_glp1[grep("coding.promoter", RE_m1m2_atac_glp1$V4),]
prom_no_m1m2_atac_glp1 = RE_no_m1m2_atac_glp1[grep("coding.promoter", RE_no_m1m2_atac_glp1$V4),]

RE_m1m2_atac_pgc = read.table("atac_cov/ele_m1m2_RE_atac_pgc.txt")
RE_no_m1m2_atac_pgc = read.table("atac_cov/ele_no_m1m2_RE_atac_pgc.txt")

prom_m1m2_atac_pgc = RE_m1m2_atac_pgc[grep("coding.promoter", RE_m1m2_atac_pgc$V4),]
prom_no_m1m2_atac_pgc = RE_no_m1m2_atac_pgc[grep("coding.promoter", RE_no_m1m2_atac_glp1$V4),]

RE_m1m2_atac_gl = read.table("atac_cov/ele_m1m2_RE_atac_adult_gl.txt")
RE_no_m1m2_atac_gl = read.table("atac_cov/ele_no_m1m2_RE_atac_adult_gl.txt")

prom_m1m2_atac_gl = RE_m1m2_atac_gl[grep("coding.promoter", RE_m1m2_atac_gl$V4),]
prom_no_m1m2_atac_gl = RE_no_m1m2_atac_gl[grep("coding.promoter", RE_no_m1m2_atac_gl$V4),]


wilcox_all_p_values = c()
wilcox_all_p_values = append(wilcox_all_p_values, wilcox.test(log2(prom_m1m2_atac_wt$V7),
                                                              log2(prom_no_m1m2_atac_wt$V7))$p.value)
wilcox_all_p_values = append(wilcox_all_p_values, wilcox.test(log2(prom_m1m2_atac_glp1$V7),
                                                              log2(prom_no_m1m2_atac_glp1$V7))$p.value)
wilcox_all_p_values = append(wilcox_all_p_values, wilcox.test(log2(prom_m1m2_atac_pgc$V7),
                                                              log2(prom_no_m1m2_atac_pgc$V7))$p.value)
wilcox_all_p_values = append(wilcox_all_p_values, wilcox.test(log2(prom_m1m2_atac_gl$V7),
                                                              log2(prom_no_m1m2_atac_gl$V7))$p.value)

adj_wilcox_all_p_values = p.adjust(wilcox_all_p_values, method = "fdr")

pdf(file = "plots/ATAC_coverage_GL_promoters_elegans.pdf", width = 10, height = 6, useDingbats = F)
par(mfrow = c(1,1), mar = c(5,5,3,3), oma = c(2,0,0,0))
k = vioplot(log2(prom_m1m2_atac_pgc$V7/2),
            log2(prom_no_m1m2_atac_pgc$V7/2),
            log2(prom_m1m2_atac_wt$V7/2),
            log2(prom_no_m1m2_atac_wt$V7/2),
            log2(prom_m1m2_atac_gl$V7/2),
            log2(prom_no_m1m2_atac_gl$V7/2),
            log2(prom_m1m2_atac_glp1$V7/2),
            log2(prom_no_m1m2_atac_glp1$V7/2), col = rep(c(t_col("darkorchid3"), t_col("grey")), 4), main = "GL-specific promoter accessibility", 
            las = 1, ylab = "log2(ATAC-seq read coverage)", names = F, tcl=FALSE, yaxt="n", adj = 0.5)
text(c(1.5, 3.5, 5.5, 7.5), max(k$upper), labels = graphical_pvalue(adj_wilcox_all_p_values), pos = 1, adj = 0.5)
axis(1, at = c(1.5, 3.5, 5.5, 7.5), labels = c("L1 PGCs", "adult wt",  "adult gl", "adult glp-1"), adj = 0.5)
axis(2, at = seq(0, 12, by=2), labels = seq(0, 12, by=2), las=1, adj = 1)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(0, -0.9, legend = c("m1m2", "no m1m2"), fill = c("darkorchid3", "grey"), bty="n", ncol = 2, xjust = 0.5, adj = 0)

dev.off()
