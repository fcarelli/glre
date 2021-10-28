source("scripts/common_R_functions.R")

# motif3 enrichment
m3_enrichment = read.table("motif_enrichment/elegans/elegans.motif3_association.txt")
t_m3_enrichment = matrix(as.numeric(c(m3_enrichment[1,], m3_enrichment[2,])), ncol=4, byrow=F)

ylim_plot = c(0, max(colSums(t_m3_enrichment)))

pdf(file = "plots/m3_overlap.pdf", width = 8, height = 6, useDingbats = FALSE )

par(mfrow=c(1,1), mar=c(4, 4, 3, 2), oma=c(2,0,0,0))
k = barplot(t_m3_enrichment, col= c("grey40", "grey90") , ylim = c(ylim_plot[1], ylim_plot[2]*1.2),
        names.arg=c("GL-specific", "non GL-specific", "GL-specific", "non GL-specific"), las=1, ylab = "n. promoters", main = "m3 overlap")
text(k[1]+(k[2]-k[1])/2, ylim_plot[2]*1.1, labels = graphical_pvalue(fisher.test(t_m3_enrichment[,c(1,2)])$p.value), pos=3)
text(k[3]+(k[4]-k[3])/2, ylim_plot[2]*1.1, labels = graphical_pvalue(fisher.test(t_m3_enrichment[,c(3,4)])$p.value), pos=3)
segments(c(k[1], k[3]), ylim_plot[2]*1.1, c(k[2], k[4]), ylim_plot[2]*1.1)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(0, -0.9, legend = c("m3", "no m3"), fill = c("grey40", "grey90"), bty="n", ncol = 2, xjust = 0.5)

dev.off()

# GC and CpG content at C. elegans promoters

pdf(file = "plots/GC_content.pdf", width = 12, height = 12, useDingbats = FALSE )
par(cex=1, mfrow = c(2,2), cex.axis=0.8)

for (species in c("elegans", "briggsae")){
  GC_GL_m1m2 = read.table(paste("RE_features/reg_elements_all.", species, ".gl_specific.promoters.m1m2.nuc", sep=""))
  GC_nonGL_m1m2 = read.table(paste("RE_features/reg_elements_all.", species, ".not_gl_specific.promoters.m1m2.nuc", sep=""))
  GC_GL_no_m1m2 = read.table(paste("RE_features/reg_elements_all.", species, ".gl_specific.promoters.no_m1m2.nuc", sep=""))
  GC_nonGL_no_m1m2 = read.table(paste("RE_features/reg_elements_all.", species, ".not_gl_specific.promoters.no_m1m2.nuc", sep=""))

  CpG_GL_m1m2 = read.table(paste("RE_features/reg_elements_all.", species, ".gl_specific.promoters.m1m2.cpg", sep=""))
  CpG_nonGL_m1m2 = read.table(paste("RE_features/reg_elements_all.", species, ".not_gl_specific.promoters.m1m2.cpg", sep=""))
  CpG_GL_no_m1m2 = read.table(paste("RE_features/reg_elements_all.", species, ".gl_specific.promoters.no_m1m2.cpg", sep=""))
  CpG_nonGL_no_m1m2 = read.table(paste("RE_features/reg_elements_all.", species, ".not_gl_specific.promoters.no_m1m2.cpg", sep=""))

boxplot(GC_GL_m1m2$V8,
        GC_nonGL_m1m2$V8,
        GC_GL_no_m1m2$V8,
        GC_nonGL_no_m1m2$V8, notch=T, col=c("darkorchid2", "grey"), las=1, ylab = "GC content",
        names = c("GL m1m2", "nonGL m1m2", "GL nomotif", "noGL nomotif"), main = paste("GC content C.", species, "promoters"))
boxplot(CpG_GL_m1m2$V2,
        CpG_nonGL_m1m2$V2,
        CpG_GL_no_m1m2$V2,
        CpG_nonGL_no_m1m2$V2, notch=T, col=c("darkorchid2", "grey"), las=1, ylab = "CpG content",
        names = c("GL m1m2", "nonGL m1m2", "GL nomotif", "noGL nomotif"), main = paste("CpG content C.", species, "promoters"))
}
dev.off()

