# plotting overlap of m1m2 and REs from current and Serizay annotation
RE_m1m2_overlap = read.table("RE_annotation_comparison/current_vs_serizay_RE_annotation.txt")

glre_stats = RE_m1m2_overlap[,c(1,2)]
glre_stats = cbind(glre_stats, glre_stats$GLRE-glre_stats$GLRE_m1m2)

prom_stats = RE_m1m2_overlap[,c(3,4)]
prom_stats = cbind(prom_stats, prom_stats$GL_promoters-prom_stats$GL_promoters_m1m2)

pdf(file = "plots/current_vs_serizay_RE_annotation.m1m2.pdf", width = 8, height = 6, useDingbats = F)
par(mfrow=c(1,2), mar = c(4,4,3,2), oma=c(2,0,0,0))
barplot(as.matrix(t(glre_stats[,c(2,3)])), beside = F, las = 1, names.arg=c("current", "Serizay", "shared"),  main = "all GLRE", col = c("darkorchid3", "grey"))
barplot(as.matrix(t(prom_stats[,c(2,3)])), beside = F, las = 1, names.arg=c("current", "Serizay", "shared"), main = "GL promoters", col = c("darkorchid3", "grey"))
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(0, -0.9, legend = c("m1m2", "no m1m2"), fill = c("darkorchid3", "grey"), bty="n", ncol = 2, xjust = 0.5, adj = 0)
dev.off()

# retrieving info on Serizay REs not annotated as GL-specific
current_only_in_serizay = read.table("RE_annotation_comparison/current_specific.GLRE.annot_in_serizay.txt")

pdf(file = "plots/current_vs_serizay_RE_annotation.current_specific_in_Serizay.pdf", width = 8, height = 6, useDingbats = F)
par(mfrow=c(1,2), mar = c(15,4,3,2), oma=c(0,0,0,0))
barplot(table(current_only_in_serizay$V1), beside = T, las = 2, main = "current-spec. in Serizay: highest signal")
barplot(table(current_only_in_serizay$V6), beside = F, las = 2, main = "current-spec. in Serizay: annotation")
dev.off()

# plotting annotation-specific REs absent in other annotation
absent_REs = read.table("RE_annotation_comparison/annot_specific.GLRE.absent_in_other_annot.txt")
absent_REs = cbind(absent_REs, absent_REs$specific-absent_REs$absent_in_other)

pdf(file = "plots/current_vs_serizay_RE_annotation.absent_in_other.pdf", width = 6, height = 6, useDingbats = F)
par(mfrow=c(1,1), mar = c(4,4,3,2), oma=c(2,0,0,0))
barplot(as.matrix(t(absent_REs[,c(2,3)])), beside = F, las = 1, main = "all GLRE", names.arg=c("current", "Serizay"))
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(0, -0.9, legend = c("absent in other annotation", "shared with other annotation"), fill = c("grey30", "grey"), bty="n", ncol = 2, xjust = 0.5, adj = 0)
dev.off()
