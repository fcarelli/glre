# statistics for shared peaks
him17_peaks_stats = read.table("HIM17_peaks/statistics.txt")

him17_peaks_stats_modified = rbind(him17_peaks_stats[2,], him17_peaks_stats[1,] - him17_peaks_stats[2,])
him17_peaks_stats_fraction = rbind(him17_peaks_stats_modified[1,]/him17_peaks_stats[1,],
                                    him17_peaks_stats_modified[2,]/him17_peaks_stats[1,])

pdf(file = "plots/HIM17_peaks_promoters_overlap.pdf", width = 6, height = 6, useDingbats = F)

par(mfrow=c(1,1), mar=c(6, 4, 3, 2), oma=c(2,0,0,0))
barplot(as.matrix(him17_peaks_stats_fraction), las=1, ylab = "fraction of promoters", col=c("forestgreen", "grey"))
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(0, -0.9, legend = "shared peaks", fill = "forestgreen", bty="n", xjust = 0.5)

dev.off()

# venn diagrams shared peaks and m1m2 pairs
library("VennDiagram")
source("scripts/common_R_functions.R")

him17_m1m2_ovlp = read.table("HIM17_peaks/m1m2_overlap.txt")
pdf(paste("plots/HIM17_peaks_m1m2_overlap.pdf", sep=""), width = 6, height = 6)
par(mfrow=c(1,1), oma=c(1,1,1,1), mar=c(4,4,3,3))
venn.plot <- draw.pairwise.venn(area1 = him17_m1m2_ovlp[1,1], area2 = him17_m1m2_ovlp[2,1], cross.area = him17_m1m2_ovlp[3,1],
                                category = c("m1m2", "HIM17_peaks"), euler.d = T, scaled = T, fill=c("grey30", "grey80"), alpha = c(0.5, 0.5), main = "ChIP-seq peaks overlap")
grid.draw(venn.plot);
dev.off()
