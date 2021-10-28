# plot fraction of divergent and tandem_m2m1 pairs at conserved and species-specific promoters
source("scripts/common_R_functions.R")
m1m2_conservation_stats=read.table("RE_conservation_all/conserved_vs_species_specific_promoters.m1m2_statistics.txt")

m1m2_conservation_stats_fraction = t(apply(m1m2_conservation_stats, 1, function(x){x/sum(x)}))

pdf(file = "plots/conserved_sp_specific_promoters.m1m2_overlap.pdf", width = 4, height = 8, useDingbats = F)

par(mfrow=c(2,1), mar=c(4, 4, 2, 2), oma=c(2,0,0,0))
barplot(m1m2_conservation_stats_fraction[c(1,2),c(2,4)], beside=T, las=1, main="C. elegans", ylab = "promoter fraction", col=c("grey30", "grey90"))
text(3.5, max(m1m2_conservation_stats_fraction[c(1,2),c(2,4)]), labels = graphical_pvalue(fisher.test(m1m2_conservation_stats[c(1,2),c(2,4)])$p.value), pos = 1)
barplot(m1m2_conservation_stats_fraction[c(3,4),c(2,4)], beside=T, las=1, main="C. briggsae", ylab = "promoter fraction", col=c("grey30", "grey90"))
text(3.5, max(m1m2_conservation_stats_fraction[c(1,2),c(2,4)]), labels = graphical_pvalue(fisher.test(m1m2_conservation_stats[c(3,4),c(2,4)])$p.value), pos = 1)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(0, -0.9, legend = c("conserved", "species-specific"), fill = c("grey30", "grey90"), bty="n", xjust = 0.5, ncol = 2)

dev.off()
