# intermotif distance in different species
m1m2_ele = read.table("motif_enrichment/elegans/elegans.m1m2_clusters.arrangements_by_distance.summary", header=T)
m1m2_bri = read.table("motif_enrichment/briggsae/briggsae.m1m2_clusters.arrangements_by_distance.summary", header=T)
m1m2_con = read.table("motif_enrichment/contortus/contortus.m1m2_clusters.arrangements_by_distance.summary", header=T)
m1m2_cey = read.table("motif_enrichment/ceylanicum/ceylanicum.m1m2_clusters.arrangements_by_distance.summary", header=T)
m1m2_pac = read.table("motif_enrichment/pacificus/pacificus.m1m2_clusters.arrangements_by_distance.summary", header=T)
m1m2_bec = read.table("motif_enrichment/becei/becei.m1m2_clusters.arrangements_by_distance.summary", header=T)
m1m2_mon = read.table("motif_enrichment/monodelphis/monodelphis.m1m2_clusters.arrangements_by_distance.summary", header=T)

pdf(file = "plots/m1m2_divergent_distance.pdf", width = 8, height = 5, useDingbats = F, pointsize = 5)
par(mfrow=c(1,1), oma=c(2,0,0,0), mar=c(8,4,3,2))
barplot(t(matrix(c(m1m2_ele$divergent_cluster,
               m1m2_bri$divergent_cluster,
               m1m2_con$divergent_cluster,
               m1m2_cey$divergent_cluster,
               m1m2_pac$divergent_cluster), ncol=5, byrow = F)), beside=T, 
        col=c("deepskyblue1", "deepskyblue4", "chocolate1", "chocolate4", "darkred"), las=1, main = "divergent pairs", ylab="nº of pairs",
        names=c(seq(10,30)), xlab = "m1-m2 distance")
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(0, -0.9, legend = c("C. elegans", "C. briggsae", "H. contortus", "A. ceylanicum", "P. pacificus"), fill = c("deepskyblue1", "deepskyblue4", "chocolate1", "chocolate4", "darkred"), bty="n", xjust = 0.5, ncol = 5)
dev.off()

pdf(file = "plots/m2m1_tandem_distance.pdf", width = 8, height = 5, useDingbats = F, pointsize = 5)
par(mfrow=c(1,1), oma=c(2,0,0,0), mar=c(8,4,3,2))
barplot(t(matrix(c(m1m2_ele$tandem_m2m1_cluster,
                   m1m2_bec$tandem_m2m1_cluster,
                   m1m2_mon$tandem_m2m1_cluster), ncol=3, byrow = F)), beside=T, 
        col=c("deepskyblue1", "grey80", "grey20"), las=1, main = "tandem_m2m1 pairs", ylab="nº of pairs",
        names=c(seq(10,30)), xlab = "m1-m2 distance")
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(0, -0.9, legend = c("C. elegans", "C. becei", "C. monodelphis"), fill = c("deepskyblue1", "grey80", "grey20"), bty="n", xjust = 0.5, ncol = 3)
dev.off()


