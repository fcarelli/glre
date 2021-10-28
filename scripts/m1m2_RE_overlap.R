# statistics on RE overlap with m1m2 pairs
elegans_stats = read.table("RE_features/elegans.motif_GL_overlap.summary")
briggsae_stats = read.table("RE_features/briggsae.motif_GL_overlap.summary")


pdf(file = "plots/RE_m1m2_overlap.pdf", width = 8, height = 8)
par(cex=1, mfrow = c(1,2), cex.axis=0.8, oma=c(2,0,0,0), mar=c(4,5,3,3))

k = barplot(c(elegans_stats[1,2]/elegans_stats[1,1],
          elegans_stats[2,2]/elegans_stats[2,1],
          briggsae_stats[1,2]/briggsae_stats[1,1],
          briggsae_stats[2,2]/briggsae_stats[2,1]), las=1, col=c("darkorchid2", "grey"),
          main = "m1m2 overlap all REs", xaxt = "n", ylab = "fraction of REs")
axis(1, at=c(k[1,1]+(k[2,1] - k[1,1])/2, k[3,1]+(k[4,1] - k[3,1])/2), labels = c("C.elegans", "C.briggsae"))          
text(k, c(elegans_stats[1,2]/elegans_stats[1,1],
          elegans_stats[2,2]/elegans_stats[2,1],
          briggsae_stats[1,2]/briggsae_stats[1,1],
          briggsae_stats[2,2]/briggsae_stats[2,1]), 
     labels=c(elegans_stats[1,2], elegans_stats[2,2], briggsae_stats[1,2], briggsae_stats[2,2]), pos=1)
barplot(c(elegans_stats[1,8]/elegans_stats[1,7],
          elegans_stats[2,8]/elegans_stats[2,7],
          briggsae_stats[1,8]/briggsae_stats[1,7],
          briggsae_stats[2,8]/briggsae_stats[2,7]), las=1, col=c("darkorchid2", "grey"),
        main = "m1m2 overlap promoters", xaxt = "n", ylab = "fraction of promoters")
axis(1, at=c(k[1,1]+(k[2,1] - k[1,1])/2, k[3,1]+(k[4,1] - k[3,1])/2), labels = c("C.elegans", "C.briggsae"))          
text(k, c(elegans_stats[1,8]/elegans_stats[1,7],
          elegans_stats[2,8]/elegans_stats[2,7],
          briggsae_stats[1,8]/briggsae_stats[1,7],
          briggsae_stats[2,8]/briggsae_stats[2,7]), 
     labels=c(elegans_stats[1,8], elegans_stats[2,8], briggsae_stats[1,8], briggsae_stats[2,8]), pos=1)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(0, -0.9, legend = c("GL-specific", "nonGL-specific"), fill = c("darkorchid2", "grey"), bty="n", ncol = 2, xjust = 0.5)
dev.off()

# stacked barplot
elegans_stats_m1m2_GL = c(elegans_stats[1,1] - elegans_stats[1,2], elegans_stats[1,2], elegans_stats[1,7]-elegans_stats[1,8], elegans_stats[1,8])
briggsae_stats_m1m2_GL = c(briggsae_stats[1,1] - briggsae_stats[1,2], briggsae_stats[1,2], briggsae_stats[1,7]-briggsae_stats[1,8], briggsae_stats[1,8])

pdf(file = "plots/RE_m1m2_overlap.ele_bri.pdf", width = 8, height = 8)
par(cex=1, mfrow = c(1,1), cex.axis=0.8, oma=c(2,0,0,0), mar=c(4,5,3,3))
k = barplot(matrix(c(elegans_stats_m1m2_GL[c(2,1)], briggsae_stats_m1m2_GL[c(2,1)],
                     elegans_stats_m1m2_GL[c(4,3)], briggsae_stats_m1m2_GL[c(4,3)]), ncol=4, byrow = F), las=1, col=c("darkorchid2", "grey"),
            main = "m1m2 overlap all REs", xaxt = "n", ylab = "fraction of REs")
axis(1, at=k, label=rep(c("C. elegans", "C. briggsae"), 2))
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(0, -0.9, legend = c("m1m2", "no m1m2"), fill = c("darkorchid2", "grey"), bty="n", ncol = 2, xjust = 0.5)
dev.off()

# same plot, with info on m1m2 orientation
pdf(file = "plots/RE_m1m2_overlap.m1m2_arrangement.pdf", width = 8, height = 8)
par(cex=1, mfrow = c(1,1), cex.axis=0.8, oma=c(0,0,0,15), mar=c(4,5,3,3))

GL_RE_m1m2_arrangement = matrix(rev(c(elegans_stats[1,7]-elegans_stats[1,8], elegans_stats[1,9], elegans_stats[1,10], elegans_stats[1,11], elegans_stats[1,12],
                                  elegans_stats[1,1]-elegans_stats[1,2], elegans_stats[1,3], elegans_stats[1,4], elegans_stats[1,5], elegans_stats[1,6])), ncol = 5, byrow = T)
k = barplot(t(GL_RE_m1m2_arrangement), beside=F, las=1, col=c("deepskyblue", "orangered", "forestgreen", "gold", "grey"),
            main = "m1m2 - RE overlap", xaxt = "n", ylab = "number of REs")
axis(1, at=c(k[1], k[2]), labels = c("all REs", "promoters"))          

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(0, 0, legend = c("no m1m2", "convergent_m1m2", "divergent_m1m2", "tandem_m1m2", "tandem_m2m1"), fill = rev(c("deepskyblue", "orangered", "forestgreen", "gold", "grey")), bty="n", ncol = 1, xjust = -1, yjust = 0)
dev.off()
