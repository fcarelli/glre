# CERP2 alignment of divergent motif pairs
library(ggplot2)

pdf(file = "plots/motif_pairs_repeat_consensi_aln_score.pdf", width = 8, height = 8)
par(mfrow=c(2,1), mar = c(5,4,3,2), oma = c(0,0,0,0))

promoters_cerp2_score = read.table("motif_repeat_consensus/CERP2_motifs.promoters.CERP2.scores")
inactive_cerp2_score = read.table("motif_repeat_consensus/CERP2_motifs.inactive.CERP2.scores")

promoters_cerp2_score$V3 = alpha("darkred", 0)
promoters_cerp2_score$V4 = "darkred"
inactive_cerp2_score$V3 = "grey"
inactive_cerp2_score$V4 = alpha("grey", 0)

all_cerp2_score = rbind(promoters_cerp2_score[,c(2,3,4)], inactive_cerp2_score[,c(2,3,4)])
all_cerp2_score = all_cerp2_score[order(all_cerp2_score$V2,decreasing=T),]

plot(all_cerp2_score$V2, col = all_cerp2_score$V3, pch = 19, main = "divergent pairs - CERP2 alignment", ylab = "hmmsearch score", las = 1, xlab = "")
points(all_cerp2_score$V2, col = all_cerp2_score$V4, pch = 19)
text(dim(all_cerp2_score)[1]*0.75, max(all_cerp2_score[,1])-(max(all_cerp2_score[,1])-min(all_cerp2_score[,1]))*0.1, pos = 4, labels = paste("promoters (N=", dim(promoters_cerp2_score)[1], ")", sep=""))
points(dim(all_cerp2_score)[1]*0.75, max(all_cerp2_score[,1])-(max(all_cerp2_score[,1])-min(all_cerp2_score[,1]))*0.10, col = "darkred", pch=19)

# CELE2 alignment of tandem_m2m1 motif pairs
promoters_cele2_score = read.table("motif_repeat_consensus/CELE2_motifs.promoters.CELE2.scores")
inactive_cele2_score = read.table("motif_repeat_consensus/CELE2_motifs.inactive.CELE2.scores")

promoters_cele2_score$V3 = alpha("darkred", 0)
promoters_cele2_score$V4 = "darkred"
inactive_cele2_score$V3 = "grey"
inactive_cele2_score$V4 = alpha("grey", 0)

all_cele2_score = rbind(promoters_cele2_score[,c(2,3,4)], inactive_cele2_score[,c(2,3,4)])
all_cele2_score = all_cele2_score[order(all_cele2_score$V2,decreasing=T),]

plot(all_cele2_score$V2, col = all_cele2_score$V3, pch = 19, main = "tandem_m2m1 pairs - CELE2 alignment", ylab = "hmmsearch score", las = 1, xlab = "")
points(all_cele2_score$V2, col = all_cele2_score$V4, pch = 19)
text(dim(all_cele2_score)[1]*0.75, max(all_cele2_score[,1])-(max(all_cele2_score[,1])-min(all_cele2_score[,1]))*0.1, pos = 4, labels = paste("promoters (N=", dim(promoters_cele2_score)[1], ")", sep=""))
points(dim(all_cele2_score)[1]*0.75, max(all_cele2_score[,1])-(max(all_cele2_score[,1])-min(all_cele2_score[,1]))*0.10, col = "darkred", pch=19)

dev.off()

