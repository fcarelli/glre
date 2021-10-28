# plot length of annotated CERP2 and CELE2 repeats
repeats_bed = read.table("species/elegans/repeats_dfam/elegans_dfam.repeats.id.bed")
CERP2 = repeats_bed[repeats_bed$V4 == "CERP2",]
CELE2 = repeats_bed[repeats_bed$V4 == "CELE2",]

CERP2_consensus_len = 328
CELE2_consensus_len = 325

pdf(file = "plots/CERP2_CELE2_length.pdf", width = 6, height = 6)
par(mfrow=c(2,1), mar=c(4,4,3,2), oma=c(0,0,0,0))
hist(CERP2$V3-CERP2$V2, breaks=100, main = "CERP2", ylab = "nº of repeats", las = 1, xlab = "repeat length", col = "grey")
abline(v=CERP2_consensus_len, col=2, lty=2, lwd=2)
hist(CELE2$V3-CELE2$V2, breaks=100, main = "CELE2", ylab = "nº of repeats", las = 1, xlab = "repeat length", col = "grey")
abline(v=CELE2_consensus_len, col=2, lty=2, lwd=2)
dev.off()

