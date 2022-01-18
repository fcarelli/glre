# plotting statistics on m1m2 arrangements and intermotif distances
source("scripts/common_R_functions.R")

elegans_m1m2_all=read.table("motif_enrichment/elegans/elegans.m1m2_clusters.arrangements_all.summary", header=T)

pdf(file = "plots/m1m2_motif_arrangements.elegans.pdf", width = 7, height = 5)

par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(4,4,3,3))
k=barplot(elegans_m1m2_all$N_clusters, las=1, names.arg=elegans_m1m2_all$motif_arrangement,
          col=c("grey90", "grey60", "grey30", "grey1"), main = "C. elegans m1m2 pairs", ylab = "nº of m1m2 pairs")

dev.off()


elegans_m1m2=read.table("motif_enrichment/elegans/elegans.m1m2_clusters.arrangements_by_distance.summary", header=T)

pdf(file = "plots/m1m2_motif_distance.elegans.pdf", width = 10, height = 5)

par(mfrow=c(1,1), oma=c(2,0,0,0), mar=c(4,4,3,3), pch=21)
k=barplot(as.matrix(t(elegans_m1m2[,c(2:5)])), beside=T, las=1,
        col=c("grey90", "grey60", "grey30", "grey1"), main = "C. elegans m1-m2 distance", ylab = "nº of m1m2 pairs")
text(seq(3, 103, by=5), -50, labels = elegans_m1m2$motif_distance, srt = 45, xpd=T, cex=0.9, adj = c(0.5, 0))
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(0, -0.9, legend = c("convergent", "divergent", "tandem_m1m2", "tandem_m2m1"), fill = c("grey90", "grey60", "grey30", "grey1"), bty="n", ncol = 4, xjust = 0.5)

dev.off()


elegans_m1m2_repeats = read.table("motif_enrichment/elegans/repeat_m1m2_overlap.summary", header=T)
elegans_m1m2_repeats_filtered = elegans_m1m2_repeats[elegans_m1m2_repeats$total_repeats > 100,]
elegans_m1m2_repeats_filtered_enriched = elegans_m1m2_repeats_filtered[which(apply(elegans_m1m2_repeats_filtered, 1, function(x){sum(as.numeric(x[c(3:6)])) > as.numeric(x[2])*0.5})),]

elegans_m1m2_repeats_filtered_enrichment_matrix=apply(elegans_m1m2_repeats_filtered_enriched, 1, function(x){as.numeric(x[c(3:6)])/as.numeric(x[2])})

pdf(file = "plots/m1m2_repeat_overlap.elegans.pdf", width = 7, height = 5)

par(mfrow=c(1,1), oma=c(2,0,0,0), mar=c(4,4,3,3), pch=21)
k=barplot(elegans_m1m2_repeats_filtered_enrichment_matrix, beside=F, las=1,
          col=c("yellow", "forestgreen", "orangered", "deepskyblue", "grey70", "grey1"), main = "C. elegans m1m2 repeat overlap", ylab = "repeat fraction", ylim=c(0, 1), 
          names.arg=paste(elegans_m1m2_repeats_filtered_enriched$repeat., " (n=", elegans_m1m2_repeats_filtered_enriched$total_repeats, ")", sep=""))
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(0, -0.9, legend = c("convergent", "divergent", "tandem_m1m2", "tandem_m2m1", "single motif","no motif"), fill = c("yellow", "forestgreen", "orangered", "deepskyblue", "grey70", "grey1"), bty="n", ncol = 3, xjust = 0.5)

dev.off()

# GL_RE/repeats overlap
elegans_GL_RE_repeats = read.table("motif_enrichment/elegans/repeat_GL_RE_overlap.summary", header=T)
elegans_GL_RE_repeats_FC=c()
elegans_GL_RE_repeats_p_value=c()

for (i in c(1:dim(elegans_GL_RE_repeats)[1])){
  elegans_GL_RE_repeats_FC <-  append(elegans_GL_RE_repeats_FC, (elegans_GL_RE_repeats[i,3]/elegans_GL_RE_repeats[i,4]) / (elegans_GL_RE_repeats[i,5]/elegans_GL_RE_repeats[i,6]))
  elegans_GL_RE_repeats_p_value <- append(elegans_GL_RE_repeats_p_value, fisher.test(matrix(c(elegans_GL_RE_repeats[i,3], elegans_GL_RE_repeats[i,4]-elegans_GL_RE_repeats[i,3], elegans_GL_RE_repeats[i,5],
                                                                                              elegans_GL_RE_repeats[i,6]-elegans_GL_RE_repeats[i,5]), ncol=2, byrow = T))$p.value)
}

elegans_GL_RE_repeats_q_value = p.adjust(elegans_GL_RE_repeats_p_value, method="fdr")
elegans_GL_RE_repeats_enrichment = cbind(elegans_GL_RE_repeats, elegans_GL_RE_repeats_FC, elegans_GL_RE_repeats_p_value, elegans_GL_RE_repeats_q_value)
elegans_GL_RE_repeats_enriched = elegans_GL_RE_repeats_enrichment[elegans_GL_RE_repeats_enrichment$elegans_GL_RE_repeats_q_value < 0.001 & elegans_GL_RE_repeats_enrichment$elegans_GL_RE_repeats_FC > 1 & is.finite(elegans_GL_RE_repeats_enrichment$elegans_GL_RE_repeats_FC),]
elegans_GL_RE_repeats_enriched = elegans_GL_RE_repeats_enriched[order(elegans_GL_RE_repeats_enriched$elegans_GL_RE_repeats_q_value),]

write.table(elegans_GL_RE_repeats_enrichment, sep="\t", file = "motif_enrichment/elegans/repeat_GL_RE_overlap.enrichment", quote = F, row.names = F, col.names=T)

q_val_col = heat.colors(n = 10)
y <- seq(ceiling(log10(min(elegans_GL_RE_repeats_enriched$elegans_GL_RE_repeats_q_value))), 1, len = 10)

elegans_GL_RE_repeats_enriched$elegans_GL_RE_repeats_enriched_col_bin = NA
for (i in seq(10, 1, -1)){
  elegans_GL_RE_repeats_enriched$elegans_GL_RE_repeats_enriched_col_bin[which(log10(elegans_GL_RE_repeats_enriched$elegans_GL_RE_repeats_q_value) <= y[i])] = q_val_col[i]
}

pdf("plots/GL_RE_repeat_overlap.elegans.pdf", width=12, height=8)

layout(matrix(1:2,nrow=1),widths=c(0.9,0.1))

par(mar=c(12,4,3,3))
k=barplot(elegans_GL_RE_repeats_enriched$elegans_GL_RE_repeats_FC, las=2,
          col=elegans_GL_RE_repeats_enriched$elegans_GL_RE_repeats_enriched_col_bin, main = "C. elegans GL-specific RE repeat overlap", ylab = "GL RE/repeats overlap enrichment", 
          names.arg=paste(elegans_GL_RE_repeats_enriched$repeat_class, " (n=", elegans_GL_RE_repeats_enriched$repeat_number, ")", sep=""))
text(k[,1], elegans_GL_RE_repeats_enriched$elegans_GL_RE_repeats_FC, labels = elegans_GL_RE_repeats_enriched$rep_in_GL_specific_peaks, pos = 1)
abline(h=1, lty=2, lwd=2, col=2)

xl <- 1
yb <- 1
xr <- 1.5
yt <- 2

par(mar=c(5.1,0.5,4.1,0.5))
plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(1,2),xaxt="n",yaxt="n",bty="n")
rect(
  xl,
  head(seq(yb,yt,(yt-yb)/10),-1),
  xr,
  tail(seq(yb,yt,(yt-yb)/10),-1),
  col=q_val_col
)

mtext(ceiling(y),side=2,at=tail(seq(yb,yt,(yt-yb)/10),-1)-0.05,las=2,cex=0.7)
mtext("log10(P-value)",side=3,las=1,cex=0.7)
dev.off()




# GL_promoters/repeats overlap
elegans_GL_promoter_repeats = read.table("motif_enrichment/elegans/repeat_GL_promoter_overlap.summary", header=T)
elegans_GL_promoter_repeats_FC=c()
elegans_GL_promoter_repeats_p_value=c()

for (i in c(1:dim(elegans_GL_promoter_repeats)[1])){
  elegans_GL_promoter_repeats_FC <-  append(elegans_GL_promoter_repeats_FC, (elegans_GL_promoter_repeats[i,3]/elegans_GL_promoter_repeats[i,4]) / (elegans_GL_promoter_repeats[i,5]/elegans_GL_promoter_repeats[i,6]))
  elegans_GL_promoter_repeats_p_value <- append(elegans_GL_promoter_repeats_p_value, fisher.test(matrix(c(elegans_GL_promoter_repeats[i,3], elegans_GL_promoter_repeats[i,4]-elegans_GL_promoter_repeats[i,3], elegans_GL_promoter_repeats[i,5],
                                                                                              elegans_GL_promoter_repeats[i,6]-elegans_GL_promoter_repeats[i,5]), ncol=2, byrow = T))$p.value)
}

elegans_GL_promoter_repeats_q_value = p.adjust(elegans_GL_promoter_repeats_p_value, method="fdr")
elegans_GL_promoter_repeats_enrichment = cbind(elegans_GL_promoter_repeats, elegans_GL_promoter_repeats_FC, elegans_GL_promoter_repeats_p_value, elegans_GL_promoter_repeats_q_value)
elegans_GL_promoter_repeats_enriched = elegans_GL_promoter_repeats_enrichment[elegans_GL_promoter_repeats_enrichment$elegans_GL_promoter_repeats_q_value < 0.001 & elegans_GL_promoter_repeats_enrichment$elegans_GL_promoter_repeats_FC > 1 & is.finite(elegans_GL_promoter_repeats_enrichment$elegans_GL_promoter_repeats_FC),]
elegans_GL_promoter_repeats_enriched = elegans_GL_promoter_repeats_enriched[order(elegans_GL_promoter_repeats_enriched$elegans_GL_promoter_repeats_q_value),]

write.table(elegans_GL_promoter_repeats_enrichment, sep="\t", file = "motif_enrichment/elegans/repeat_GL_promoter_overlap.enrichment", quote = F, row.names = F, col.names=T)

q_val_col = heat.colors(n = 10)
y <- seq(ceiling(log10(min(elegans_GL_promoter_repeats_enriched$elegans_GL_promoter_repeats_q_value))), 1, len = 10)

elegans_GL_promoter_repeats_enriched$elegans_GL_promoter_repeats_enriched_col_bin = NA
for (i in seq(10, 1, -1)){
  elegans_GL_promoter_repeats_enriched$elegans_GL_promoter_repeats_enriched_col_bin[which(log10(elegans_GL_promoter_repeats_enriched$elegans_GL_promoter_repeats_q_value) <= y[i])] = q_val_col[i]
}


pdf("plots/GL_promoter_repeat_overlap.elegans.pdf", width=12, height=8)

layout(matrix(1:2,nrow=1),widths=c(0.9,0.1))

par(mar=c(12,4,3,3))
k=barplot(elegans_GL_promoter_repeats_enriched$rep_in_GL_specific_peaks, las=2,
          col=elegans_GL_promoter_repeats_enriched$elegans_GL_promoter_repeats_enriched_col_bin, main = "C. elegans GL-specific promoter repeat overlap", ylab = "GL promoter/repeats overlap enrichment",
          names.arg=paste(elegans_GL_promoter_repeats_enriched$repeat_class, " (n=", elegans_GL_promoter_repeats_enriched$repeat_number, ")", sep=""))

xl <- 1
yb <- 1
xr <- 1.5
yt <- 2

par(mar=c(5.1,0.5,4.1,0.5))
plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(1,2),xaxt="n",yaxt="n",bty="n")
rect(
  xl,
  head(seq(yb,yt,(yt-yb)/10),-1),
  xr,
  tail(seq(yb,yt,(yt-yb)/10),-1),
  col=q_val_col
)

mtext(ceiling(y),side=2,at=tail(seq(yb,yt,(yt-yb)/10),-1)-0.05,las=2,cex=0.7)
mtext("log10(P-value)",side=3,las=1,cex=0.7)
dev.off()


#pdf("plots/GL_promoter_repeat_overlap.elegans.pdf", width=12, height=8)
#layout(matrix(1:2,nrow=1),widths=c(0.9,0.1))
#par(mar=c(12,4,3,3))
#k=barplot(elegans_GL_promoter_repeats_enriched$elegans_GL_promoter_repeats_FC, las=2,
#          col=elegans_GL_promoter_repeats_enriched$elegans_GL_promoter_repeats_enriched_col_bin, main = "C. elegans GL-specific promoter repeat overlap", ylab = "GL promoter/repeats overlap enrichment",
#          names.arg=paste(elegans_GL_promoter_repeats_enriched$repeat_class, " (n=", elegans_GL_promoter_repeats_enriched$repeat_number, ")", sep=""))
#text(k[,1], elegans_GL_promoter_repeats_enriched$elegans_GL_promoter_repeats_FC, labels = elegans_GL_promoter_repeats_enriched$rep_in_GL_specific_peaks, pos = 1)
#abline(h=1, lty=2, lwd=2, col=2)
#xl <- 1
#yb <- 1
#xr <- 1.5
#yt <- 2
#par(mar=c(5.1,0.5,4.1,0.5))
#plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(1,2),xaxt="n",yaxt="n",bty="n")
#rect(
#  xl,
#  head(seq(yb,yt,(yt-yb)/10),-1),
#  xr,
#  tail(seq(yb,yt,(yt-yb)/10),-1),
#  col=q_val_col
#)
#mtext(ceiling(y),side=2,at=tail(seq(yb,yt,(yt-yb)/10),-1)-0.05,las=2,cex=0.7)
#mtext("log10(P-value)",side=3,las=1,cex=0.7)
#dev.off()
