# plot identity heatmap of HIM-17 orthologs similarity
library(gplots)
library(RColorBrewer)
library(ape)

par(mfrow=c(1,2))
him17_aln = read.table("TF_evolution/orthologs/HIM17.all.dist", row.names = 1)
names(him17_aln) = rownames(him17_aln)

him17_aln_mtx <- as.matrix(him17_aln)
pdf(file = "plots/HIM17_identity.heat.pdf", width = 10, height = 10)
heatmap.2(him17_aln_mtx, dendrogram = "row", breaks = 50, cellnote = ceiling(him17_aln_mtx), notecol=1, notecex = 1, 
          trace="none", srtCol=45, keysize=0.9, density.info="none", margins=c(5,8), key.title=NA, key.xlab = "% identity",
          key.par=list(mar=c(5,3,3,1)), lmat = rbind(c(0,3),c(2,1),c(4,0)), lhei = c(0.3,4,0.7), lwid = c(1.5,4))
dev.off()


# plot pairwise identity of HIM-17 orthologs between C. elegans and other related species
ele_ino_ortho = read.table("data/external_data/elegans.inopinata.1to1.txt")
ele_bri_ortho = read.table("data/external_data/elegans.briggsae.1to1.txt")
ele_nig_ortho = read.table("data/external_data/elegans.nigoni.1to1.txt")
ele_rem_ortho = read.table("data/external_data/elegans.remanei.1to1.txt")

other_species = c("inopinata", "briggsae", "nigoni", "remanei")
HIM17_ortho = c()
for (species in other_species){
  ele_sp_HIM17 = read.table(paste("TF_evolution/orthologs/HIM17.", species, ".aln", sep=""))
  if (ele_sp_HIM17[1,11] < 0.0001){
    HIM17_ortho = append(HIM17_ortho, ele_sp_HIM17[1,3])}
  else{
    HIM17_ortho = append(HIM17_ortho, NA)}
}

pdf(file = "plots/HIM17_orthologs_identity.pdf", width = 8, height = 8, useDingbats = FALSE)
par(mfrow=c(1,2), mar=c(6,4,3,2))
k=boxplot(ele_ino_ortho$V3,
        ele_bri_ortho$V3,
        ele_nig_ortho$V3,
        ele_rem_ortho$V3, notch=T, las = 1, xaxt = "n", 
        col = "grey", main = "C. elegans 1-to-1 orthologs identity", ylab = "%id")
axis(1, at = c(1:4), labels = FALSE)
text(c(1:4), par("usr")[3] - 0.2, labels = c("C. inopinata", "C. briggsae", "C. nigoni", "C. remanei"), srt = 45, pos = 1, xpd = TRUE, offset = 2)
points(c(1:4), HIM17_ortho, pch=17, col=3, cex=2)


bri_nig_ortho = read.table("data/external_data/briggsae.nigoni.1to1.txt")
bri_rem_ortho = read.table("data/external_data/briggsae.remanei.1to1.txt")
bri_ele_ortho = read.table("data/external_data/briggsae.elegans.1to1.txt")
bri_ino_ortho = read.table("data/external_data/briggsae.inopinata.1to1.txt")


other_species = c("nigoni", "remanei", "elegans", "inopinata")
HIM17_ortho = c()

HIM17_ortho = append(HIM17_ortho, bri_nig_ortho[bri_nig_ortho$V1 == "WBGene00031009", 3])
HIM17_ortho = append(HIM17_ortho, bri_rem_ortho[bri_rem_ortho$V1 == "WBGene00031009", 3])
HIM17_ortho = append(HIM17_ortho, bri_ele_ortho[bri_ele_ortho$V1 == "WBGene00031009", 3])
HIM17_ortho = append(HIM17_ortho, bri_ino_ortho[bri_ino_ortho$V1 == "WBGene00031009", 3])

k=boxplot(bri_nig_ortho$V3,
          bri_rem_ortho$V3,
          bri_ele_ortho$V3,
          bri_ino_ortho$V3, notch=T, las = 1, xaxt = "n", 
          col = "grey", main = "C. briggsae 1-to-1 orthologs identity", ylab = "%id")
axis(1, at = c(1:4), labels = FALSE)
text(c(1:4), par("usr")[3] - 0.2, labels = c("C. nigoni", "C. remanei", "C. elegans", "C. inopinata"), srt = 45, pos = 1, xpd = TRUE, offset = 2)
points(c(1:4), HIM17_ortho, pch=17, col=3, cex=2)

dev.off()

# plot THAP domain score and location across different species
all_species = c("elegans", "inopinata", "briggsae", "nigoni", "sinica", "remanei", "becei", "panamensis", "afra", "uteleia", "quiockensis", "virilis", "bovis", "plicata", "monodelphis", "contortus", 
                "placei", "polygyrus", "brasiliensis", "viviparus", "caninum", "ceylanicum", "tipulae", "rhodensis")
all_species_full = c("C. elegans", "C. inopinata", "C. briggsae", "C. nigoni", "C. sinica", "C. remanei", "C. becei", "C. panamensis", "C. afra", "C. uteleia", "C. quiockensis", "C. virilis", "C. bovis", "C. plicata", "C. monodelphis", "H. contortus",
                     "H. placei", "H. polygyrus", "N. brasiliensis", "D. viviparus", "A. caninum", "A. ceylanicum", "O. tipulae", "A. rhodensis")

cat("(((((((((C. elegans,C. inopinata),(((C. briggsae,C. nigoni),C. sinica),C. remanei)),((C. becei,C. panamensis),C.afra)),",
    "C. uteleia),C. quiockensis),C. virilis),(C.bovis,C. plicata)), C.monodelphis), (((((H. contortus, H. placei), (H. polygyrus, N. brasiliensis)), D. viviparus), (A. caninum, A. ceylanicum)), (O. tipulae, A. rhodensis)));", 
    file = "ex.tre", sep = "\n")

min_score = 0
max_score = 0
prot_length = c()
for (species in c(1: length(all_species))){
  thap = read.table(paste("TF_evolution/THAP_conservation/", all_species[species], ".THAP.parsed", sep = ""))
  min_score = min(min_score, thap$V3)
  max_score = max(max_score, thap$V3)
  prot_length=append(prot_length, thap[1,6])
}

q_val_col = rev(brewer.pal(10, "RdBu"))
y <- seq(min_score, max_score, len = 10)


pdf(file = "plots/THAP_domains_HIM17.pdf", width = 10, height = 7)
layout(matrix(1:3,nrow=1),widths=c(0.2,0.7,0.1))

tree.worms <- read.tree("ex.tre")
par(mar = c(5, 4, 3, 1))
plot(tree.worms, show.tip.label = FALSE)

unlink("ex.tre") # delete the file "ex.tre"

par(mar = c(5, 6, 3, 3))
plot(-100, -10, xlim = c(0, max(prot_length)), ylim = c(1, length(all_species)), 
     yaxt = "n", ylab = "", xlab = "HIM-17 AA position", bty="n", main = "THAP domains in HIM-17 orthologs")
segments(1, c(1: length(all_species)), prot_length, c(1: length(all_species)), lty = 2)
axis(2, at = c(1:length(all_species)), labels = all_species_full, las = 2, font=3, hadj = 1)
for (species in c(1: length(all_species))){
  thap = read.table(paste("TF_evolution/THAP_conservation/", all_species[species], ".THAP.parsed", sep = ""))
  thap$color = NA
  for (i in seq(10, 1, -1)){
    thap$color[which(thap$V3 <= y[i])] = q_val_col[i]
  }
  rect(xleft = thap$V4, ybottom = species-0.4, xright = thap$V5, ytop = species+0.4, col = 0, lwd = 0)
  rect(xleft = thap$V4, ybottom = species-0.4, xright = thap$V5, ytop = species+0.4, col = thap$color)
}


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

mtext(ceiling(y),side=2,at=tail(seq(yb,yt,(yt-yb)/10),-1)-0.05,las=2,cex=0.8)
mtext("hmmsearch score",side=3,las=1,cex=0.7)
dev.off()


