# plot re-annotated CERP2/CELE2 numbers in different species alongside their m1m2 class overlap
all_species = c("elegans", "inopinata", "briggsae", "nigoni", "remanei", "bovis", "contortus", "tipulae", "pacificus")
CERP2 = matrix(rep(NA, 36), ncol = 9, dimnames = list(c("divergent", "tandem_m2m1", "any motif", "no_motif"), 
                                                      all_species))

for (species in c(1:length(all_species))){
  infile=read.table(paste("CELE2_CERP2_other_species/", all_species[species], ".summary", sep=""))
  CERP2[c(1:2),species] = infile[c(3,5),2]
  CERP2[3,species] = infile[6,2] - sum(infile[c(2:5),2])
  CERP2[4,species] = infile[1,2] - infile[6,2]
}

CELE2 = matrix(rep(NA, 36), ncol = 9, dimnames = list(c("divergent", "tandem_m2m1", "any_motif", "no_motif"), 
                                                      all_species))

for (species in c(1:length(all_species))){
  infile=read.table(paste("CELE2_CERP2_other_species/", all_species[species], ".summary", sep=""))
  CELE2[c(1:2),species] = infile[c(9,11),2]
  CELE2[3,species] = infile[12,2] - sum(infile[c(8:11),2])
  CELE2[4,species] = infile[7,2] - infile[12,2]
}

pdf(file = "plots/CERP2_CELE2_n_other_species.pdf", width = 8, height = 8, useDingbats = FALSE)
par(oma=c(2.2,0,0,0), mar=c(5,4,3,2), mfrow=c(1,2))
barplot(CERP2, beside=F, col=c("forestgreen", "deepskyblue", "grey30", "grey80"), las=2, main = "CERP2", ylab = "repeat number")
barplot(CELE2, beside=F, col=c("forestgreen", "deepskyblue", "grey30", "grey80"), las=2, main = "CELE2", ylab = "repeat number")
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(0, -0.9, legend = c("divergent", "tandem_m2m1", "any motif", "no motif"),
       fill = c("forestgreen", "deepskyblue", "grey30", "grey80"), bty="n", ncol = 4, xjust = 0.5)
dev.off()
