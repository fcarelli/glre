# plotting statistics on m1m2 arrangements and intermotif distances for different species
source("scripts/common_R_functions.R")

pdf(file = "plots/m1m2_motif_arrangements.all_species.pdf", width = 10, height = 10)

all_species_full_names = c("C. elegans", "C. inopinata", "C. briggsae", "C. nigoni", "C. remanei", "C. becei", "C. quiockensis", "C. bovis", "C. monodelphis", "H. contortus", "A. ceylanicum", "O. tipulae", "H. bacteriophora", "P. pacificus", "P. redivivus")
all_species = c("elegans", "inopinata", "briggsae", "nigoni", "remanei", "becei", "quiockensis", "bovis", "monodelphis", "contortus", "ceylanicum", "tipulae", "bacteriophora", "pacificus", "redivivus")
par(mfrow=c(3,5), oma=c(0,0,0,0), mar=c(4,4,3,3))

for (species in c(1:length(all_species))){
  m1m2_all=read.table(paste("motif_enrichment/", all_species[species], "/", all_species[species], ".m1m2_clusters.arrangements_all.summary", sep=""), header=T)
  chr_size=read.table(paste("species/", all_species[species], "/genome/", all_species[species], ".chrom.sizes.txt", sep=""))
  genome_size=sum(chr_size$V2)
  k=barplot(m1m2_all$N_clusters/genome_size*1000000, las=1, names.arg=m1m2_all$motif_arrangement,
            col=c("grey90", "grey60", "grey30", "grey1"), main = all_species_full_names[species], ylab = "m1m2 pairs/Mbp", ylim=c(0, 50))
}
dev.off()


# density plot for Caenorhabditis only - divergent_m1m2

all_species_full_names = c("C. elegans", "C. inopinata", "C. briggsae", "C. nigoni", "C. remanei", "C. becei", "C. quiockensis", "C. bovis", "C. monodelphis")
all_species = c("elegans", "inopinata", "briggsae", "nigoni", "remanei", "becei", "quiockensis", "bovis", "monodelphis")

m1m2_dist_ele=read.table("motif_enrichment/elegans/elegans.m1m2_clusters.arrangements_by_distance.summary", header=T)

m1m2_dist_all = data.frame(row.names = m1m2_dist_ele$motif_distance, elegans = m1m2_dist_ele$divergent_cluster)
for (species in c(2:length(all_species))){
  species_name=all_species[species]
  m1m2_dist=read.table(paste("motif_enrichment/", all_species[species], "/", all_species[species], ".m1m2_clusters.arrangements_by_distance.summary", sep=""), header=T)
  m1m2_dist_all = cbind(m1m2_dist_all, m1m2_dist$divergent_cluster)
}
names(m1m2_dist_all) = all_species

pdf(file = "plots/m1m2_motif_distance_divergent.caenorhabditis.pdf", width = 8, height = 6)

par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(4,4,3,3))

k = c()
dist_bins = c(10:30)
for (i in c(1:length(dist_bins))){
  k = append(k, rep(dist_bins[i], m1m2_dist_all[i, 1]))
}

plot(density(k, bw=0.5), lwd = 2, col = "deepskyblue4", ylim = c(0, 0.4), main = "divergent_m1m2 Caenorhabditis", xlab = "m1-m2 distance (bp)", las = 1)

for (species in c(2: length(all_species))){
  k = c()
  dist_bins = c(10:30)
  for (i in c(1:length(dist_bins))){
    k = append(k, rep(dist_bins[i], m1m2_dist_all[i, species]))
  }
  lines(density(k, bw=0.5), lwd = 2, col = t_col("grey", 70), ylim = c(0, 0.4))
}

k = c()
dist_bins = c(10:30)
for (i in c(1:length(dist_bins))){
  k = append(k, rep(dist_bins[i], m1m2_dist_all[i, 6]))
}

lines(density(k, bw=0.5), lwd = 2, lty=2, col = "forestgreen", ylim = c(0, 0.4))

k = c()
dist_bins = c(10:30)
for (i in c(1:length(dist_bins))){
  k = append(k, rep(dist_bins[i], m1m2_dist_all[i, length(all_species)]))
}

lines(density(k, bw=0.5), lwd = 2, lty=2, col = "deepskyblue", ylim = c(0, 0.4))

dev.off()


# density plot for C. elegans, A. ceylanicum and P. pacificus - divergent_m1m2

all_species_full_names = c("C. elegans", "A. ceylanicum", "P. pacificus")
all_species = c("elegans", "ceylanicum", "pacificus")

m1m2_dist_ele=read.table("motif_enrichment/elegans/elegans.m1m2_clusters.arrangements_by_distance.summary", header=T)

m1m2_dist_all = data.frame(row.names = m1m2_dist_ele$motif_distance, elegans = m1m2_dist_ele$divergent_cluster)
for (species in c(2:length(all_species))){
  species_name=all_species[species]
  m1m2_dist=read.table(paste("motif_enrichment/", all_species[species], "/", all_species[species], ".m1m2_clusters.arrangements_by_distance.summary", sep=""), header=T)
  m1m2_dist_all = cbind(m1m2_dist_all, m1m2_dist$divergent_cluster)
}
names(m1m2_dist_all) = all_species

pdf(file = "plots/m1m2_motif_distance_divergent.other_species.pdf", width = 8, height = 6)

par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(4,4,3,3))

k = c()
dist_bins = c(10:30)
for (i in c(1:length(dist_bins))){
  k = append(k, rep(dist_bins[i], m1m2_dist_all[i, 1]))
}

plot(density(k, bw=0.5), lwd = 2, col = "deepskyblue4", ylim = c(0, 0.55), main = "divergent_m1m2 nematodes", xlab = "m1-m2 distance (bp)", las = 1)


k = c()
dist_bins = c(10:30)
for (i in c(1:length(dist_bins))){
  k = append(k, rep(dist_bins[i], m1m2_dist_all[i, 2]))
}

lines(density(k, bw=0.5), lwd = 2, lty=2, col = "orangered", ylim = c(0, 0.4))

k = c()
dist_bins = c(10:30)
for (i in c(1:length(dist_bins))){
  k = append(k, rep(dist_bins[i], m1m2_dist_all[i, 3]))
}

lines(density(k, bw=0.5), lwd = 2, lty=2, col = "black", ylim = c(0, 0.4))
dev.off()


# density plot for Caenorhabditis only - tandem_m2m1

all_species_full_names = c("C. elegans", "C. inopinata", "C. briggsae", "C. nigoni", "C. remanei", "C. becei", "C. quiockensis", "C. bovis", "C. monodelphis")
all_species = c("elegans", "inopinata", "briggsae", "nigoni", "remanei", "becei", "quiockensis", "bovis", "monodelphis")

m1m2_dist_ele=read.table("motif_enrichment/elegans/elegans.m1m2_clusters.arrangements_by_distance.summary", header=T)

m1m2_dist_all = data.frame(row.names = m1m2_dist_ele$motif_distance, elegans = m1m2_dist_ele$tandem_m2m1_cluster)
for (species in c(2:length(all_species))){
  species_name=all_species[species]
  m1m2_dist=read.table(paste("motif_enrichment/", all_species[species], "/", all_species[species], ".m1m2_clusters.arrangements_by_distance.summary", sep=""), header=T)
  m1m2_dist_all = cbind(m1m2_dist_all, m1m2_dist$tandem_m2m1_cluster)
}
names(m1m2_dist_all) = all_species

pdf(file = "plots/m1m2_motif_distance_tandem_m2m1.caenorhabditis.pdf", width = 8, height = 6)

par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(4,4,3,3))

k = c()
dist_bins = c(10:30)
for (i in c(1:length(dist_bins))){
  k = append(k, rep(dist_bins[i], m1m2_dist_all[i, 1]))
}

plot(density(k, bw=0.5), lwd = 2, col = "deepskyblue4", ylim = c(0, 0.3), main = "tandem_m2m1 Caenorhabditis", xlab = "m1-m2 distance (bp)", las = 1)

for (species in c(2: length(all_species))){
  k = c()
  dist_bins = c(10:30)
  for (i in c(1:length(dist_bins))){
    k = append(k, rep(dist_bins[i], m1m2_dist_all[i, species]))
  }
  lines(density(k, bw=0.5), lwd = 2, col = t_col("grey", 70), ylim = c(0, 0.4))
}

k = c()
dist_bins = c(10:30)
for (i in c(1:length(dist_bins))){
  k = append(k, rep(dist_bins[i], m1m2_dist_all[i, 6]))
}

lines(density(k, bw=0.5), lwd = 2, lty=2, col = "forestgreen", ylim = c(0, 0.4))

k = c()
dist_bins = c(10:30)
for (i in c(1:length(dist_bins))){
  k = append(k, rep(dist_bins[i], m1m2_dist_all[i, length(all_species)]))
}

lines(density(k, bw=0.5), lwd = 2, lty=2, col = "deepskyblue", ylim = c(0, 0.4))

dev.off()



