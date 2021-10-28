# dinucleotide frequency
library(grDevices)
repeat_class = c("inactive", "GLRE")

pdf(file = "plots/TT_periodicity_CERP2_elements.pdf", width = 12, height = 4, useDingbats = F)

par(mfrow=c(2,3), oma=c(0,0,0,0), mar=c(4,3,3,3))
for (rep_class in repeat_class){
  GLRE_CERP2_TT = read.table(paste("intact_repeats/elegans/elegans_CERP2.", rep_class, ".for_clustering.TT_freq", sep=""), row.names = 1)
  set.seed(12)
  kmeansObj = kmeans(GLRE_CERP2_TT, 2)
  
  seq_clustered_names = names(kmeansObj$cluster[kmeansObj$cluster == 2])
  seq_clustered = GLRE_CERP2_TT[row.names(GLRE_CERP2_TT) %in% seq_clustered_names,]
  
  seq_to_recluster_names = names(kmeansObj$cluster[kmeansObj$cluster == 1])
  seq_to_recluster = GLRE_CERP2_TT[row.names(GLRE_CERP2_TT) %in% seq_to_recluster_names,]
  
  reclustered_seq = seq_clustered
  all_ccf = c()
  
  for (i in c(1:dim(seq_to_recluster)[1])){
    ccfvalues = ccf(as.numeric(seq_to_recluster[i,]),as.numeric(kmeansObj$centers[2,]), lag.max = 5, plot = F)
    all_ccf = append(all_ccf, min(ccfvalues$acf))
    lowest_ccf = as.numeric(ccfvalues$lag[which(ccfvalues$acf == min(ccfvalues$acf))])
    if (lowest_ccf >= 0){
      shifted_seq = as.vector(c(rep(0, lowest_ccf), seq_to_recluster[i,][1:(length(seq_to_recluster[i,])-lowest_ccf)]))
    }
    else{
      shifted_seq = as.vector(c(seq_to_recluster[i,][(abs(lowest_ccf)+1):(length(seq_to_recluster[i,]))], rep(0, abs(lowest_ccf))))
    }
    names(shifted_seq) = names(seq_to_recluster[i,])
    reclustered_seq = rbind(reclustered_seq, shifted_seq)
  }
  
  reclustered_kmeansObj = kmeans(reclustered_seq, 2)
  
  image(t(GLRE_CERP2_TT)[, nrow(GLRE_CERP2_TT):1], yaxt = "n", main = "Original Data", ylab = paste(rep_class, "m1m2 pairs"), col=(c("white", "dodgerblue4")))
  image(t(reclustered_seq)[, nrow(reclustered_seq):1], yaxt = "n", main = "Clustered Data", col=(c("white", "dodgerblue4")))
  segments(0, kmeansObj$size[1]/dim(reclustered_seq)[1], 1, kmeansObj$size[1]/dim(reclustered_seq)[1])
  image(t(reclustered_seq)[, order(reclustered_kmeansObj$cluster)], yaxt = "n", main = "Re-clustered Data", col=(c("white", "dodgerblue4")))
  segments(0, reclustered_kmeansObj$size[1]/dim(reclustered_seq)[1], 1, reclustered_kmeansObj$size[1]/dim(reclustered_seq)[1])
  
}
dev.off()

pdf(file = "plots/TT_periodicity_CELE2_elements.pdf", width = 12, height = 4, useDingbats = F)

par(mfrow=c(2,3), oma=c(0,0,0,0), mar=c(4,3,3,3))
for (rep_class in repeat_class){
  GLRE_CERP2_TT = read.table(paste("intact_repeats/elegans/elegans_CELE2.", rep_class, ".for_clustering.TT_freq", sep=""), row.names = 1)
  set.seed(12)
  kmeansObj = kmeans(GLRE_CERP2_TT, 2)
  
  seq_clustered_names = names(kmeansObj$cluster[kmeansObj$cluster == 2])
  seq_clustered = GLRE_CERP2_TT[row.names(GLRE_CERP2_TT) %in% seq_clustered_names,]
  
  seq_to_recluster_names = names(kmeansObj$cluster[kmeansObj$cluster == 1])
  seq_to_recluster = GLRE_CERP2_TT[row.names(GLRE_CERP2_TT) %in% seq_to_recluster_names,]
  
  reclustered_seq = seq_clustered
  all_ccf = c()
  
  for (i in c(1:dim(seq_to_recluster)[1])){
    ccfvalues = ccf(as.numeric(seq_to_recluster[i,]),as.numeric(kmeansObj$centers[2,]), lag.max = 5, plot = F)
    all_ccf = append(all_ccf, min(ccfvalues$acf))
    lowest_ccf = as.numeric(ccfvalues$lag[which(ccfvalues$acf == min(ccfvalues$acf))])
    if (lowest_ccf >= 0){
      shifted_seq = as.vector(c(rep(0, lowest_ccf), seq_to_recluster[i,][1:(length(seq_to_recluster[i,])-lowest_ccf)]))
    }
    else{
      shifted_seq = as.vector(c(seq_to_recluster[i,][(abs(lowest_ccf)+1):(length(seq_to_recluster[i,]))], rep(0, abs(lowest_ccf))))
    }
    names(shifted_seq) = names(seq_to_recluster[i,])
    reclustered_seq = rbind(reclustered_seq, shifted_seq)
  }
  
  reclustered_kmeansObj = kmeans(reclustered_seq, 2)
  
  image(t(GLRE_CERP2_TT)[, nrow(GLRE_CERP2_TT):1], yaxt = "n", main = "Original Data", ylab = paste(rep_class, "m1m2 pairs"), col=(c("white", "dodgerblue4")))
  image(t(reclustered_seq)[, nrow(reclustered_seq):1], yaxt = "n", main = "Clustered Data", col=(c("white", "dodgerblue4")))
  segments(0, kmeansObj$size[1]/dim(reclustered_seq)[1], 1, kmeansObj$size[1]/dim(reclustered_seq)[1])
  image(t(reclustered_seq)[, order(reclustered_kmeansObj$cluster)], yaxt = "n", main = "Re-clustered Data", col=(c("white", "dodgerblue4")))
  segments(0, reclustered_kmeansObj$size[1]/dim(reclustered_seq)[1], 1, reclustered_kmeansObj$size[1]/dim(reclustered_seq)[1])
  
}
dev.off()

repeat_class = c("intact", "inactive", "GLRE")
pdf(file = "plots/TT_periodicity_divergent_m1m2_elements.pdf", width = 12, height = 6, useDingbats = F)

par(mfrow=c(3,3), oma=c(0,0,0,0), mar=c(4,3,3,3))
for (rep_class in repeat_class){
  GLRE_CERP2_TT = read.table(paste("intact_repeats/elegans/elegans_divergent_m1m2.", rep_class, ".for_clustering.TT_freq", sep=""), row.names = 1)
  set.seed(12)
  kmeansObj = kmeans(GLRE_CERP2_TT, 2)
  
  seq_clustered_names = names(kmeansObj$cluster[kmeansObj$cluster == 2])
  seq_clustered = GLRE_CERP2_TT[row.names(GLRE_CERP2_TT) %in% seq_clustered_names,]
  
  seq_to_recluster_names = names(kmeansObj$cluster[kmeansObj$cluster == 1])
  seq_to_recluster = GLRE_CERP2_TT[row.names(GLRE_CERP2_TT) %in% seq_to_recluster_names,]
  
  reclustered_seq = seq_clustered
  all_ccf = c()
  
  for (i in c(1:dim(seq_to_recluster)[1])){
    ccfvalues = ccf(as.numeric(seq_to_recluster[i,]),as.numeric(kmeansObj$centers[2,]), lag.max = 5, plot = F)
    all_ccf = append(all_ccf, min(ccfvalues$acf))
    lowest_ccf = as.numeric(ccfvalues$lag[which(ccfvalues$acf == min(ccfvalues$acf))])
    if (lowest_ccf >= 0){
      shifted_seq = as.vector(c(rep(0, lowest_ccf), seq_to_recluster[i,][1:(length(seq_to_recluster[i,])-lowest_ccf)]))
    }
    else{
      shifted_seq = as.vector(c(seq_to_recluster[i,][(abs(lowest_ccf)+1):(length(seq_to_recluster[i,]))], rep(0, abs(lowest_ccf))))
    }
    names(shifted_seq) = names(seq_to_recluster[i,])
    reclustered_seq = rbind(reclustered_seq, shifted_seq)
  }
  
  reclustered_kmeansObj = kmeans(reclustered_seq, 2)
  
  image(t(GLRE_CERP2_TT)[, nrow(GLRE_CERP2_TT):1], yaxt = "n", main = "Original Data", ylab = paste(rep_class, "m1m2 pairs"), col=(c("white", "dodgerblue4")))
  image(t(reclustered_seq)[, nrow(reclustered_seq):1], yaxt = "n", main = "Clustered Data", col=(c("white", "dodgerblue4")))
  segments(0, kmeansObj$size[1]/dim(reclustered_seq)[1], 1, kmeansObj$size[1]/dim(reclustered_seq)[1])
  image(t(reclustered_seq)[, order(reclustered_kmeansObj$cluster)], yaxt = "n", main = "Re-clustered Data", col=(c("white", "dodgerblue4")))
  segments(0, reclustered_kmeansObj$size[1]/dim(reclustered_seq)[1], 1, reclustered_kmeansObj$size[1]/dim(reclustered_seq)[1])
  
}
dev.off()

pdf(file = "plots/TT_periodicity_tandem_m2m1_elements.pdf", width = 12, height = 6, useDingbats = F)

par(mfrow=c(3,3), oma=c(0,0,0,0), mar=c(4,3,3,3))
for (rep_class in repeat_class){
  GLRE_CERP2_TT = read.table(paste("intact_repeats/elegans/elegans_tandem_m2m1.", rep_class, ".for_clustering.TT_freq", sep=""), row.names = 1)
  set.seed(12)
  kmeansObj = kmeans(GLRE_CERP2_TT, 2)
  
  seq_clustered_names = names(kmeansObj$cluster[kmeansObj$cluster == 2])
  seq_clustered = GLRE_CERP2_TT[row.names(GLRE_CERP2_TT) %in% seq_clustered_names,]
  
  seq_to_recluster_names = names(kmeansObj$cluster[kmeansObj$cluster == 1])
  seq_to_recluster = GLRE_CERP2_TT[row.names(GLRE_CERP2_TT) %in% seq_to_recluster_names,]
  
  reclustered_seq = seq_clustered
  all_ccf = c()
  
  for (i in c(1:dim(seq_to_recluster)[1])){
    ccfvalues = ccf(as.numeric(seq_to_recluster[i,]),as.numeric(kmeansObj$centers[2,]), lag.max = 5, plot = F)
    all_ccf = append(all_ccf, min(ccfvalues$acf))
    lowest_ccf = as.numeric(ccfvalues$lag[which(ccfvalues$acf == min(ccfvalues$acf))])
    if (lowest_ccf >= 0){
      shifted_seq = as.vector(c(rep(0, lowest_ccf), seq_to_recluster[i,][1:(length(seq_to_recluster[i,])-lowest_ccf)]))
    }
    else{
      shifted_seq = as.vector(c(seq_to_recluster[i,][(abs(lowest_ccf)+1):(length(seq_to_recluster[i,]))], rep(0, abs(lowest_ccf))))
    }
    names(shifted_seq) = names(seq_to_recluster[i,])
    reclustered_seq = rbind(reclustered_seq, shifted_seq)
  }
  
  reclustered_kmeansObj = kmeans(reclustered_seq, 2)
  
  image(t(GLRE_CERP2_TT)[, nrow(GLRE_CERP2_TT):1], yaxt = "n", main = "Original Data", ylab = paste(rep_class, "m1m2 pairs"), col=(c("white", "dodgerblue4")))
  image(t(reclustered_seq)[, nrow(reclustered_seq):1], yaxt = "n", main = "Clustered Data", col=(c("white", "dodgerblue4")))
  segments(0, kmeansObj$size[1]/dim(reclustered_seq)[1], 1, kmeansObj$size[1]/dim(reclustered_seq)[1])
  image(t(reclustered_seq)[, order(reclustered_kmeansObj$cluster)], yaxt = "n", main = "Re-clustered Data", col=(c("white", "dodgerblue4")))
  segments(0, reclustered_kmeansObj$size[1]/dim(reclustered_seq)[1], 1, reclustered_kmeansObj$size[1]/dim(reclustered_seq)[1])
  
}
dev.off()

#BiocManager::install("periodicDNA")
#library(ggplot2)
#library(magrittr)
#library(periodicDNA)

#pdf(file = "plots/TT_periodicity_CERP2_inactive.periodicDNA.pdf", width = 8, height = 6, useDingbats = F)
#inactive_m1m2 = import("intact_repeats/elegans/elegans_CERP2.inactive.for_clustering.fa", format = "fasta")

#inactive_periodicity_result <- getPeriodicity(
#  inactive_m1m2,
#  motif = 'TT', 
#  BPPARAM = setUpBPPARAM(6)
#)
#plotPeriodicityResults(inactive_periodicity_result)
#dev.off()

#pdf(file = "plots/TT_periodicity_CERP2_GLRE.periodicDNA.pdf", width = 8, height = 6, useDingbats = F)
#GLRE_m1m2 = import("intact_repeats/elegans/elegans_CERP2.GLRE.for_clustering.fa", format = "fasta")

#GLRE_periodicity_result <- getPeriodicity(
#  GLRE_m1m2,
#  motif = 'TT', 
#  BPPARAM = setUpBPPARAM(6)
#)
#plotPeriodicityResults(GLRE_periodicity_result)
#dev.off()

#intact_m1m2 = import("intact_repeats/elegans/elegans_CERP2.intact.for_clustering.fa", format = "fasta")
#pdf(file = "plots/TT_periodicity_CERP2_intact.periodicDNA.pdf", width = 8, height = 6, useDingbats = F)

#intact_periodicity_result <- getPeriodicity(
#  intact_m1m2,
#  motif = 'TT', 
#  BPPARAM = setUpBPPARAM(6)
#)
#plotPeriodicityResults(intact_periodicity_result)

#dev.off()
