# dinucleotide frequency
library(ggplot2)
library(magrittr)
library(periodicDNA)

motif_types = c("GLpromoter", "inactive", "intact")
repeat_types = c("CERP2", "CELE2")
for (motif_type in motif_types){
    for (repeat_type in repeat_types){
        m1m2_type = import(paste("intact_repeats/elegans/elegans_", repeat_type, ".", motif_type, ".for_clustering.fa", sep=""), format = "fasta")
        m1m2_type_periodicity_result <- getPeriodicity(
            m1m2_type,
            motif = 'TT', 
            BPPARAM = setUpBPPARAM(6), range_spectrum = seq(1, 90)
        )
        par(mar=c(4,4,3,3), oma=c(0,0,0,0), mfrow=c(1,1))
        pdf(file = paste("plots/TT_periodicity_", repeat_type, "_", motif_type, ".periodicDNA.pdf", sep=""), width = 10, height = 6, useDingbats = F)
        print(plotPeriodicityResults(m1m2_type_periodicity_result))
        dev.off()
    }
}
