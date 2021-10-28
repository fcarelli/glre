library(DiffBind)
library(rtracklayer)

samplesID_all=list.files("relmapping/atac/tg_se.bwa_se_ce11.rm_unmapped.rm_contigs_ce11.rm_q10")
Peaks="relmapping/annot_ce/reg_elements_ce.bed"

samplesID=sapply(strsplit(samplesID_all, split = "_"), function(x) paste(x[3], x[4], x[6], sep="_"))
factors=rep("worm", length(samplesID_all))
condition=sapply(strsplit(samplesID_all, split = "_"), function(x) x[4])
treatment=rep("normal", length(samplesID_all))
replicate=rep(c(1,2), 2)
bamReads=paste("relmapping/atac/tg_se.bwa_se_ce11.rm_unmapped.rm_contigs_ce11.rm_q10", samplesID_all, sep="/")

samples=data.frame("SampleID"=samplesID, "Tissue"=samplesID, "Factor"=factors, "Condition"=condition, "Treatment"=treatment, "Replicate"=replicate, "bamReads"=bamReads, "Peaks"=rep(Peaks, 4), "PeakCaller"=rep("bed", 4))


diff_test <- dba(sampleSheet=samples)
diff_test <- dba.count(diff_test)
diff_test <- dba.contrast(diff_test, categories = DBA_CONDITION, minMembers=2)
diff_test <- dba.analyze(diff_test)

diff_test_ya_glp1 <- dba.report(diff_test, th=0.1)

diff_test_ya_glp1.DF=mcols(diff_test_ya_glp1)
k = which(diff_test_ya_glp1.DF$Fold < -2 & diff_test_ya_glp1.DF$FDR < 0.01)
diff_test_ya_glp1.DF.logF_filtered = diff_test_ya_glp1[k]


df <- data.frame(seqnames=seqnames(diff_test_ya_glp1.DF.logF_filtered),
                 starts=start(diff_test_ya_glp1.DF.logF_filtered)-1,
                 ends=end(diff_test_ya_glp1.DF.logF_filtered),
		 lfc=elementMetadata(diff_test_ya_glp1.DF.logF_filtered)$Fold,
		 fdr=elementMetadata(diff_test_ya_glp1.DF.logF_filtered)$FDR)
write.table(df, "RE_annotation/reg_elements_all.elegans.gl_specific.diffbind.bed", quote = F, sep = "\t", row.names = F, col.names = F)

