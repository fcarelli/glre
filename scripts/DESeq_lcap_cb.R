library("DESeq2")
library("GenomicAlignments")
library("rtracklayer")
library("GenomicFeatures")

dir.create("gene_annotation")

#Import the gene model and the bam files that I want the counts to be done on 
sample_ID_list=data.frame("sample" = c("annot_cb_lcap_wt_ya_rep1", "annot_cb_lcap_wt_ya_rep2", "annot_cb_lcap_glp1_ya_rep1", "annot_cb_lcap_glp1_ya_rep2"), "strain" = c("wt", "wt", "glp1", "glp1"))

strains=levels(sample_ID_list$strain)

bam_files=c()
for (i in c(1:dim(sample_ID_list)[1])){
  bam_files = c(bam_files, paste("relmapping/lcap/trim20.bwa_pe_cb4.rm_unmapped_pe.rm_contigs_cb4.rm_q10/", sample_ID_list[i,1], ".trim20.bwa_pe_cb4.rm_unmapped_pe.rm_contigs_cb4.rm_q10.bam", sep=""))
}

bf = BamFileList(bam_files )
gtf = makeTxDbFromGFF("species/briggsae/gene_annotation/briggsae.annotations.coding_genes.gtf", format="gtf")
exonsByGene <- exonsBy( gtf, by="gene" )
genehits <- summarizeOverlaps(exonsByGene, bf, mode="Union", singleEnd=F, ignore.strand=F, preprocess.reads=invertStrand, fragments=T)
colData(genehits)$strain = factor(sample_ID_list$strain)
de <- DESeqDataSet(genehits, ~strain)
de1 <- DESeq(de)

de_est <- estimateSizeFactors(de)
normalized_counts <- counts(de_est, normalized=TRUE)
row.names(normalized_counts) = substr(row.names(normalized_counts),6,nchar(row.names(normalized_counts)))
write.table(normalized_counts, file="gene_annotation/briggsae.normalized_counts_DESeq.txt", sep="\t", quote=F, col.names=NA)

all_strains=unique(sample_ID_list$strain)
mutant="glp1"
res <- results(de1, c('strain', mutant, 'wt'))
row.names(res) = substr(row.names(res),6,nchar(row.names(res)))
downreg=res[res$log2FoldChange < -2 & res$padj < 0.001 & !is.na(res$padj),]
upreg=res[res$log2FoldChange > 2 & res$padj < 0.001 & !is.na(res$padj),]

write.table(downreg, file = paste("gene_annotation/briggsae.genes_downregulated_", mutant, "_vs_wt.txt", sep=""), quote=F, col.names = F, row.names = T, sep="\t")
write.table(upreg, file = paste("gene_annotation/briggsae.genes_upregulated_", mutant, "_vs_wt.txt", sep=""), quote=F, col.names = F, row.names = T, sep="\t")
write.table(res, file = paste("gene_annotation/briggsae.genes_DESeq_", mutant, "_vs_wt.txt", sep=""), quote=F, row.names = T, sep="\t")

