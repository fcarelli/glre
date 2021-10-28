library("DESeq2")
library("GenomicAlignments")
library("rtracklayer")
library("GenomicFeatures")
library(vidger)

dir.create("DE_analysis")

#Import the gene model and the bam files that I want the counts to be done on 
samples_list = lapply(strsplit(list.files("data/elegans/rnaseq/fastq/"), "[.]"), '[[', 1)
sample_ID_list=data.frame("sample" = unlist(samples_list), "strain" = unlist(lapply(strsplit(unlist(samples_list), "_"), '[[', 2)))

strains=levels(sample_ID_list$strain)

bam_files=c()
for (i in c(1:dim(sample_ID_list)[1])){
  bam_files = c(bam_files, paste("data/elegans/rnaseq/alignment/", sample_ID_list[i,1], ".Aligned.sortedByCoord.out.bam", sep=""))
}

bf = BamFileList(bam_files )
#gtf = import("species/elegans/gene_annotation/elegans.annotations.coding_genes.gtf")
#gtf = reduce(split(gtf, gtf$gene_id))
gtf = makeTxDbFromGFF("species/elegans/gene_annotation/elegans.annotations.coding_genes.gtf", format="gtf")
exonsByGene <- exonsBy( gtf, by="gene" )
#genehits <- summarizeOverlaps(gtf, bf, mode="Union", ignore.strand=F, preprocess.reads=invertStrand)
genehits <- summarizeOverlaps(exonsByGene, bf, mode="Union", singleEnd=T, ignore.strand=F, preprocess.reads=invertStrand)
colData(genehits)$strain = factor(sample_ID_list$strain)
de <- DESeqDataSet(genehits, ~strain)
de1 <- DESeq(de)

de_est <- estimateSizeFactors(de)
normalized_counts <- counts(de_est, normalized=TRUE)
row.names(normalized_counts) = substr(row.names(normalized_counts),6,nchar(row.names(normalized_counts)))
write.table(normalized_counts, file="DE_analysis/normalized_counts_DESeq.txt", sep="\t", quote=F, col.names=NA)

all_strains=unique(sample_ID_list$strain)
all_mutants=all_strains[all_strains != "wt"]
for (mutant in all_mutants){
  res <- results(de1, c('strain', mutant, 'wt'))
  row.names(res) = substr(row.names(res),6,nchar(row.names(res)))
  downreg=res[res$log2FoldChange < 0 & res$padj < 0.001 & !is.na(res$padj),]
  upreg=res[res$log2FoldChange > 0 & res$padj < 0.001 & !is.na(res$padj),]
  write.table(downreg, file = paste("DE_analysis/genes_downregulated_", mutant, "_vs_N2.txt", sep=""), quote=F, col.names = F, row.names = T, sep="\t")
  write.table(upreg, file = paste("DE_analysis/genes_upregulated_", mutant, "_vs_N2.txt", sep=""), quote=F, col.names = F, row.names = T, sep="\t")
  write.table(res, file = paste("DE_analysis/genes_DESeq_", mutant, "_vs_N2.txt", sep=""), quote=F, row.names = T, sep="\t")
}

# volcanoplot

pdf(file = "plots/him17_vs_wt_volcano.pdf", width = 4, height = 2)
vsVolcano(
  x = 'him17', y = 'wt',
  data = de1, d.factor = 'strain', type = 'deseq',
  padj = 0.001, x.lim = c(-4, 4), lfc = NULL, title = TRUE,
  legend = TRUE, grid = TRUE, data.return = FALSE,
  xaxis.text.size = 10, yaxis.text.size = 10, xaxis.title.size = 10, yaxis.title.size = 10,
  main.title.size = 10, legend.text.size = 10, highlight = NULL
)
dev.off()
