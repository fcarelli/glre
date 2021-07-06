rule modern_enrichment:
  input:
    "data/external_data/annotatedPeak.bed",
    "RE_annotation/reg_elements_all.elegans.gl_specific.bed",
    "RE_annotation/reg_elements_all.elegans.not_gl_specific.bed",
    "motif_enrichment/elegans/elegans.m1m2_clusters.bed",
  output:
    "modern_analysis/annotatedPeak.total_peaks_TF.txt",
    "modern_analysis/annotatedPeak.peaks_on_GLRE.txt",
    "modern_analysis/annotatedPeak.peaks_on_nonGLRE.txt",
    "modern_analysis/annotatedPeak.peaks_GLRE.enrichment",
    "modern_analysis/annotatedPeak.peaks_on_GLRE.m1m2.txt",
    "modern_analysis/annotatedPeak.peaks_on_GLRE.no_m1m2.txt",
    "modern_analysis/annotatedPeak.peaks_GLRE.m1m2.enrichment",
  shell:
    '''
    grep young_adult {input[0]} | cut -f 4 | uniq -c | awk '{{print $2 "\t" $1}}' | sort -k 1,1 > {output[0]}
    grep young_adult {input[0]} | sed 's/^...//' | cut -f 1,2,3,4 | intersectBed -a stdin -b {input[1]} -u | cut -f 4 | uniq -c | awk '{{print $2 "\t" $1}}' | sort -k 1,1 | join - {output[0]} | awk 'BEGIN{{OFS="\t";}}{{print $1, $2, $3, $2/$3}}' > {output[1]}
    grep young_adult {input[0]} | sed 's/^...//' | cut -f 1,2,3,4 | intersectBed -a stdin -b {input[2]} -u | cut -f 4 | uniq -c | awk '{{print $2 "\t" $1}}' | sort -k 1,1 | join - {output[0]} | awk 'BEGIN{{OFS="\t";}}{{print $1, $2, $3, $2/$3}}' > {output[2]}
    nonGLRE_n=$(wc -l < {input[2]})
    GLRE_n=$(wc -l < {input[1]})
    join {output[1]} {output[2]} | awk -v glre="$GLRE_n" -v non_glre="$nonGLRE_n" 'BEGIN{{OFS="\t";}}{{print $1, $3, $2, $5, $4, $7, ($2/glre)/($5/non_glre)}}' | sort -k7,7rn > {output[3]}
    grep young_adult {input[0]} | sed 's/^...//' | cut -f 1,2,3,4 | intersectBed -a stdin -b {input[1]} -u | intersectBed -a stdin -b {input[3]} -u | cut -f 4 | uniq -c | awk '{{print $2 "\t" $1}}' | sort -k 1,1 | join - {output[1]} | awk 'BEGIN{{OFS="\t";}}{{print $1, $2, $3, $2/$3}}' > {output[4]}
    grep young_adult {input[0]} | sed 's/^...//' | cut -f 1,2,3,4 | intersectBed -a stdin -b {input[1]} -u | intersectBed -a stdin -b {input[3]} -v | cut -f 4 | uniq -c | awk '{{print $2 "\t" $1}}' | sort -k 1,1 | join - {output[1]} | awk 'BEGIN{{OFS="\t";}}{{print $1, $2, $3, $2/$3}}' > {output[5]}
    m1m2_GLRE_n=$(intersectBed -a {input[1]} -b {input[3]} -u | wc -l)
    non_m1m2_GLRE_n=$(intersectBed -a {input[1]} -b {input[3]} -v | wc -l)
    join {output[4]} {output[5]} | awk -v m1m2="$m1m2_GLRE_n" -v non_m1m2="$non_m1m2_GLRE_n" 'BEGIN{{OFS="\t";}}{{print $1, $3, $2, $5, $4, $7, ($2/m1m2)/($5/non_m1m2)}}' | sort -k7,7rn > {output[6]}
    '''



rule samples_ce_r1:
  input:
    lambda wildcards: config["chip_samples"][wildcards.sample]
  output:
    'data/elegans/chipseq/fastq/{sample}.r1.fq.gz'
  shell:
    '''
    cp {input} {output}
    '''

rule trim_adapt:
  input:
    expand('data/elegans/chipseq/fastq/{sample}.r1.fq.gz', sample=config["chip_samples"])
  output:
    "data/elegans/chipseq/fastq_trimmed/{sample}.r1.fq.gz_trimming_report.txt",
    "data/elegans/chipseq/fastq_trimmed/{sample}.r1_trimmed.fq.gz"
  params:
    "data/elegans/chipseq/fastq_trimmed"
  resources:
    cpus=4
  shell:
    '''
    trim_galore --gzip --trim-n -j {resources.cpus} -o {params} {input}
    '''


rule bwa:
  input:
    "bwa_idx/elegans.amb",
    "bwa_idx/elegans.ann",
    "bwa_idx/elegans.bwt",
    "bwa_idx/elegans.pac",
    "bwa_idx/elegans.sa",
    "species/elegans/genome/elegans.fa",
    "data/elegans/chipseq/fastq_trimmed/{sample}.r1_trimmed.fq.gz"
  output:
    "data/elegans/chipseq/aligned/{sample}.bam"
  resources:
    cpus=8
  params:
    'bwa_idx/elegans'
  shell:
    '''
    bwa mem -t {resources.cpus} {params} {input[6]} | samtools view -@ {resources.cpus} -bT {input[5]} - > {output}.temp
    samtools sort -@ {resources.cpus} {output}.temp > {output}
    rm {output}.temp
    samtools index {output}
    '''

rule macs2:
  input:
    "data/elegans/chipseq/aligned/{sample}.bam",
    "species/elegans/genome/elegans.chrom.sizes.txt"
  output:
    "data/elegans/chipseq/peaks/{sample}/{sample}_peaks.narrowPeak",
    "data/elegans/chipseq/peaks/{sample}/{sample}_model.r",
    "data/elegans/chipseq/peaks/{sample}/{sample}_peaks.xls",
    "data/elegans/chipseq/peaks/{sample}/{sample}_summits.bed",
    "data/elegans/chipseq/peaks/{sample}/{sample}_treat_pileup.bdg",
    temp("data/elegans/chipseq/peaks/{sample}/{sample}_treat_pileup.sorted.bdg"),
    "data/elegans/chipseq/peaks/{sample}/{sample}_treat_pileup.sorted.bw"
  params:
    "data/elegans/chipseq/peaks/{sample}"
  shell:
    '''
    macs2 callpeak --bdg --SPMR --gsize ce --nolambda -n {wildcards.sample} --outdir {params} -t {input[0]}
    sort -k 1,1 -k2,2n {output[4]} > {output[5]}
    bedGraphToBigWig {output[5]} {input[1]} {output[6]}
    '''

rule yapc:
  input:
    "data/elegans/chipseq/peaks/{sample}_rep1/{sample}_rep1_treat_pileup.sorted.bw",
    "data/elegans/chipseq/peaks/{sample}_rep2/{sample}_rep2_treat_pileup.sorted.bw",
  output:
    "data/elegans/chipseq/yapc/{sample}.smooth_100_yapc_0.001.bed",
    "data/elegans/chipseq/yapc/{sample}.smooth_100_yapc_0.005.bed",
    "data/elegans/chipseq/yapc/{sample}.smooth_100_yapc_0.01.bed",
    "data/elegans/chipseq/yapc/{sample}.smooth_100_yapc_0.05.bed",
    "data/elegans/chipseq/yapc/{sample}.smooth_100_yapc_coverage.bw",
    "data/elegans/chipseq/yapc/{sample}.smooth_100_yapc_d2smooth.bw",
    "data/elegans/chipseq/yapc/{sample}.smooth_100_yapc_peaksall.bed",
    "data/elegans/chipseq/yapc/{sample}.smooth_100_yapc_peaksall.tsv",
    "data/elegans/chipseq/yapc/{sample}.smooth_100_yapc.tsv"
  params:
    "data/elegans/chipseq/yapc/{sample}.smooth_100_yapc"
  shell:
    '''
    yapc --smoothing-window-width 100 {params} {wildcards.sample} {input[0]} {input[1]}
    '''


rule chip_stats:
  input:
    "data/elegans/chipseq/yapc/elegans_HIM17.smooth_100_yapc_0.001.bed",
    "data/elegans/chipseq/yapc/elegans_XND1.smooth_100_yapc_0.001.bed",
    "species/elegans/repeats_dfam/elegans_dfam.repeats.id.bed",
    "motif_enrichment/elegans/elegans.m1m2_clusters.bed",
  output:
    "chip_stats/chip_stats.txt",
    "plots/HIM17_XND1_peak_overlap.pdf"
  shell:
    '''
    echo -e "factor\tpeaks\toverlap_with_other_factor\toverlap_with_CERP2\tfraction_of_CERP2\toverlap_with_CELE2\tfraction_of_CELE2\toverlap_with_m1m2\tfraction_of_m1m2" > {output}
    HIM17_peak_n=$(wc -l < {input[0]})
    XND1_peak_n=$(wc -l < {input[1]})
    peak_overlap=$(intersectBed -a {input[0]} -b {input[1]} -u	| wc -l)
    CERP2_n=$(grep CERP2 {input[2]} | wc -l)
    HIM17_CERP2_overlap=$(grep CERP2 {input[2]} | intersectBed -a stdin -b {input[0]} -u | wc -l)
    HIM17_CERP2_fraction=$(awk -v ovlp="$HIM17_CERP2_overlap" -v rep="$CERP2_n" 'BEGIN{{print ovlp/rep}}')
    XND1_CERP2_overlap=$(grep CERP2 {input[2]} | intersectBed -a stdin -b {input[1]} -u | wc -l)
    XND1_CERP2_fraction=$(awk -v ovlp="$XND1_CERP2_overlap" -v rep="$CERP2_n" 'BEGIN{{print ovlp/rep}}')
    CELE2_n=$(grep CELE2 {input[2]} | wc -l)
    HIM17_CELE2_overlap=$(grep CELE2 {input[2]} | intersectBed -a stdin -b {input[0]} -u | wc -l)
    HIM17_CELE2_fraction=$(awk -v ovlp="$HIM17_CELE2_overlap" -v rep="$CELE2_n" 'BEGIN{{print ovlp/rep}}')
    XND1_CELE2_overlap=$(grep CELE2 {input[2]} | intersectBed -a stdin -b {input[1]} -u | wc -l)
    XND1_CELE2_fraction=$(awk -v ovlp="$XND1_CELE2_overlap" -v rep="$CELE2_n" 'BEGIN{{print ovlp/rep}}')
    m1m2_n=$(wc -l < {input[3]})
    HIM17_m1m2_overlap=$(intersectBed -a {input[3]} -b {input[0]} -u | wc -l)
    HIM17_m1m2_fraction=$(awk -v ovlp="$HIM17_m1m2_overlap" -v rep="$m1m2_n" 'BEGIN{{print ovlp/rep}}')
    XND1_m1m2_overlap=$(intersectBed -a {input[3]} -b {input[1]} -u | wc -l)
    XND1_m1m2_fraction=$(awk -v ovlp="$XND1_m1m2_overlap" -v rep="$m1m2_n" 'BEGIN{{print ovlp/rep}}')
    echo -e "HIM-17_peak_n\t"$HIM17_peak_n"\t"$peak_overlap"\t"$HIM17_CERP2_overlap"\t"$HIM17_CERP2_fraction"\t"$HIM17_CELE2_overlap"\t"$HIM17_CELE2_fraction"\t"$HIM17_m1m2_overlap"\t"$HIM17_m1m2_fraction >> {output}
    echo -e "XND1-17_peak_n\t"$XND1_peak_n"\t"$peak_overlap"\t"$XND1_CERP2_overlap"\t"$XND1_CERP2_fraction"\t"$XND1_CELE2_overlap"\t"$XND1_CELE2_fraction"\t"$XND1_m1m2_overlap"\t"$XND1_m1m2_fraction >> {output}
    Rscript scripts/HIM17_XND1_peak_overlap.R
    '''


rule chip_corr:
  input:
    "data/elegans/chipseq/yapc/elegans_HIM17.smooth_100_yapc_coverage.bw",
    "data/elegans/chipseq/yapc/elegans_XND1.smooth_100_yapc_coverage.bw",
    "data/elegans/chipseq/yapc/elegans_HIM17.smooth_100_yapc_0.001.bed",
    "data/elegans/chipseq/yapc/elegans_XND1.smooth_100_yapc_0.001.bed",
  output:
    "data/elegans/chipseq/correlation/HIM17_XND1_loci.bed",
    "data/elegans/chipseq/correlation/corr_HIM17_XND1",
    "data/elegans/chipseq/correlation/corr_HIM17_XND1.spearman.scatter.pdf",
  resources:
    cpus=8
  shell:
    '''
    cat {input[2]} {input[3]} | cut -f 1,2,3 | sort -k 1,1 -k2,2n | grep -v "\#" | mergeBed -i stdin > {output[0]}
    multiBigwigSummary BED-file --bwfiles {input[0]} {input[1]} --outFileName {output[1]} --BED {output[0]} --labels HIM17 XND1 --numberOfProcessors {resources.cpus}
    plotCorrelation --corData {output[1]} --corMethod spearman --whatToPlot scatterplot --plotFile {output[2]} --plotTitle str1.spearman --plotFileFormat pdf --plotHeight 6 --plotWidth 6
    '''


rule chip_motif_stats:
  input:
    "data/elegans/chipseq/yapc/{sample}.smooth_100_yapc_0.001.bed",
    "motif_enrichment/elegans/elegans.m1m2_clusters.bed",
    "RE_annotation/reg_elements_all.elegans.bed"
  output:
    "chip_stats/{sample}.motif_stats.txt",
  shell:
    '''
    echo -e "factor\tmotif_RE\tmotif_only\tRE_only\tno_overlap" > {output[0]}
    motif_RE=$(intersectBed -a {input[0]} -b {input[1]} -u | intersectBed -a stdin -b {input[2]} -u | wc -l)
    motif_noRE=$(intersectBed -a {input[0]} -b {input[1]} -u | intersectBed -a stdin -b {input[2]} -v | wc -l)
    nomotif_RE=$(intersectBed -a {input[0]} -b {input[1]} -v | intersectBed -a stdin -b {input[2]} -u | wc -l)
    nomotif_noRE=$(intersectBed -a {input[0]} -b {input[1]} -v | intersectBed -a stdin -b {input[2]} -v | wc -l)
    echo -e {wildcards.sample}"\t"$motif_RE"\t"$motif_noRE"\t"$nomotif_RE"\t"$nomotif_noRE >> {output}
    '''

rule chip_motif_enrichment:
  input:
    "species/elegans/genome/elegans.fa",
    "data/elegans/chipseq/yapc/{sample}.smooth_100_yapc_0.001.bed",
    "motif_enrichment/elegans/elegans.m1.bed",
    "motif_enrichment/elegans/elegans.m2.bed",
  output:
    temp("data/elegans/chipseq/yapc/{sample}.smooth_100_yapc_0.001.fa"),
    directory("meme/{sample}.enrichment"),
    temp("data/elegans/chipseq/yapc/{sample}.smooth_100_yapc_0.001.nomotif.fa"),
    directory("meme/{sample}.nomotif.enrichment"),
  resources:
    ntasks=8
  shell:
    '''
    fastaFromBed -fi {input[0]} -bed {input[1]} -fo {output[0]}
    meme -oc {output[1]} -objfun de -dna -revcomp -nmotifs 5 -p {resources.ntasks} {output[0]}
    intersectBed -a {input[1]} -b {input[2]} {input[3]} -v | fastaFromBed -fi {input[0]} -bed stdin -fo {output[2]}
    meme -oc {output[3]} -objfun de -dna -revcomp -nmotifs 5 -p {resources.ntasks} {output[2]}
    '''

rule star_idx:
  input:
    "species/elegans/genome/elegans.fa",
    "species/elegans/gene_annotation/elegans.annotations.genes.gtf"
  output:
    directory("star_idx/elegans")
  resources:
    cpus=8
  shell:
    '''
    STAR --genomeSAindexNbases 12 --sjdbGTFfile {input[1]} --runMode genomeGenerate --runThreadN {resources.cpus} --genomeDir {output} --genomeFastaFiles {input[0]}
    '''

rule samples_rnaseq:
  input:
    lambda wildcards: config["rnaseq_samples"][wildcards.sample]
  output:
    'data/elegans/rnaseq/fastq/{sample}.r1.fq.gz'
  shell:
    '''
    cp {input} {output}
    '''

rule star:
  input:
    "star_idx/elegans",
    "data/elegans/rnaseq/fastq/{sample}.r1.fq.gz",
    "species/elegans/genome/elegans.chrom.sizes.txt"
  output:
    "data/elegans/rnaseq/alignment/{sample}.Aligned.sortedByCoord.out.bam",
    "data/elegans/rnaseq/alignment/{sample}.Signal.UniqueMultiple.str1.out.wig",
    "data/elegans/rnaseq/alignment/{sample}.Signal.UniqueMultiple.str2.out.wig"
  params:
    "data/elegans/rnaseq/alignment/{sample}."
  resources:
    cpus=8
  shell:
    '''
    STAR --readFilesCommand zcat --runThreadN {resources.cpus} --genomeDir {input[0]} --readFilesIn {input[1]} --outFileNamePrefix {params} --outSAMtype BAM SortedByCoordinate --outSAMattrIHstart 0 --outWigType wiggle --twopassMode Basic
    '''


rule bigwig_generator:
  input:
    "data/elegans/rnaseq/alignment/{sample}.Signal.UniqueMultiple.str1.out.wig",
    "data/elegans/rnaseq/alignment/{sample}.Signal.UniqueMultiple.str2.out.wig",
    "species/elegans/genome/elegans.chrom.sizes.txt"
  output:
    "data/elegans/rnaseq/alignment/{sample}.Signal.UniqueMultiple.str1.out.bw",
    "data/elegans/rnaseq/alignment/{sample}.Signal.UniqueMultiple.str2.out.bw"
  shell:
    '''
    wigToBigWig {input[0]} {input[2]} {output[0]}
    wigToBigWig {input[1]} {input[2]} {output[1]}
    '''


rule merged_rnaseq_tracks:
  input:
    "data/elegans/rnaseq/alignment/elegans_{sample}_rep1.Signal.UniqueMultiple.str{strand}.out.bw",
    "data/elegans/rnaseq/alignment/elegans_{sample}_rep2.Signal.UniqueMultiple.str{strand}.out.bw",
  output:
    "data/elegans/rnaseq/tracks/elegans_{sample}_merged.Signal.UniqueMultiple.str{strand}.out.bw",
  resources:
    cpus=6
  shell:
    '''
    bigwigCompare --bigwig1 {input[0]} --bigwig2 {input[1]} --operation mean -p {resources.cpus} -o {output[0]} -of bigwig -bs 1
    '''


rule kallisto_index:
  input:
    "species/elegans/genome/elegans.fa",
    "species/elegans/gene_annotation/elegans.gene_annotation.coding_gene.bed",
  output:
    "species/elegans/gene_annotation/elegans.gene_annotation.coding_gene.fa",
    "kallisto_index/elegans"
  shell:
    '''
    fastaFromBed -fi {input[0]} -bed {input[1]} -nameOnly -split -s | awk 'BEGIN{{FS="(";}}{{print $1}}' > {output[0]}
    kallisto index -i {output[1]} {output[0]}
    '''


rule kallisto:
  input:
    "kallisto_index/elegans",
    "data/elegans/rnaseq/fastq/{sample}.r1.fq.gz",
    "species/elegans/gene_annotation/elegans.gene_annotation.gene_transcript_id.txt",
  output:
    "data/elegans/rnaseq/kallisto/{sample}/abundance.h5",
    "data/elegans/rnaseq/kallisto/{sample}/abundance.tsv",
    "data/elegans/rnaseq/kallisto/{sample}/run_info.json",
    "data/elegans/rnaseq/kallisto/{sample}/abundance.genes.txt",
  resources:
    cpus=10
  params:
    "data/elegans/rnaseq/kallisto/{sample}"
  shell:
    '''
    kallisto quant -i {input[0]} -b 100 --single -o {params} -l 200 -s 20 -t {resources.cpus} --rf-stranded {input[1]}
    python scripts/sum_TPM_per_gene.py {input[2]} {output[1]} {output[3]}
    sort -k 1,1 {output[3]} > {output[3]}.sorted
    mv {output[3]}.sorted {output[3]}
    '''

rule kallisto_final_table:
  input:
    kallisto_out=expand("data/elegans/rnaseq/kallisto/{sample}/abundance.genes.txt", sample=[item for item in config["rnaseq_samples"] if "adult" not in item]),
  output:
    'data/elegans/rnaseq/gene_expression.kallisto.txt',
  params:
    filenames=[item for item in config["rnaseq_samples"] if "adult" not in item],
  shell:
    '''
    echo {params.filenames} | tr " " "\t" > {output[0]}
    paste {input.kallisto_out} | awk '{{ for (i=3;i<=NF;i+=2) $i="" }} 1' | tr " " "\t" >> {output[0]}
    '''

#rule sleuth_analysis:
#  input:
#    kallisto_out=expand('data/elegans/rnaseq/kallisto/{sample}/abundance.h5', sample=config["rnaseq_samples"]),
#    kallisto_table='data/elegans/rnaseq/gene_expression.kallisto.txt',
#  output:
#    'sleuth_DE/shared_DE_genes.LRT_0.01.txt',
#    'sleuth_DE/him17_specific_DE_genes.LRT_0.01.txt',
#    'sleuth_DE/xnd1_specific_DE_genes.LRT_0.01.txt',
#    'sleuth_DE/shared_DE_genes.him17_xnd1_down.LRT_0.01.txt',
#    'sleuth_DE/shared_DE_genes.him17_xnd1_up.LRT_0.01.txt',
#    'sleuth_DE/him17_specific_DE_genes.him17_down.LRT_0.01.txt',
#    'sleuth_DE/him17_specific_DE_genes.him17_up.LRT_0.01.txt',
#    'sleuth_DE/xnd1_specific_DE_genes.xnd1_down.LRT_0.01.txt',
#    'sleuth_DE/xnd1_specific_DE_genes.xnd1_up.LRT_0.01.txt',
#  shell:
#    '''
#    Rscript scripts/sleuth_mutants.R
#    sort -k 1,1 {output[0]} | join - {input.kallisto_table} | awk '{{if(($8 + $9)/2 < ($6 + $7)/2 && ($10 + $11)/2 < ($6 + $7)/2) print $1}}' > {output[3]}
#    sort -k 1,1 {output[0]} | join - {input.kallisto_table} | awk '{{if(($8 + $9)/2 > ($6 + $7)/2 && ($10 + $11)/2 > ($6 + $7)/2) print $1}}' > {output[4]}
#    sort -k 1,1 {output[1]} | join - {input.kallisto_table} | awk '{{if(($8 + $9)/2 < ($6 + $7)/2) print $1}}' > {output[5]}
#    sort -k 1,1 {output[1]} | join - {input.kallisto_table} | awk '{{if(($8 + $9)/2 > ($6 + $7)/2) print $1}}' > {output[6]}
#    sort -k 1,1 {output[2]} | join - {input.kallisto_table} | awk '{{if(($10 + $11)/2 < ($6 + $7)/2) print $1}}' > {output[7]}
#    sort -k 1,1 {output[2]} | join - {input.kallisto_table} | awk '{{if(($10 + $11)/2 > ($6 + $7)/2) print $1}}' > {output[8]}
#    '''


rule deseq: 
  input:
    expand("data/elegans/rnaseq/alignment/{sample}.Aligned.sortedByCoord.out.bam", sample=config["rnaseq_samples"]),
    "species/elegans/gene_annotation/elegans.annotations.coding_genes.gtf"
  output:
    "DE_analysis/genes_DESeq_him17_vs_N2.txt",
    "DE_analysis/genes_downregulated_him17_vs_N2.txt",
    "DE_analysis/genes_upregulated_him17_vs_N2.txt",
    "DE_analysis/genes_DESeq_xnd1_vs_N2.txt",
    "DE_analysis/genes_downregulated_xnd1_vs_N2.txt",
    "DE_analysis/genes_upregulated_xnd1_vs_N2.txt",
    "DE_analysis/normalized_counts_DESeq.txt"
  shell:
    '''
    Rscript scripts/DESeq_mutants.R
    '''


rule chip_targets:
  input:
    "promoter_annotation/promoters_all.elegans.bed",
    "data/elegans/chipseq/yapc/{sample}.smooth_100_yapc_0.001.bed",
  output:
    "DE_analysis/{sample}.bound_genes.txt"
  shell:
    '''
    intersectBed -a {input[0]} -b {input[1]} -u | cut -f 4 | tr ',' '\n' | sort | uniq | sed '1d' > {output[0]}
    '''

rule direct_targets:
  input:
    "DE_analysis/genes_downregulated_him17_vs_N2.txt",
    "DE_analysis/genes_downregulated_xnd1_vs_N2.txt",
    "DE_analysis/elegans_HIM17.bound_genes.txt",
    "DE_analysis/elegans_XND1.bound_genes.txt"
  output:
    "DE_analysis/genes_downregulated_him17_vs_N2.direct.txt",
    "DE_analysis/genes_downregulated_him17_vs_N2.not_direct.txt",
    "DE_analysis/genes_downregulated_xnd1_vs_N2.direct.txt",
    "DE_analysis/genes_downregulated_xnd1_vs_N2.not_direct.txt"
  shell:
    '''
    join {input[0]} {input[2]} | awk '{{print $1}}' > {output[0]}
    join -v 1 {input[0]} {input[2]} | awk '{{print $1}}' > {output[1]}
    join {input[1]} {input[3]} | awk '{{print $1}}' > {output[2]}
    join -v 1 {input[1]} {input[3]} | awk '{{print $1}}' > {output[3]}
    '''

#rule xnd1_other_repeats_bound:
#  input:
#    'heatmaps_int_files/HIM17_XND1_CELE2_CERP2_m1m2_RE.heatmap.bed',
#    'RE_annotation/reg_elements_all.elegans.bed',
#    'species/elegans/repeats_dfam/elegans_dfam.repeats.id.bed',
#  output:
#    temp('XND1_other_repeats/XND1_other_repeats.bound'),
#    'XND1_other_repeats/XND1_other_repeats.summary',
#    'plots/XND1_binding_other_repeats.pdf'
#  shell:
#    '''
#    sed '1d' {input[0]} | awk 'BEGIN{{OFS="\t";}}{{if ($13 == "cluster_3") print $1, $2, $3}}' | sort -k 1,1 -k2,2n | intersectBed -a stdin -b {input[1]} -v | intersectBed -b stdin -a {input[2]} -u | cut -f 4 | sort | uniq -c | awk '{{print $2 "\t" $1}}' > {output[0]}
#    cut -f 4 {input[2]} | sort | uniq -c | awk '{{print $2 "\t" $1}}' | join - {output[0]} | awk 'BEGIN{{OFS="\t";}}{{if ($3 > 50) print $1, $2, $3, $3/$2}}' | sort -k4,4rn > {output[1]}
#    Rscript scripts/XND1_other_repeats.R
#    '''

rule xnd1_other_repeats_bound:
  input:
    'heatmaps_int_files/XND1_only.bed',
    'RE_annotation/reg_elements_all.elegans.bed',
    'species/elegans/repeats_dfam/elegans_dfam.repeats.id.bed',
  output:
    temp('XND1_other_repeats/XND1_other_repeats.bound'),
    'XND1_other_repeats/XND1_other_repeats.summary',
    'plots/XND1_binding_other_repeats.pdf'
  shell:
    '''
    intersectBed -a {input[0]} -b {input[1]} -v | intersectBed -b stdin -a {input[2]} -u | cut -f 4 | sort | uniq -c | awk '{{print $2 "\t" $1}}' > {output[0]}
    cut -f 4 {input[2]} | sort | uniq -c | awk '{{print $2 "\t" $1}}' | join - {output[0]} | awk 'BEGIN{{OFS="\t";}}{{if ($3 > 50) print $1, $2, $3, $3/$2}}' | sort -k4,4rn > {output[1]}
    Rscript scripts/XND1_other_repeats.R
    '''

rule him17_xnd1_shared_statistics:
  input:
    'data/elegans/chipseq/yapc/elegans_HIM17.smooth_100_yapc_0.001.bed',
    'data/elegans/chipseq/yapc/elegans_XND1.smooth_100_yapc_0.001.bed',
    'RE_annotation/reg_elements_all.elegans.gl_specific.bed',
    'RE_annotation/reg_elements_all.elegans.not_gl_specific.bed',
    'motif_enrichment/elegans/elegans.m1m2_clusters.bed',
  output:
    'shared_peaks/shared_peaks.bed',
    'shared_peaks/statistics.txt',
    'shared_peaks/m1m2_overlap.txt'
    'plots/shared_peaks_promoters_overlap.pdf',
    'plots/shared_peaks_m1m2_overlap.pdf',
  shell:
    '''
    intersectBed -a {input[0]} -b {input[1]} > {output[0]}
    echo -e "GL_promoters_m1m2\tGL_promoters_not_m1m2\tnotGL_promoters_m1m2\tnot_GL_promoters_not_m1m2" > {output[1]}
    GL_m1m2=$(grep coding_promoter {input[2]} | intersectBed -a stdin -b {input[4]} -u | wc -l)
    GL_no_m1m2=$(grep coding_promoter {input[2]} | intersectBed -a stdin -b {input[4]} -v | wc -l)
    no_GL_m1m2=$(grep coding_promoter {input[3]} | intersectBed -a stdin -b {input[4]} -u | wc -l)
    no_GL_no_m1m2=$(grep coding_promoter {input[3]} | intersectBed -a stdin -b {input[4]} -v | wc -l)
    echo -e "total\t"$GL_m1m2"\t"$GL_no_m1m2"\t"$no_GL_m1m2"\t"$no_GL_no_m1m2 >> {output[1]}
    GL_m1m2=$(grep coding_promoter {input[2]} | intersectBed -a {output[0]} -b stdin -u | intersectBed -a stdin -b {input[4]} -u | wc -l)
    GL_no_m1m2=$(grep coding_promoter {input[2]} | intersectBed -a {output[0]} -b stdin -u | intersectBed -a stdin -b {input[4]} -v | wc -l)
    no_GL_m1m2=$(grep coding_promoter {input[3]} | intersectBed -a {output[0]} -b stdin -u | intersectBed -a stdin -b {input[4]} -u | wc -l)
    no_GL_no_m1m2=$(grep coding_promoter {input[3]} | intersectBed -a {output[0]} -b stdin -u | intersectBed -a stdin -b {input[4]} -v | wc -l)
    echo -e "shared_peaks\t"$GL_m1m2"\t"$GL_no_m1m2"\t"$no_GL_m1m2"\t"$no_GL_no_m1m2 >> {output[1]}
    echo "loci" > {output[2]}
    m1m2_n=$(wc -l < {input[4]})
    shared_p_n=$(wc -l {output[0]})
    overlap_m1m2_shared=$(intersectBed -a {input[4]} -b {output{0}} -u | wc -l)
    echo -e "m1m2_pairs\t"$m1m2_n >> {output[2]}
    echo -e "shared_peaks\t"$shared_p_n >> {output[2]}
    echo -e "overlap\t"$overlap_m1m2_shared >> {output[2]}
    Rscript scripts/shared_peaks_stats.R
    '''

