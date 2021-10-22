rule modencode_him17_chip:
  input:
    'bwa_idx/elegans.amb',
    'bwa_idx/elegans.ann',
    'bwa_idx/elegans.bwt',
    'bwa_idx/elegans.pac',
    'bwa_idx/elegans.sa',
    'species/elegans/genome/elegans.fa',
  output:
    'data/external_data/chipseq/fastq/him17_rep1.fq.gz',
    'data/external_data/chipseq/alignment/him17_rep1.bam',
    'data/external_data/chipseq/peaks/him17_rep1/him17_rep1_peaks.narrowPeak',
    'data/external_data/chipseq/peaks/him17_rep1/him17_rep1_model.r',
    'data/external_data/chipseq/peaks/him17_rep1/him17_rep1_peaks.xls',
    'data/external_data/chipseq/peaks/him17_rep1/him17_rep1_summits.bed',
    'data/external_data/chipseq/peaks/him17_rep1/him17_rep1_treat_pileup.bdg',
  resources:
    cpus=8
  params:
    'bwa_idx/elegans',
    'data/external_data/chipseq/peaks/'
  shell:
    '''
    wget ftp://data.modencode.org/all_files/cele-raw-4/3916_SDQ0801_HIM17_FEM2_AD_r1_export.fq.gz -O {output[0]}
    bwa mem -t {resources.cpus} {params[0]} {output[0]} | samtools view -@ {resources.cpus} -q 10 -bT {input[5]} - > {output[1]}.temp
    samtools sort -@ {resources.cpus} {output[1]}.temp > {output[1]}
    rm {output[1]}.temp
    mkdir {params[1]}him17_rep1
    macs2 callpeak --bdg --SPMR --gsize ce -n him17_rep1 --outdir {params[1]}him17_rep1 -t {output[1]}
    '''


rule modern_him17_peak_merge:
  input:
    'data/external_data/annotatedPeak.bed',
    'data/external_data/chipseq/peaks/him17_rep1/him17_rep1_peaks.narrowPeak',
  output:
    'data/external_data/annotatedPeak_him17.bed',
  shell:
    '''
    grep young_adult {input[0]} | sed 's/^...//' | cut -f 1,2,3,4 > {output[0]}
    awk 'BEGIN{{OFS="\t";}}{{print $1, $2, $3, "him17_young_adult"}}' {input[1]} >> {output[0]}
    '''


rule modern_enrichment:
  input:
    'data/external_data/annotatedPeak_him17.bed',
    'RE_annotation/reg_elements_all.elegans.gl_specific.bed',
    'RE_annotation/reg_elements_all.elegans.not_gl_specific.bed',
    'motif_enrichment/elegans/elegans.m1m2_clusters.bed',
  output:
    'modern_analysis/annotatedPeak.total_peaks_TF.txt',
    'modern_analysis/annotatedPeak.peaks_on_GLRE.txt',
    'modern_analysis/annotatedPeak.peaks_on_nonGLRE.txt',
    'modern_analysis/annotatedPeak.peaks_GLRE.enrichment',
    'modern_analysis/annotatedPeak.peaks_on_GLRE.m1m2.txt',
    'modern_analysis/annotatedPeak.peaks_on_GLRE.no_m1m2.txt',
    'modern_analysis/annotatedPeak.peaks_GLRE.m1m2.enrichment',
    temp('modern_analysis/gl_coding_promoters.bed'),
    'modern_analysis/annotatedPeak.peaks_on_GLpromoters.txt',
    'modern_analysis/annotatedPeak.peaks_on_GLpromoters.m1m2.txt',
    'modern_analysis/annotatedPeak.peaks_on_GLpromoters.no_m1m2.txt',
    'modern_analysis/annotatedPeak.peaks_GLpromoters.m1m2.enrichment',
  shell:
    '''
    cut -f 4 {input[0]} | uniq -c | awk '{{print $2 "\t" $1}}' | sort -k 1,1 > {output[0]}
    intersectBed -a {input[0]} -b {input[1]} -u | cut -f 4 | uniq -c | awk '{{print $2 "\t" $1}}' | sort -k 1,1 | join - {output[0]} | awk 'BEGIN{{OFS="\t";}}{{print $1, $2, $3, $2/$3}}' > {output[1]}
    intersectBed -a {input[0]} -b {input[2]} -u | cut -f 4 | uniq -c | awk '{{print $2 "\t" $1}}' | sort -k 1,1 | join - {output[0]} | awk 'BEGIN{{OFS="\t";}}{{print $1, $2, $3, $2/$3}}' > {output[2]}
    nonGLRE_n=$(wc -l < {input[2]})
    GLRE_n=$(wc -l < {input[1]})
    join {output[1]} {output[2]} | awk -v glre="$GLRE_n" -v non_glre="$nonGLRE_n" 'BEGIN{{OFS="\t";}}{{print $1, $3, $2, $5, $4, $7, ($2/glre)/($5/non_glre)}}' | sort -k7,7rn > {output[3]}
    intersectBed -a {input[0]} -b {input[1]} -u | intersectBed -a stdin -b {input[3]} -u | cut -f 4 | uniq -c | awk '{{print $2 "\t" $1}}' | sort -k 1,1 | join - {output[1]} | awk 'BEGIN{{OFS="\t";}}{{print $1, $2, $3, $2/$3}}' > {output[4]}
    intersectBed -a {input[0]} -b {input[1]} -u | intersectBed -a stdin -b {input[3]} -v | cut -f 4 | uniq -c | awk '{{print $2 "\t" $1}}' | sort -k 1,1 | join - {output[1]} | awk 'BEGIN{{OFS="\t";}}{{print $1, $2, $3, $2/$3}}' > {output[5]}
    m1m2_GLRE_n=$(intersectBed -a {input[1]} -b {input[3]} -u | wc -l)
    non_m1m2_GLRE_n=$(intersectBed -a {input[1]} -b {input[3]} -v | wc -l)
    join {output[4]} {output[5]} | awk -v m1m2="$m1m2_GLRE_n" -v non_m1m2="$non_m1m2_GLRE_n" 'BEGIN{{OFS="\t";}}{{print $1, $3, $2, $5, $4, $7, ($2/m1m2)/($5/non_m1m2)}}' | sort -k7,7rn > {output[6]}
    grep coding_promoter {input[1]} > {output[7]}
    intersectBed -a {input[0]} -b {output[7]} -u | cut -f 4 | uniq -c | awk '{{print $2 "\t" $1}}' | sort -k 1,1 | join - {output[0]} | awk 'BEGIN{{OFS="\t";}}{{print $1, $2, $3, $2/$3}}' > {output[8]}
    intersectBed -a {input[0]} -b {output[7]} -u | intersectBed -a stdin -b {input[3]} -u | cut -f 4 | uniq -c | awk '{{print $2 "\t" $1}}' | sort -k 1,1 | join - {output[8]} | awk 'BEGIN{{OFS="\t";}}{{print $1, $2, $3, $2/$3}}' > {output[9]}
    intersectBed -a {input[0]} -b {output[7]} -u | intersectBed -a stdin -b {input[3]} -v | cut -f 4 | uniq -c | awk '{{print $2 "\t" $1}}' | sort -k 1,1 | join - {output[8]} | awk 'BEGIN{{OFS="\t";}}{{print $1, $2, $3, $2/$3}}' > {output[10]}
    m1m2_GL_prom_n=$(intersectBed -a {output[7]} -b {input[3]} -u | wc -l)
    non_m1m2_GL_prom_n=$(intersectBed -a {output[7]} -b {input[3]} -v | wc -l)
    join {output[9]} {output[10]} | awk -v m1m2="$m1m2_GL_prom_n" -v non_m1m2="$non_m1m2_GL_prom_n" 'BEGIN{{OFS="\t";}}{{print $1, $3, $2, $5, $4, $7, ($2/m1m2)/($5/non_m1m2)}}' | sort -k7,7rn > {output[11]}
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
    'data/elegans/chipseq/fastq_trimmed/{sample}.r1.fq.gz_trimming_report.txt',
    'data/elegans/chipseq/fastq_trimmed/{sample}.r1_trimmed.fq.gz'
  params:
    'data/elegans/chipseq/fastq_trimmed'
  resources:
    cpus=4
  shell:
    '''
    trim_galore --gzip --trim-n -j {resources.cpus} -o {params} {input}
    '''
rule bwa:
  input:
    'bwa_idx/elegans.amb',
    'bwa_idx/elegans.ann',
    'bwa_idx/elegans.bwt',
    'bwa_idx/elegans.pac',
    'bwa_idx/elegans.sa',
    'species/elegans/genome/elegans.fa',
    'data/elegans/chipseq/fastq_trimmed/{sample}.r1_trimmed.fq.gz'
  output:
    'data/elegans/chipseq/aligned/{sample}.mapq10.bam'
  resources:
    cpus=8
  params:
    'bwa_idx/elegans'
  shell:
    '''
    bwa mem -t {resources.cpus} {params} {input[6]} | samtools view -@ {resources.cpus} -q 10 -bT {input[5]} - > {output}.temp
    samtools sort -@ {resources.cpus} {output}.temp > {output}
    rm {output}.temp
    samtools index {output}
    '''


rule macs2:
  input:
    'data/elegans/chipseq/aligned/{sample}.mapq10.bam',
    'species/elegans/genome/elegans.chrom.sizes.txt'
  output:
    'data/elegans/chipseq/peaks/{sample}/{sample}_peaks.narrowPeak',
    'data/elegans/chipseq/peaks/{sample}/{sample}_model.r',
    'data/elegans/chipseq/peaks/{sample}/{sample}_peaks.xls',
    'data/elegans/chipseq/peaks/{sample}/{sample}_summits.bed',
    'data/elegans/chipseq/peaks/{sample}/{sample}_treat_pileup.bdg',
    temp('data/elegans/chipseq/peaks/{sample}/{sample}_treat_pileup.sorted.bdg'),
    'data/elegans/chipseq/peaks/{sample}/{sample}_treat_pileup.sorted.bw'
  params:
    'data/elegans/chipseq/peaks/{sample}'
  shell:
    '''
    macs2 callpeak --bdg --SPMR --gsize ce --nolambda -n {wildcards.sample} --outdir {params} -t {input[0]}
    sort -k 1,1 -k2,2n {output[4]} > {output[5]}
    bedGraphToBigWig {output[5]} {input[1]} {output[6]}
    '''


rule yapc:
  input:
    'data/elegans/chipseq/peaks/{sample}_rep1/{sample}_rep1_treat_pileup.sorted.bw',
    'data/elegans/chipseq/peaks/{sample}_rep2/{sample}_rep2_treat_pileup.sorted.bw',
  output:
    'data/elegans/chipseq/yapc/{sample}.smooth_100_yapc_0.001.bed',
    'data/elegans/chipseq/yapc/{sample}.smooth_100_yapc_0.005.bed',
    'data/elegans/chipseq/yapc/{sample}.smooth_100_yapc_0.01.bed',
    'data/elegans/chipseq/yapc/{sample}.smooth_100_yapc_0.05.bed',
    'data/elegans/chipseq/yapc/{sample}.smooth_100_yapc_coverage.bw',
    'data/elegans/chipseq/yapc/{sample}.smooth_100_yapc_d2smooth.bw',
    'data/elegans/chipseq/yapc/{sample}.smooth_100_yapc_peaksall.bed',
    'data/elegans/chipseq/yapc/{sample}.smooth_100_yapc_peaksall.tsv',
    'data/elegans/chipseq/yapc/{sample}.smooth_100_yapc.tsv',
    'data/elegans/chipseq/yapc/{sample}.smooth_100_yapc_0.00001.bed'
  params:
    'data/elegans/chipseq/yapc/{sample}.smooth_100_yapc'
  shell:
    '''
    yapc --smoothing-window-width 100 {params} {wildcards.sample} {input[0]} {input[1]}
    sample_name=$(head -n 1 {output[8]} | cut -f 10)
    sed '1d' {output[8]} | awk -v var=$sample_name 'BEGIN{{OFS="\t";}}{{if ($10 >= 5) print $1, $2, $3, "Name=;" var "=" $10, 0, ".", $2, $3, "#0072b2"}}' > {output[9]}
    '''


rule chip_stats:
  input:
    'data/elegans/chipseq/yapc/elegans_HIM17.smooth_100_yapc_0.00001.bed',
    'species/elegans/repeats_dfam/elegans_dfam.repeats.id.bed',
    'motif_enrichment/elegans/elegans.m1m2_clusters.bed',
  output:
    'chip_stats/chip_stats.txt',
  shell:
    '''
    echo -e "factor\tpeaks\toverlap_with_CERP2\tfraction_of_CERP2\toverlap_with_CELE2\tfraction_of_CELE2\toverlap_with_m1m2\tfraction_of_m1m2" > {output}
    HIM17_peak_n=$(wc -l < {input[0]})
    CERP2_n=$(grep CERP2 {input[1]} | wc -l)
    HIM17_CERP2_overlap=$(grep CERP2 {input[1]} | intersectBed -a stdin -b {input[0]} -u | wc -l)
    HIM17_CERP2_fraction=$(awk -v ovlp="$HIM17_CERP2_overlap" -v rep="$CERP2_n" 'BEGIN{{print ovlp/rep}}')
    CELE2_n=$(grep CELE2 {input[1]} | wc -l)
    HIM17_CELE2_overlap=$(grep CELE2 {input[1]} | intersectBed -a stdin -b {input[0]} -u | wc -l)
    HIM17_CELE2_fraction=$(awk -v ovlp="$HIM17_CELE2_overlap" -v rep="$CELE2_n" 'BEGIN{{print ovlp/rep}}')
    m1m2_n=$(wc -l < {input[3]})
    HIM17_m1m2_overlap=$(intersectBed -a {input[2]} -b {input[0]} -u | wc -l)
    HIM17_m1m2_fraction=$(awk -v ovlp="$HIM17_m1m2_overlap" -v rep="$m1m2_n" 'BEGIN{{print ovlp/rep}}')
    echo -e "HIM-17_peak_n\t"$HIM17_peak_n"\t"$HIM17_CERP2_overlap"\t"$HIM17_CERP2_fraction"\t"$HIM17_CELE2_overlap"\t"$HIM17_CELE2_fraction"\t"$HIM17_m1m2_overlap"\t"$HIM17_m1m2_fraction >> {output}
    '''


rule chip_motif_enrichment:
  input:
    'species/elegans/genome/elegans.fa',
    'data/elegans/chipseq/yapc/{sample}.smooth_100_yapc_0.00001.bed',
    'motif_enrichment/elegans/elegans.m1.bed',
    'motif_enrichment/elegans/elegans.m2.bed',
  output:
    temp('data/elegans/chipseq/yapc/{sample}.smooth_100_yapc_0.001.fa'),
    directory('meme/{sample}.enrichment'),
    temp('data/elegans/chipseq/yapc/{sample}.smooth_100_yapc_0.001.nomotif.fa'),
    directory('meme/{sample}.nomotif.enrichment'),
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
    'species/elegans/genome/elegans.fa',
    'species/elegans/gene_annotation/elegans.annotations.genes.gtf'
  output:
    directory('star_idx/elegans')
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
    'star_idx/elegans',
    'data/elegans/rnaseq/fastq/{sample}.r1.fq.gz',
    'species/elegans/genome/elegans.chrom.sizes.txt'
  output:
    'data/elegans/rnaseq/alignment/{sample}.Aligned.sortedByCoord.out.bam',
    'data/elegans/rnaseq/alignment/{sample}.Signal.UniqueMultiple.str1.out.wig',
    'data/elegans/rnaseq/alignment/{sample}.Signal.UniqueMultiple.str2.out.wig'
  params:
    'data/elegans/rnaseq/alignment/{sample}.'
  resources:
    cpus=8
  shell:
    '''
    STAR --readFilesCommand zcat --runThreadN {resources.cpus} --genomeDir {input[0]} --readFilesIn {input[1]} --outFileNamePrefix {params} --outSAMtype BAM SortedByCoordinate --outSAMattrIHstart 0 --outWigType wiggle --twopassMode Basic
    '''


rule bigwig_generator:
  input:
    'data/elegans/rnaseq/alignment/{sample}.Signal.UniqueMultiple.str1.out.wig',
    'data/elegans/rnaseq/alignment/{sample}.Signal.UniqueMultiple.str2.out.wig',
    'species/elegans/genome/elegans.chrom.sizes.txt'
  output:
    'data/elegans/rnaseq/alignment/{sample}.Signal.UniqueMultiple.str1.out.bw',
    'data/elegans/rnaseq/alignment/{sample}.Signal.UniqueMultiple.str2.out.bw'
  shell:
    '''
    wigToBigWig {input[0]} {input[2]} {output[0]}
    wigToBigWig {input[1]} {input[2]} {output[1]}
    '''


rule merged_rnaseq_tracks:
  input:
    'data/elegans/rnaseq/alignment/elegans_{sample}_rep1.Signal.UniqueMultiple.str{strand}.out.bw',
    'data/elegans/rnaseq/alignment/elegans_{sample}_rep2.Signal.UniqueMultiple.str{strand}.out.bw',
  output:
    'data/elegans/rnaseq/tracks/elegans_{sample}_merged.Signal.UniqueMultiple.str{strand}.out.bw',
  resources:
    cpus=6
  shell:
    '''
    bigwigCompare --bigwig1 {input[0]} --bigwig2 {input[1]} --operation mean -p {resources.cpus} -o {output[0]} -of bigwig -bs 1
    '''


rule kallisto_index:
  input:
    'species/elegans/genome/elegans.fa',
    'species/elegans/gene_annotation/elegans.gene_annotation.coding_gene.bed',
  output:
    'species/elegans/gene_annotation/elegans.gene_annotation.coding_gene.fa',
    'kallisto_index/elegans'
  shell:
    '''
    fastaFromBed -fi {input[0]} -bed {input[1]} -nameOnly -split -s | awk 'BEGIN{{FS="(";}}{{print $1}}' > {output[0]}
    kallisto index -i {output[1]} {output[0]}
    '''


rule kallisto:
  input:
    'kallisto_index/elegans',
    'data/elegans/rnaseq/fastq/{sample}.r1.fq.gz',
    'species/elegans/gene_annotation/elegans.gene_annotation.gene_transcript_id.txt',
  output:
    'data/elegans/rnaseq/kallisto/{sample}/abundance.h5',
    'data/elegans/rnaseq/kallisto/{sample}/abundance.tsv',
    'data/elegans/rnaseq/kallisto/{sample}/run_info.json',
    'data/elegans/rnaseq/kallisto/{sample}/abundance.genes.txt',
  resources:
    cpus=10
  params:
    'data/elegans/rnaseq/kallisto/{sample}'
  shell:
    '''
    kallisto quant -i {input[0]} -b 100 --single -o {params} -l 200 -s 20 -t {resources.cpus} --rf-stranded {input[1]}
    python scripts/sum_TPM_per_gene.py {input[2]} {output[1]} {output[3]}
    sort -k 1,1 {output[3]} > {output[3]}.sorted
    mv {output[3]}.sorted {output[3]}
    '''


rule kallisto_final_table:
  input:
    kallisto_out=expand('data/elegans/rnaseq/kallisto/{sample}/abundance.genes.txt', sample=[item for item in config["rnaseq_samples"] if "adult" not in item]),
  output:
    'data/elegans/rnaseq/gene_expression.kallisto.txt',
  params:
    filenames=[item for item in config["rnaseq_samples"] if "adult" not in item],
  shell:
    '''
    echo {params.filenames} | tr " " "\t" > {output[0]}
    paste {input.kallisto_out} | awk '{{ for (i=3;i<=NF;i+=2) $i="" }} 1' | tr " " "\t" >> {output[0]}
    '''


rule deseq:
  input:
    expand('data/elegans/rnaseq/alignment/{sample}.Aligned.sortedByCoord.out.bam', sample=config["rnaseq_samples"]),
    'species/elegans/gene_annotation/elegans.annotations.coding_genes.gtf'
  output:
    'DE_analysis/genes_DESeq_him17_vs_N2.txt',
    'DE_analysis/genes_downregulated_him17_vs_N2.txt',
    'DE_analysis/genes_upregulated_him17_vs_N2.txt',
    'DE_analysis/normalized_counts_DESeq.txt',
    'plots/him17_vs_wt_volcano.pdf',
  shell:
    '''
    Rscript scripts/DESeq_mutants.R
    '''


rule chip_targets:
  input:
    'promoter_annotation/promoters_all.elegans.bed',
    'data/elegans/chipseq/yapc/{sample}.smooth_100_yapc_0.00001.bed',
  output:
    'DE_analysis/{sample}.bound_genes.txt'
  shell:
    '''
    intersectBed -a {input[0]} -b {input[1]} -u | cut -f 4 | tr ',' '\n' | sort | uniq | sed '1d' > {output[0]}
    '''


rule direct_targets:
  input:
    'DE_analysis/genes_downregulated_him17_vs_N2.txt',
    'DE_analysis/elegans_HIM17.bound_genes.txt',
  output:
    'DE_analysis/genes_downregulated_him17_vs_N2.direct.txt',
    'DE_analysis/genes_downregulated_him17_vs_N2.not_direct.txt',
  shell:
    '''
    join {input[0]} {input[1]} | awk '{{print $1}}' > {output[0]}
    join -v 1 {input[0]} {input[1]} | awk '{{print $1}}' > {output[1]}
    '''


rule him17_statistics:
  input:
    'data/elegans/chipseq/yapc/elegans_HIM17.smooth_100_yapc_0.00001.bed',
    'RE_annotation/reg_elements_all.elegans.gl_specific.bed',
    'RE_annotation/reg_elements_all.elegans.not_gl_specific.bed',
    'motif_enrichment/elegans/elegans.m1m2_clusters.bed',
  output:
    'HIM17_peaks/statistics.txt',
    'HIM17_peaks/m1m2_overlap.txt',
    'plots/HIM17_peaks_promoters_overlap.pdf',
    'plots/HIM17_peaks_m1m2_overlap.pdf',
  shell:
    '''
    echo -e "GL_promoters_m1m2\tGL_promoters_not_m1m2\tnotGL_promoters_m1m2\tnot_GL_promoters_not_m1m2" > {output[0]}
    GL_m1m2=$(grep coding_promoter {input[1]} | intersectBed -a stdin -b {input[3]} -u | wc -l)
    GL_no_m1m2=$(grep coding_promoter {input[1]} | intersectBed -a stdin -b {input[3]} -v | wc -l)
    no_GL_m1m2=$(grep coding_promoter {input[2]} | intersectBed -a stdin -b {input[3]} -u | wc -l)
    no_GL_no_m1m2=$(grep coding_promoter {input[2]} | intersectBed -a stdin -b {input[3]} -v | wc -l)
    echo -e "total\t"$GL_m1m2"\t"$GL_no_m1m2"\t"$no_GL_m1m2"\t"$no_GL_no_m1m2 >> {output[0]}
    GL_m1m2=$(grep coding_promoter {input[1]} | intersectBed -a stdin -b {input[0]} -u | intersectBed -a stdin -b {input[3]} -u | wc -l)
    GL_no_m1m2=$(grep coding_promoter {input[1]} | intersectBed -a stdin -b {input[0]} -u | intersectBed -a stdin -b {input[3]} -v | wc -l)
    no_GL_m1m2=$(grep coding_promoter {input[2]} | intersectBed -a stdin -b {input[0]} -u | intersectBed -a stdin -b {input[3]} -u | wc -l)
    no_GL_no_m1m2=$(grep coding_promoter {input[2]} | intersectBed -a stdin -b {input[0]} -u | intersectBed -a stdin -b {input[3]} -v | wc -l)
    echo -e "HIM17_peaks\t"$GL_m1m2"\t"$GL_no_m1m2"\t"$no_GL_m1m2"\t"$no_GL_no_m1m2 >> {output[0]}
    echo "loci" > {output[1]}
    m1m2_n=$(wc -l < {input[3]})
    him17_p_n=$(wc -l < {input[0]})
    overlap_m1m2=$(intersectBed -a {input[0]} -b {input[3]} -u | wc -l)
    echo -e "m1m2_pairs\t"$m1m2_n >> {output[1]}
    echo -e "him17_peaks\t"$him17_p_n >> {output[1]}
    echo -e "overlap\t"$overlap_m1m2 >> {output[1]}
    more {output[0]}
    more {output[1]}
    Rscript scripts/HIM17_peaks_stats.R
    '''
