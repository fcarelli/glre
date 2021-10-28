rule misc_R_plots:
  input:
    expand('GL_genes_exp/promoters_not_gl_specific.elegans.no_m1m2.unique_promoter.{sample}.exp', sample=STAGES),
    'data/external_data/rnaseq/elegans/gene_expression.pgc.txt',
    'data/external_data/rnaseq/elegans/gene_expression.stages.txt',
    'promoter_annotation/promoters_gl_specific.elegans.m1m2.unique_promoter.genes',
    'motif_enrichment/elegans/elegans.motif3_association.txt',
    'RE_features/reg_elements_all.elegans.gl_specific.promoters.m1m2.nuc',
    'RE_features/reg_elements_all.briggsae.gl_specific.promoters.m1m2.nuc',
    'RE_features/elegans.motif_GL_overlap.summary',
    'RE_features/briggsae.motif_GL_overlap.summary',
    expand('motif_enrichment/{sample}/{sample}.m1m2_clusters.arrangements_all.summary', sample=ALL_SPECIES),
    'motif_enrichment/elegans/elegans.m1m2_clusters.arrangements_by_distance.summary',
    'motif_enrichment/elegans/repeat_m1m2_overlap.summary',
    'motif_enrichment/elegans/repeat_GL_RE_overlap.summary',
    'HIM17_qPCR/HIM17_pct_input.txt',
    'DE_analysis/genes_DESeq_him17_vs_N2.txt',
    'gene_annotation/elegans.GL_specific_genes.m1m2_unique.genes',
    'DE_analysis/elegans_HIM17.bound_genes.txt',
    'species/elegans/repeats_dfam/elegans_dfam.repeats.id.bed',
    'DE_analysis/genes_downregulated_him17_vs_N2.direct.txt',
    expand('CELE2_CERP2_other_species/{sample}.summary', sample=REPEAT_SEARCH_SPECIES),
  output:
    'plots/expression_unique_GL_promoter_genes.pdf',
    'plots/GC_content.pdf',
    'plots/RE_m1m2_overlap.m1m2_arrangement.pdf',
    'plots/m1m2_motif_distance_divergent.other_species.pdf',
    'plots/m1m2_motif_arrangements.all_species.pdf',
    'plots/m1m2_motif_arrangements.elegans.pdf',
    'plots/qPCR_T05F1.2_fold_change.pdf',
    'plots/qPCR_HIM17_pct_input.pdf',
    'plots/mut_vs_wt_expression.m1m2_promoters.pdf',
    'plots/CERP2_CELE2_length.pdf',
    'plots/CERP2_CELE2_n_other_species.pdf',
    'plots/GO_enrichment_direct_targets.pdf',
  shell:
    '''
    Rscript scripts/exp_analysis_m1m2_genes.R
    Rscript scripts/GC_analysis.R
    Rscript scripts/m1m2_RE_overlap.R
    Rscript scripts/motif_by_species_statistics.R
    Rscript scripts/motif_RE_statistics.R
    Rscript scripts/qpcr_enrichment.R
    Rscript scripts/mutant_wt_expression.R
    Rscript scripts/CERP2_CELE2_len.R
    Rscript scripts/CERP2_CELE2_other_species.R
    Rscript scripts/GO_cluster_profiler.R
    '''

rule him17_peaks_heatmap:
  input:
    "data/elegans/chipseq/yapc/elegans_HIM17.smooth_100_yapc_0.00001.bed",
    "species/elegans/repeats_dfam/elegans_dfam.repeats.id.bed",
    "motif_enrichment/elegans/elegans.m1m2_clusters.bed",
    "RE_annotation/reg_elements_all.elegans.gl_specific.bed",
    "RE_annotation/reg_elements_all.elegans.not_gl_specific.bed",
    "species/elegans/genome/elegans.chrom.sizes.txt",
    "data/elegans/chipseq/yapc/elegans_HIM17.smooth_100_yapc_coverage.bw",
  output:
    "heatmaps_int_files/HIM17_peaks.all.bed",
    temp("heatmaps_int_files/CELE2_CERP2.bg"),
    temp("heatmaps_int_files/m1m2.bg"),
    temp("heatmaps_int_files/RE_GL_specific.bg"),
    temp("heatmaps_int_files/RE_nonGL_specific.bg"),
    "heatmaps_int_files/CELE2_CERP2.bw",
    "heatmaps_int_files/m1m2.bw",
    "heatmaps_int_files/RE_GL_specific.bw",
    "heatmaps_int_files/RE_nonGL_specific.bw",
    "heatmaps_int_files/HIM17_CELE2_CERP2_m1m2_RE.heatmap.bed",
    "heatmaps_int_files/HIM17_CELE2_CERP2_m1m2_RE.sorted.mat.gz",
    "plots/HIM17_CELE2_CERP2_m1m2_RE.heatmap.repeat_sorted.pdf",
    "heatmaps_int_files/HIM17_CELE2_CERP2_m1m2_RE.sorted.clusters",
  resources:
    cpus=10
  shell:
    '''
    cut -f 1,2,3 {input[0]} | sort -k 1,1 -k2,2n > {output[0]}
    grep 'CERP2\|CELE2' {input[1]} | genomeCoverageBed -i stdin -g {input[5]} -bga | LC_COLLATE=C sort -k 1,1 -k2,2n > {output[1]}
    genomeCoverageBed -i {input[2]} -g {input[5]} -bga | LC_COLLATE=C sort -k 1,1 -k2,2n > {output[2]}
    genomeCoverageBed -i {input[3]} -g {input[5]} -bga | LC_COLLATE=C sort -k 1,1 -k2,2n > {output[3]}
    genomeCoverageBed -i {input[4]} -g {input[5]} -bga | LC_COLLATE=C sort -k 1,1 -k2,2n > {output[4]}
    bedGraphToBigWig {output[1]} {input[5]} {output[5]}
    bedGraphToBigWig {output[2]} {input[5]} {output[6]}
    bedGraphToBigWig {output[3]} {input[5]} {output[7]}
    bedGraphToBigWig {output[4]} {input[5]} {output[8]}
    intersectBed -a {output[0]} -b {input[2]} -u > {output[0]}.m1m2
    intersectBed -a {output[0]} -b {input[2]} -v > {output[0]}.no_m1m2
    echo -e "#chrom\tstart\tend" > {output[9]}
    for cluster in m1m2 no_m1m2
    do
    intersectBed -a {output[0]}.$cluster -b RE_annotation/reg_elements_all.elegans.gl_specific.bed -u | awk '{{print $0 "\t1" }}' > {output[0]}.$cluster.glre
    intersectBed -a {output[0]}.$cluster -b RE_annotation/reg_elements_all.elegans.gl_specific.bed -v | awk '{{print $0 "\t0" }}' >> {output[0]}.$cluster.glre
    intersectBed -a {output[0]}.$cluster.glre -b RE_annotation/reg_elements_all.elegans.not_gl_specific.bed -u | awk '{{print $0 "\t1" }}' > {output[0]}.$cluster.glre.re
    intersectBed -a {output[0]}.$cluster.glre -b RE_annotation/reg_elements_all.elegans.not_gl_specific.bed -v | awk '{{print $0 "\t0" }}' >> {output[0]}.$cluster.glre.re
    grep 'CERP2\|CELE2' {input[1]} | intersectBed -a {output[0]}.$cluster.glre.re -b stdin -u | awk '{{print $0 "\t1" }}' > {output[0]}.$cluster.bed.glre.re.TE
    grep 'CERP2\|CELE2' {input[1]} | intersectBed -a {output[0]}.$cluster.glre.re -b stdin -v | awk '{{print $0 "\t0" }}' >> {output[0]}.$cluster.bed.glre.re.TE
    sort -k4,4rn -k5,5rn -k6,6rn {output[0]}.$cluster.bed.glre.re.TE | cut -f 1,2,3 >> {output[9]}
    echo -e "#chrom\tstart\tend" >> {output[9]}
    done
    computeMatrix scale-regions -R {output[9]} -S {input[6]} {output[5]} {output[6]} {output[7]} {output[8]} --beforeRegionStartLength 500 --afterRegionStartLength 500 --regionBodyLength 200 --startLabel peak_start --endLabel peak_end --binSize 10 --samplesLabel HIM17 CELE2_CERP2 m1m2 GLRE nonGLRE -p {resources.cpus} -o {output[10]} 
    plotHeatmap -m {output[10]} -out {output[11]} --sortRegions keep --averageTypeSummaryPlot mean --missingDataColor 0.5 --colorList 'white,forestgreen' 'white,black' 'white,darkorange' 'white,purple' 'white,grey' --plotTitle "HIM-17 peaks" -min 0 -max 60 1 1 1 1 --whatToShow 'heatmap and colorbar' --interpolationMethod gaussian --outFileSortedRegions {output[12]}
    '''


rule phyloP_download:
  output:
    'data/external_data/phyloP_26way.bw'
  shell:
    '''
    wget ftp://hgdownload.soe.ucsc.edu/goldenPath/ce11/phyloP26way/ce11.phyloP26way.bw -O {output[0]}
    '''


rule m1m2_conservation_phylop:
  input:
    'motif_enrichment/elegans/elegans.m1m2_clusters.bed',
    'RE_annotation/reg_elements_all.elegans.bed',
    'RE_annotation/reg_elements_all.elegans.gl_specific.bed',
    'data/external_data/phyloP_26way.bw',
  output:
    'motif_conservation/elegans.canonical_divergent_m1m2.bed',
    'motif_conservation/elegans.canonical_tandem_m2m1.bed',
    'motif_conservation/elegans.divergent.all_RE.bed',
    'motif_conservation/elegans.divergent.GLRE.bed',
    'motif_conservation/elegans.tandem_m2m1.all_RE.bed',
    'motif_conservation/elegans.tandem_m2m1.GLRE.bed',
    'motif_conservation/elegans.divergent.all_RE.heat.pdf',
    'motif_conservation/elegans.divergent.GLRE.heat.pdf',
    'motif_conservation/elegans.tandem_m2m1.all_RE.heat.pdf',
    'motif_conservation/elegans.tandem_m2m1.GLRE.heat.pdf',
  resources:
    cpus=6
  shell:
    '''
    grep divergent {input[0]} | awk 'BEGIN{{FS="_|bp";}}{{if ($2 > 10 && $2 < 20) print}}' > {output[0]}
    grep tandem_m2m1 {input[0]} | awk 'BEGIN{{FS="_|bp";}}{{if ($2 > 20) print}}' > {output[1]}
    intersectBed -a {input[1]} -b {output[1]} -v | intersectBed -a {output[0]} -b stdin -u > {output[2]}
    intersectBed -a {input[2]} -b {output[1]} -v | intersectBed -a {output[0]} -b stdin -u > {output[3]}
    intersectBed -a {input[1]} -b {output[0]} -v | intersectBed -a {output[1]} -b stdin -u > {output[4]}
    intersectBed -a {input[2]} -b {output[0]} -v | intersectBed -a {output[1]} -b stdin -u > {output[5]}
    computeMatrix scale-regions -R {output[2]} -S {input[3]} --beforeRegionStartLength 50 --afterRegionStartLength 50 --regionBodyLength 50 --startLabel motif_start --endLabel motif_end --binSize 1 --samplesLabel divergent.all_RE -o {output[6]}.mat.gz -p {resources.cpus}
    computeMatrixOperations filterStrand -m {output[6]}.mat.gz -o {output[6]}.for.mat.gz --strand +
    computeMatrixOperations filterStrand -m {output[6]}.mat.gz -o {output[6]}.rev.mat.gz --strand -
    computeMatrixOperations rbind -m {output[6]}.for.mat.gz {output[6]}.rev.mat.gz -o {output[6]}.stranded.mat.gz
    plotHeatmap -m {output[6]}.stranded.mat.gz -out {output[6]} --averageTypeSummaryPlot mean --missingDataColor 0.5 --plotTitle "phyloP divergent_m1m2 RE" --colorMap Blues --heatmapHeight 12 --heatmapWidth 6 --yAxisLabel "divergent_m2m1" --xAxisLabel "" --startLabel "m1" --endLabel "m2"
    computeMatrix scale-regions -R {output[3]} -S {input[3]} --beforeRegionStartLength 50 --afterRegionStartLength 50 --regionBodyLength 50 --startLabel motif_start --endLabel motif_end --binSize 1 --samplesLabel divergent.GLRE -o {output[7]}.mat.gz -p {resources.cpus}
    computeMatrixOperations filterStrand -m {output[7]}.mat.gz -o {output[7]}.for.mat.gz --strand +
    computeMatrixOperations filterStrand -m {output[7]}.mat.gz -o {output[7]}.rev.mat.gz --strand -
    computeMatrixOperations rbind -m {output[7]}.for.mat.gz {output[7]}.rev.mat.gz -o {output[7]}.stranded.mat.gz
    plotHeatmap -m {output[7]}.stranded.mat.gz -out {output[7]} --averageTypeSummaryPlot mean --missingDataColor 0.5 --plotTitle "phyloP divergent_m1m2 GL RE" --colorMap Purples --heatmapHeight 12 --heatmapWidth 6 --yAxisLabel "divergent_m2m1" --xAxisLabel "" --startLabel "m1" --endLabel "m2"
    computeMatrix scale-regions -R {output[4]} -S {input[3]} --beforeRegionStartLength 50 --afterRegionStartLength 50 --regionBodyLength 50 --startLabel motif_start --endLabel motif_end --binSize 1 --samplesLabel tandem_m1m2.all_RE -o {output[8]}.mat.gz -p {resources.cpus}
    computeMatrixOperations filterStrand -m {output[8]}.mat.gz -o {output[8]}.for.mat.gz --strand +
    computeMatrixOperations filterStrand -m {output[8]}.mat.gz -o {output[8]}.rev.mat.gz --strand -
    computeMatrixOperations rbind -m {output[8]}.for.mat.gz {output[8]}.rev.mat.gz -o {output[8]}.stranded.mat.gz
    plotHeatmap -m {output[8]}.stranded.mat.gz -out {output[8]} --averageTypeSummaryPlot mean --missingDataColor 0.5 --plotTitle "phyloP tandem_m1m2 all RE" --colorMap Blues --heatmapHeight 12 --heatmapWidth 6 --yAxisLabel "tandem_m2m1" --xAxisLabel "" --startLabel "m1" --endLabel "m2"
    computeMatrix scale-regions -R {output[5]} -S {input[3]} --beforeRegionStartLength 50 --afterRegionStartLength 50 --regionBodyLength 50 --startLabel motif_start --endLabel motif_end --binSize 1 --samplesLabel tandem_m1m2.GLRE -o {output[9]}.mat.gz -p {resources.cpus}
    computeMatrixOperations filterStrand -m {output[9]}.mat.gz -o {output[9]}.for.mat.gz --strand +
    computeMatrixOperations filterStrand -m {output[9]}.mat.gz -o {output[9]}.rev.mat.gz --strand -
    computeMatrixOperations rbind -m {output[9]}.for.mat.gz {output[9]}.rev.mat.gz -o {output[9]}.stranded.mat.gz
    plotHeatmap -m {output[9]}.stranded.mat.gz -out {output[9]} --averageTypeSummaryPlot mean --missingDataColor 0.5 --plotTitle "phyloP tandem_m1m2 GL RE" --colorMap Purples --heatmapHeight 12 --heatmapWidth 6 --yAxisLabel "tandem_m2m1" --xAxisLabel "" --startLabel "m1" --endLabel "m2"
    '''


rule RE_ATAC_heatmap:
  input:
    'RE_annotation/reg_elements_all.elegans.gl_specific.bed',
    'RE_annotation/reg_elements_all.elegans.not_gl_specific.bed',
    'relmapping/processed_tracks/elegans_atac_ya_wt.bw',
    'relmapping/processed_tracks/elegans_atac_ya_glp1.bw',
    'RE_annotation/reg_elements_all.briggsae.gl_specific.bed',
    'RE_annotation/reg_elements_all.briggsae.not_gl_specific.bed',
    'relmapping/processed_tracks/briggsae_atac_ya_wt.bw',
    'relmapping/processed_tracks/briggsae_atac_ya_glp1.bw',
  output:
    'heatmaps_int_files/RE_ATAC.elegans.mat.gz',
    'plots/RE_ATAC.elegans.heatmap.pdf',
    'heatmaps_int_files/RE_ATAC.briggsae.mat.gz',
    'plots/RE_ATAC.briggsae.heatmap.pdf',
  resources:
    cpus=10
  shell:
    '''
    computeMatrix scale-regions -R {input[0]} {input[1]} -S {input[2]} {input[3]} --beforeRegionStartLength 500 --afterRegionStartLength 500 --regionBodyLength 200 --startLabel peak_start --endLabel peak_end --binSize 10 --samplesLabel wt glp-1 -o {output[0]} -p {resources.cpus}
    plotHeatmap -m {output[0]} -out {output[1]} --averageTypeSummaryPlot mean --missingDataColor 0.5 --colorList 'white,purple' 'white,black' --plotTitle "elegans RE" -min 0 -max 10 10 --whatToShow 'heatmap and colorbar' --heatmapHeight 12 --heatmapWidth 6 --xAxisLabel "" --startLabel peak_start --endLabel peak_end
    computeMatrix scale-regions -R {input[4]} {input[5]} -S {input[6]} {input[7]} --beforeRegionStartLength 500 --afterRegionStartLength 500 --regionBodyLength 200 --startLabel peak_start --endLabel peak_end --binSize 10 --samplesLabel wt glp-1 -o {output[2]} -p {resources.cpus}
    plotHeatmap -m {output[2]} -out {output[3]} --averageTypeSummaryPlot mean --missingDataColor 0.5 --colorList 'white,purple' 'white,black' --plotTitle "briggsae RE" -min 0 -max 10 10 --whatToShow 'heatmap and colorbar' --heatmapHeight 12 --heatmapWidth 6 --xAxisLabel "" --startLabel peak_start --endLabel peak_end
    '''


rule RE_m1m2_ATAC_heatmap:
  input:
    'RE_annotation/reg_elements_all.elegans.gl_specific.m1m2.bed',
    'RE_annotation/reg_elements_all.elegans.not_gl_specific.m1m2.bed',
    'relmapping/processed_tracks/elegans_atac_ya_wt.bw',
    'relmapping/processed_tracks/elegans_atac_ya_glp1.bw',
  output:
    'heatmaps_int_files/RE_m1m2_ATAC.mat.gz',
    'plots/RE_m1m2_ATAC.heatmap.pdf',
  resources:
    cpus=10
  shell:
    '''
    computeMatrix scale-regions -R {input[0]} {input[1]} -S {input[2]} {input[3]} --beforeRegionStartLength 500 --afterRegionStartLength 500 --regionBodyLength 200 --startLabel peak_start --endLabel peak_end --binSize 10 --samplesLabel wt glp-1 -o {output[0]} -p {resources.cpus}
    plotHeatmap -m {output[0]} -out {output[1]} --averageTypeSummaryPlot mean --missingDataColor 0.5 --colorList 'white,purple' 'white,black' --plotTitle "elegans RE" -min 0 -max 20 20 --whatToShow 'heatmap and colorbar' --heatmapHeight 12 --heatmapWidth 6 --xAxisLabel "" --startLabel peak_start --endLabel peak_end
    '''


rule m1m2_coding_gene_profiles:
  input:
    'motif_enrichment/{sample}/{sample}.m1m2_clusters.bed',
    'species/{sample}/genome/{sample}.chrom.sizes.txt',
    'species/{sample}/gene_annotation/{sample}.gene_annotation.coding_gene_longest_isoform.bed',
  output:
    temp('heatmaps_int_files/{sample}.m1m2.bg'),
    'heatmaps_int_files/{sample}.m1m2.bw',
    'heatmaps_int_files/{sample}.m1m2.mat.gz',
    'plots/{sample}.m1m2.all_genes.profile.pdf',
    temp('heatmaps_int_files/{sample}.divergent_m1m2.bg'),
    'heatmaps_int_files/{sample}.divergent_m1m2.bw',
    'heatmaps_int_files/{sample}.divergent.m1m2.mat.gz',
    'plots/{sample}.divergent_m1m2.all_genes.profile.pdf',
  resources:
    cpus=10
  shell:
    '''
    genomeCoverageBed -i {input[0]} -g {input[1]} -bga | LC_COLLATE=C sort -k 1,1 -k2,2n > {output[0]}
    bedGraphToBigWig {output[0]} {input[1]} {output[1]}
    computeMatrix scale-regions -R {input[2]} -S {output[1]} --beforeRegionStartLength 2000 --afterRegionStartLength 1000 --regionBodyLength 2000 --binSize 50 --samplesLabel m1m2  -p {resources.cpus} -o {output[2]}
    plotProfile -m {output[2]} -out {output[3]} --plotType fill --colors dodgerblue --plotTitle "{wildcards.sample} m1m2 coding genes" --plotHeight 5 --plotWidth 8 --yMin 0 --yMax 0.007
    grep divergent {input[0]} | genomeCoverageBed -i stdin -g {input[1]} -bga | LC_COLLATE=C sort -k 1,1 -k2,2n > {output[4]}
    bedGraphToBigWig {output[4]} {input[1]} {output[5]}
    computeMatrix scale-regions -R {input[2]} -S {output[5]} --beforeRegionStartLength 2000 --afterRegionStartLength 1000 --regionBodyLength 2000 --binSize 10 --samplesLabel m1m2  -p {resources.cpus} -o {output[6]}
    plotProfile -m {output[6]} -out {output[7]} --plotType fill --colors dodgerblue --plotTitle "{wildcards.sample} divergent m1m2 coding genes" --plotHeight 5 --plotWidth 8 --yMin 0 --yMax 0.007
    '''


rule m1_m2_coding_gene_profiles:
  input:
    'motif_enrichment/{sample}/{sample}.m3.bed',
    'species/{sample}/genome/{sample}.chrom.sizes.txt',
    'species/{sample}/gene_annotation/{sample}.gene_annotation.coding_gene_longest_isoform.bed',
  output:
    temp('heatmaps_int_files/{sample}.m3.bg'),
    'heatmaps_int_files/{sample}.m3.bw',
    'heatmaps_int_files/{sample}.m3.mat.gz',
    'plots/{sample}.m3.all_genes.profile.pdf',
  resources:
    cpus=10
  shell:
    '''
    genomeCoverageBed -i {input[0]} -g {input[1]} -bga | LC_COLLATE=C sort -k 1,1 -k2,2n > {output[0]}
    bedGraphToBigWig {output[0]} {input[1]} {output[1]}
    computeMatrix scale-regions -R {input[2]} -S {output[1]} --beforeRegionStartLength 2000 --afterRegionStartLength 1000 --regionBodyLength 2000 --binSize 50 --samplesLabel m1m2  -p {resources.cpus} -o {output[2]}
    plotProfile -m {output[2]} -out {output[3]} --plotType fill --colors dodgerblue --plotTitle "{wildcards.sample} m3 coding genes" --plotHeight 5 --plotWidth 8
    '''


rule RE_m1m2_external_ATAC_heatmap:
  input:
    'RE_annotation/reg_elements_all.elegans.gl_specific.m1m2.bed',
    'RE_annotation/reg_elements_all.elegans.gl_specific.bed',
    'relmapping/processed_tracks/elegans_atac_ya_wt.bw',
    'relmapping/processed_tracks/elegans_atac_ya_glp1.bw',
    'data/external_data/atac/yapc/elegans.pgc.smooth_100_yapc_coverage.bw',
    'RE_annotation/reg_elements_all.elegans.not_gl_specific.m1m2.bed',
    'data/external_data/atac/yapc/elegans.adult_gl.smooth_100_yapc_coverage.bw',
  output:
    temp('RE_annotation/reg_elements_all.elegans.gl_specific.no_m1m2.bed'),
    'heatmaps_int_files/RE_m1m2_external_ATAC.mat.gz',
    'plots/RE_m1m2_external_ATAC.heatmap.pdf',
  resources:
    cpus=10
  shell:
    '''
    intersectBed -a {input[1]} -b {input[0]} -v > {output[0]}
    computeMatrix scale-regions -R {input[0]} {output[0]} {input[5]} -S {input[2]} {input[3]} {input[4]} {input[6]} --beforeRegionStartLength 500 --afterRegionStartLength 500 --regionBodyLength 200 --startLabel peak_start --endLabel peak_end --binSize 10 --samplesLabel wt glp-1 pgc adult_gl -o {output[1]} -p {resources.cpus}
    plotHeatmap -m {output[1]} -out {output[2]} --averageTypeSummaryPlot median --missingDataColor 0.5 --colorList 'white,purple' 'white,black' 'white,darkblue' 'white,royalblue' --plotTitle "elegans RE" -min 0 -max 20 20 10 10 --heatmapHeight 15 --heatmapWidth 8 --xAxisLabel "" --startLabel peak_start --endLabel peak_end
    '''


rule GLRE_m1m2_external_ATAC_heatmap:
  input:
    'RE_annotation/reg_elements_all.elegans.gl_specific.m1m2.bed',
    'data/external_data/atac/yapc/elegans.pgc.smooth_100_yapc_coverage.bw',
    'data/external_data/atac/yapc/elegans.adult_gl.smooth_100_yapc_coverage.bw',
  output:
    temp('RE_annotation/reg_elements_all.elegans.gl_specific.m1m2.promoters.bed'),
    'heatmaps_int_files/GLRE_m1m2_external_ATAC.mat.gz',
    'plots/GL_promoters_m1m2_external_ATAC.heatmap.pdf',
  resources:
    cpus=10
  shell:
    '''
    grep coding_promoter {input[0]} > {output[0]}
    computeMatrix scale-regions -R {output[0]} -S {input[1]} {input[2]} --beforeRegionStartLength 500 --afterRegionStartLength 500 --regionBodyLength 200 --startLabel peak_start --endLabel peak_end --binSize 10 --samplesLabel pgc adult_gl -o {output[1]} -p {resources.cpus}
    plotHeatmap -m {output[1]} -out {output[2]} --kmeans 3 --averageTypeSummaryPlot median --missingDataColor 0.5 --colorList 'white,darkblue' 'white,royalblue' --plotTitle "elegans m1m2 GL promoters" -min 0 -max 10 10 --heatmapHeight 10 --heatmapWidth 6 --xAxisLabel "" --startLabel peak_start --endLabel peak_end
    '''
