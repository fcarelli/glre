SPECIES = ["elegans", "inopinata", "briggsae", "remanei", "nigoni", "contortus", "tipulae", "bacteriophora", "ceylanicum", "pacificus", "redivivus"]
ADDITIONAL_SPECIES = ["becei", "bovis", "japonica", "monodelphis", "plicata", "quiockensis", "uteleia", "exspectatus"]
ALL_SPECIES = ["elegans", "inopinata", "briggsae", "remanei", "nigoni", "contortus", "tipulae", "bacteriophora", "ceylanicum", "pacificus", "redivivus", "becei", "bovis", "japonica", "monodelphis", "plicata", "quiockensis", "uteleia", "exspectatus"]
REPEAT_SEARCH_SPECIES = ["elegans", "inopinata", "briggsae", "remanei", "nigoni", "bovis", "contortus", "tipulae", "pacificus"]
CAENORHABDITIS=["elegans", "inopinata", "briggsae", "remanei", "nigoni", "becei", "bovis", "japonica", "monodelphis", "plicata", "quiockensis", "uteleia"]
NOT_ELEGANS = ["inopinata", "briggsae", "remanei", "nigoni", "contortus", "tipulae", "bacteriophora", "ceylanicum", "pacificus", "redivivus", "becei", "bovis", "japonica", "monodelphis", "plicata", "quiockensis", "uteleia", "exspectatus"]
SPECIES_MOTIF_PLOTS = ["elegans", "contortus", "pacificus", "briggsae", "nigoni", "inopinata", "remanei", "ceylanicum"]
THAP_ANALYSIS_SPECIES = ["elegans", "inopinata", "briggsae", "remanei", "nigoni", "contortus", "tipulae", "bacteriophora", "ceylanicum", "pacificus", "redivivus", "becei", "bovis", "japonica", "monodelphis", "plicata", "quiockensis", "uteleia", "exspectatus", "placei", "caninum", "polygyrus", "brasiliensis", "viviparus"]

STAGES = ["embryo_series_1", "embryo_series_2", "embryo_series_3", "embryo_series_4", "postembryonic"]
EXTERNAL_ATAC = ["pgc", "adult_gl"]
REPS = ["rep1", "rep2"]
RNASEQ_SAMPLES=["wt", "him17", "xnd1"]
STRANDS=["1", "2"]

configfile: 'workflow/config.yaml'

include: 'workflow/species_annotation.snakefile'
include: 'workflow/RE_annotation.snakefile'
include: 'workflow/repeat_analysis.snakefile'
include: 'workflow/absense.snakefile'
include: 'workflow/him17_xnd1.snakefile'
include: 'workflow/briggsae_analysis.snakefile'
include: 'workflow/results_plots.snakefile'


rule all:
  input:
    expand('species/{sample}/genome/{sample}.chrom.sizes.txt', sample=ALL_SPECIES),
    'species/gene_annotation/elegans.operons.bed',
    expand('species/{sample}/protein_sequences/{sample}.protein.longest_isoform.fa', sample=SPECIES),
    expand('species/{sample}/gene_annotation/{sample}.gene_annotation.coding_gene_longest_isoform.bed', sample=SPECIES),
    expand('species/{sample}/protein_sequences/{sample}.protein.longest_isoform.fa', sample=ALL_SPECIES),
    #expand('species/{sample}/repeatmodeler/consensi.fa', sample=REPEAT_SEARCH_SPECIES),
    expand('motif_enrichment/{sample}/{sample}.m1m2_clusters.arrangements', sample=ALL_SPECIES),
    'RE_shuffling/coding_prom_GL.shuffling.results',
    expand('data/external_data/atac/peaks/elegans.{sample}.{rep}/elegans.{sample}.{rep}_treat_pileup.sorted.bw', sample=EXTERNAL_ATAC, rep=REPS),
    expand('data/external_data/atac/yapc/elegans.{sample}.smooth_100_yapc_0.001.bed', sample=EXTERNAL_ATAC),
    'plots/ATAC_coverage_GLRE_elegans.pdf',
    'RE_annotation/reg_elements_all.elegans.gl_specific.m1m2.bed',
    'RE_features/reg_elements_all.elegans.not_gl_specific.promoters.no_m1m2.cpg',
    'motif_enrichment/elegans/elegans.motif3_association.txt',
    expand('motif_enrichment/{sample}/{sample}.m3.bed', sample=ALL_SPECIES),
    expand('GL_genes_exp/promoters_not_gl_specific.elegans.no_m1m2.unique_promoter.{sample}.exp', sample=STAGES),
    'GL_genes_exp/promoters_gl_specific.elegans.m1m2.unique_promoter.embryo_germline.exp',
    expand('motif_enrichment/{sample}/{sample}.m1m2_clusters.bw', sample=ALL_SPECIES),
    'motif_enrichment/elegans/repeat_m1m2_overlap.summary',
    'motif_enrichment/elegans/repeat_GL_promoter_overlap.summary',
    'gene_annotation/elegans.GL_specific_genes.m1m2_unique.genes',
    'gene_annotation/briggsae.genes_DESeq_glp1_vs_wt.txt',
    'data/external_data/boeck_unified_exp_data.txt',
    'meme/elegans.GL_specific_m1m2_vs_nonGL_specific_m1m2',
    'plots/current_vs_serizay_RE_annotation.heatmap.pdf',
    'data/external_data/essential_genes.txt',
    'de_novo_mites/de_novo_mites.all.aln',
    expand('TF_evolution/orthologs/HIM17.{sample}.aln', sample = ALL_SPECIES),
    expand('TF_evolution/orthologs/XND1.{sample}.aln', sample = ALL_SPECIES),
    expand('TF_evolution/THAP_conservation/{sample}.THAP.parsed', sample = THAP_ANALYSIS_SPECIES),
    expand('TF_evolution/DBD_conservation/XND1.{sample}.DBD.parsed', sample = ALL_SPECIES),
    expand('TF_evolution/DBD_conservation/HIM17.{sample}.DBD.parsed', sample = ALL_SPECIES),
    'TF_evolution/orthologs/HIM17.caenorhabditis.pacificus.hits',
    'data/external_data/elegans.tipulae.1to1.txt',
    'absense_run/concatenated_orthologs.dist',
    'absense_run/HIM17.absense.pdf',
    'absense_run/XND1.absense.pdf',
    'plots/HIM17_identity.heat.pdf',
    'plots/XND1_identity.heat.pdf',
    'plots/HIM17_XND1_orthologs_identity.pdf',
    #expand('species/{sample}/mitehunter/{sample}_mitehunter.mite.all.fa', sample=REPEAT_SEARCH_SPECIES),
    'modern_analysis/annotatedPeak.peaks_GLRE.m1m2.enrichment',
    expand('data/elegans/chipseq/yapc/{sample}.smooth_100_yapc_coverage.bw', sample=config["chip_sample_names"]),
    'chip_stats/chip_stats.txt',
    'data/elegans/chipseq/correlation/corr_HIM17_XND1.spearman.scatter.pdf',
    expand('data/elegans/rnaseq/fastq/{sample}.r1.fq.gz', sample=config["rnaseq_samples"]),
    expand('chip_stats/{sample}.motif_stats.txt', sample=config["chip_sample_names"]),
    expand('meme/{sample}.nomotif.enrichment', sample=config["chip_sample_names"]),
    expand('data/elegans/rnaseq/alignment/{sample}.Signal.UniqueMultiple.str1.out.bw', sample=config["rnaseq_samples"]),
    expand('data/elegans/rnaseq/tracks/elegans_{sample}_merged.Signal.UniqueMultiple.str{strand}.out.bw', sample=RNASEQ_SAMPLES, strand=STRANDS),
    expand('data/elegans/rnaseq/kallisto/{sample}/abundance.genes.txt', sample=config["rnaseq_samples"]),
    expand('data/elegans/rnaseq/kallisto/{sample}/abundance.h5', sample=config["rnaseq_samples"]),
    'data/elegans/rnaseq/gene_expression.kallisto.txt',
    'DE_analysis/normalized_counts_DESeq.txt',
    'DE_analysis/genes_downregulated_xnd1_vs_N2.not_direct.txt',
    'RE_features/briggsae.motif_GL_overlap.summary',
    'promoter_annotation/promoters_all.briggsae.bed',
    'meme/briggsae.GL_specific_vs_nonGL_specific',
    'meme/briggsae.GL_specific_m1m2_vs_nonGL_specific_m1m2',
    'RE_features/reg_elements_all.briggsae.gl_specific.promoters.m1m2.nuc',
    'species/elegans/orthologs/one2one.briggsae.orthologs.txt',
    'RE_conservation/conservation_statistics.txt',
    'RE_conservation_all/conservation_statistics.txt',
    'species/briggsae/repeatmasker/rm_repeatmasker/briggsae.fa.out',
    'species/briggsae/repeatmasker/mh_repeatmasker/briggsae.fa.out',
    'RE_conservation/promoters_GL_m1m2.briggsae.by_type.summary',
    'RE_conservation/promoters_GL_m1m2.elegans.by_type.summary',
    expand('data/external_data/rnaseq/elegans/kallisto/{sample}/abundance.genes.txt', sample=config['elegans_external_rnaseq_samples']),
    expand('data/external_data/rnaseq/briggsae/kallisto/{sample}/abundance.genes.txt', sample=config['briggsae_external_rnaseq_samples']),
    'data/external_data/rnaseq/briggsae/gene_expression.stages.txt',
    'data/external_data/rnaseq/elegans/gene_expression.stages.txt',
    'data/external_data/rnaseq/elegans/gene_expression.pgc.txt',
    'data/external_data/rnaseq/contortus/gene_expression.stages.txt',
    'plots/HIM17_XND1_CELE2_CERP2_m1m2_RE.heatmap.pdf',
    'heatmaps_int_files/HIM17_XND1_CELE2_CERP2_m1m2_RE.sorted.clusters',
    'motif_conservation/elegans.tandem_m2m1.GLRE.heat.pdf',
    'plots/RE_ATAC.elegans.heatmap.pdf',
    'plots/RE_m1m2_ATAC.heatmap.pdf',
    'test_him17_cov/promoters_elegans_HIM17.heatmap.pdf',
    'plots/HIM17_XND1.all_peaks.heatmap.pdf',
    'plots/RE_m1m2_external_ATAC.heatmap.pdf',
    expand('plots/{sample}.m1m2.all_genes.profile.pdf', sample=SPECIES_MOTIF_PLOTS),
    'plots/GL_promoters_m1m2_external_ATAC.heatmap.pdf',
    'plots/HIM17_XND1.all_promoters.heatmap.pdf',
    expand('plots/{sample}.m3.all_genes.profile.pdf', sample=SPECIES_MOTIF_PLOTS),
    'plots/qPCR_HIM17_XND1_pct_input.pdf',
    'plots/mut_vs_wt_expression.m1m2_promoters.pdf',
    'intact_repeats/elegans/elegans_CERP2.inactive.clw.hmm',
    'intact_repeats/elegans/elegans_CELE2.intact.clw.hmm',
    'intact_repeats/inactive_CERP2_elegans.hmm',
    'intact_repeats/elegans/elegans_CERP2.no_HIM17.clw.hmm',
    'plots/TT_periodicity_CERP2_elements.pdf',
    'plots/TT_periodicity_CELE2_elements.pdf',
    'plots/XND1_binding_other_repeats.pdf',
    'plots/m1m2_permutation_test.pdf',
    #expand('de_novo_mites/{sample}_mitehunter.mite.m1m2_pairs.canonical.CERP2.aln', sample=REPEAT_SEARCH_SPECIES),
    'plots/HIM17_XND1_TE_m1m2_RE.heatmap.repeat_sorted.pdf',
    'shared_peaks/statistics.txt',
    'RE_conservation_all/conserved_vs_species_specific_promoters.m1m2_statistics.txt',
    expand('species/{sample}/repeats_dfam/{sample}_dfam.repeats.id.bed', sample=REPEAT_SEARCH_SPECIES),
    


#    expand('intact_repeats/{sample}/{sample}_CERP2_all.clw.hmm', sample=ALL_SPECIES),
#    'intact_repeats/inactive_CERP2_elegans.hmm',
#    'intact_repeats/elegans/elegans_CERP2.inactive.clw.hmm',
    

    






#'gene_annotation/elegans.germline_expressed.serizay.bed'
