SPECIES = ["elegans", "inopinata", "briggsae", "remanei", "nigoni", "contortus", "tipulae", "bacteriophora", "ceylanicum", "pacificus", "redivivus"]
ADDITIONAL_SPECIES = ["becei", "bovis", "japonica", "monodelphis", "plicata", "quiockensis", "uteleia", "exspectatus"]
ALL_SPECIES = ["elegans", "inopinata", "briggsae", "remanei", "nigoni", "contortus", "tipulae", "bacteriophora", "ceylanicum", "pacificus", "redivivus", "becei", "bovis", "japonica", "monodelphis", "plicata", "quiockensis", "uteleia", "exspectatus"]
REPEAT_SEARCH_SPECIES = ["elegans", "inopinata", "briggsae", "remanei", "nigoni", "bovis", "contortus", "pacificus"]

STAGES = ["embryo_series_1", "embryo_series_2", "embryo_series_3", "embryo_series_4", "postembryonic"]
EXTERNAL_ATAC = ["pgc", "adult_gl"]
REPS = ["rep1", "rep2"]


configfile: 'workflow/config.yaml'

include: 'workflow/species_annotation.snakefile'
include: 'workflow/RE_annotation.snakefile'


rule all:
  input:
    expand('species/{sample}/genome/{sample}.chrom.sizes.txt', sample=ALL_SPECIES),
    'species/gene_annotation/elegans.operons.bed',
    expand('species/{sample}/protein_sequences/{sample}.protein.longest_isoform.fa', sample=SPECIES),
    expand('species/{sample}/gene_annotation/{sample}.gene_annotation.coding_gene_longest_isoform.bed', sample=SPECIES),
    expand('species/{sample}/protein_sequences/{sample}.protein.longest_isoform.fa', sample=ALL_SPECIES),
    expand('species/{sample}/repeatmodeler/consensi.fa', sample=REPEAT_SEARCH_SPECIES),
    expand('motif_enrichment/{sample}/{sample}.m1m2_clusters.arrangements', sample=ALL_SPECIES),
    expand('data/external_data/atac/peaks/elegans.{sample}.{rep}/elegans.{sample}.{rep}_treat_pileup.sorted.bw', sample=EXTERNAL_ATAC, rep=REPS),
    expand('data/external_data/atac/yapc/elegans.{sample}.smooth_100_yapc_0.001.bed', sample=EXTERNAL_ATAC),
    'plots/ATAC_coverage_GLRE_elegans.pdf',
    'RE_annotation/reg_elements_all.elegans.gl_specific.m1m2.bed',
    'meme/elegans.GL_specific_m1m2_vs_nonGL_specific_m1m2',
    'RE_features/reg_elements_all.elegans.not_gl_specific.promoters.no_m1m2.cpg',
    'motif_enrichment/elegans/elegans.motif3_association.txt',
    expand('GL_genes_exp/promoters_not_gl_specific.elegans.no_m1m2.unique_promoter.{sample}.exp', sample=STAGES),
    'GL_genes_exp/promoters_gl_specific.elegans.m1m2.unique_promoter.embryo_germline.exp',
    expand('motif_enrichment/{sample}/{sample}.m1m2_clusters.bw', sample=ALL_SPECIES),
    'motif_enrichment/elegans/repeat_m1m2_overlap.summary',
    'motif_enrichment/elegans/repeat_GL_RE_overlap.summary',
    'gene_annotation/elegans.GL_specific_genes.m1m2_unique.genes',
    'gene_annotation/briggsae.genes_DESeq_glp1_vs_wt.txt',
    'data/external_data/boeck_unified_exp_data.txt'


#expand('species/{sample}/mitehunter/{sample}_mitehunter.mite.all.fa', sample=REPEAT_SEARCH_SPECIES),
#'gene_annotation/elegans.germline_expressed.serizay.bed'
