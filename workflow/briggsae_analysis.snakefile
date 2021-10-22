rule m1m2_enrichment_non_caenorhabditis:
  input:
    'data/external_data/elegans.{sample}.1to1.txt',
    'species/{sample}/gene_annotation/{sample}.gene_annotation.coding_gene_longest_isoform.bed',
    'species/{sample}/genome/{sample}.chrom.sizes.txt',
    'gene_annotation/elegans.GL_specific_genes.bed',
    'gene_annotation/elegans.soma_specific_genes.bed',
    'species/elegans/genome/elegans.chrom.sizes.txt',
    'motif_enrichment/elegans/elegans.m1m2_clusters.bed',
    'motif_enrichment/{sample}/{sample}.m1m2_clusters.bed',
  output:
    'species/{sample}/gene_annotation/{sample}.gene_annotation.coding_gene_longest_isoform.sorted.bed',
    temp('motif_conservation_strongylida/elegans.{sample}.ortho_promoters.bed'),
    temp('motif_conservation_strongylida/{sample}.elegans.ortho_promoters.bed'),
    'motif_conservation_strongylida/{sample}.motif_orthology.txt'
  shell:
    '''
    sort -k 4,4 {input[1]} > {output[0]}
    join -1 4 -2 1 {input[3]} {input[0]} | awk 'BEGIN{{OFS="\t";}}{{print $2, $3, $4, $1, $5, $6}}' | flankBed -l 500 -r 0 -s -i stdin -g {input[5]} | cut -f 1,2,3,4,5,6 | slopBed -l 0 -r 100 -s -i stdin -g {input[5]} > {output[1]}
    join -1 1 -2 4 {input[0]} {input[3]} | awk '{{print $2}}' | sort | join -1 1 -2 4 - {output[0]} | awk 'BEGIN{{OFS="\t";}}{{print $2, $3, $4, $1, $5, $6}}' | flankBed -l 500 -r 0 -s -i stdin -g {input[2]} | slopBed -l 0 -r 100 -s -i stdin -g {input[2]} > {output[2]}
    echo -e "orthologs\twith_m1m2" > {output[3]}
    ele_ortho_gl=$(wc -l < {output[1]})
    ele_ortho_m1m2_gl=$(intersectBed -a {output[1]} -b {input[6]} -u | wc -l)
    other_ortho_gl=$(wc -l < {output[2]})
    other_ortho_m1m2_gl=$(intersectBed -a {output[2]} -b {input[7]} -u | wc -l)
    join -1 4 -2 1 {input[4]} {input[0]} | awk 'BEGIN{{OFS="\t";}}{{print $2, $3, $4, $1, $5, $6}}' | flankBed -l 500 -r 0 -s -i stdin -g {input[5]} | cut -f 1,2,3,4,5,6 | slopBed -l 0 -r 100 -s -i stdin -g {input[5]} > {output[1]}
    join -1 1 -2 4 {input[0]} {input[4]} | awk '{{print $2}}' | sort | join -1 1 -2 4 - {output[0]} | awk 'BEGIN{{OFS="\t";}}{{print $2, $3, $4, $1, $5, $6}}' | flankBed -l 500 -r 0 -s -i stdin -g {input[2]} | slopBed -l 0 -r 100 -s -i stdin -g {input[2]} > {output[2]}
    ele_ortho_soma=$(wc -l < {output[1]})
    ele_ortho_m1m2_soma=$(intersectBed -a {output[1]} -b {input[6]} -u | wc -l)
    other_ortho_soma=$(wc -l < {output[2]})
    other_ortho_m1m2_soma=$(intersectBed -a {output[2]} -b {input[7]} -u | wc -l)
    echo -e "elegans_gl\t"$ele_ortho_gl"\t"$ele_ortho_m1m2_gl"\n"{wildcards.sample}"_gl\t"$other_ortho_gl"\t"$other_ortho_m1m2_gl"\nelegans_soma\t"$ele_ortho_soma"\t"$ele_ortho_m1m2_soma"\n"{wildcards.sample}"_soma\t"$other_ortho_soma"\t"$other_ortho_m1m2_soma >> {output[3]}
    '''


rule cb_gl_elements_annotation:
  input:
    "relmapping/annot_cb/reg_elements_cb.bed"
  output:
    parsed_RE = "RE_annotation/reg_elements_all.briggsae.bed",
    GL_RE_diffbind = "RE_annotation/reg_elements_all.briggsae.gl_specific.diffbind.bed",
    GL_RE = "RE_annotation/reg_elements_all.briggsae.gl_specific.bed",
    nonGL_RE = "RE_annotation/reg_elements_all.briggsae.not_gl_specific.bed"
  conda:
    "env/diffbind.yaml"
  shell:
    '''
    sed '1d' {input} | awk -F'[\t=;]' 'BEGIN{{OFS="\t";}}{{print $1, $2, $3, "ce.re" NR "." $7, $(NF-4), $(NF-3)}}' > {output.parsed_RE}
    Rscript scripts/diffbind_cb_glre.R
    intersectBed -a {output.parsed_RE} -b {output.GL_RE_diffbind} -u > {output.GL_RE}
    intersectBed -a {output.parsed_RE} -b {output.GL_RE_diffbind} -v > {output.nonGL_RE}
    '''


rule cb_gl_motif_enrichment:
  input:
    genome = "species/briggsae/genome/briggsae.fa",
    GL_RE = "RE_annotation/reg_elements_all.briggsae.gl_specific.bed",
    nonGL_RE = "RE_annotation/reg_elements_all.briggsae.not_gl_specific.bed"
  output:
    GL_RE_fasta = temp("RE_annotation/reg_elements_all.briggsae.gl_specific.fasta"),
    nonGL_RE_fasta = temp("RE_annotation/reg_elements_all.briggsae.not_gl_specific.fasta"),
    meme_results = directory("meme/briggsae.GL_specific_vs_nonGL_specific")
  resources:
    ntasks=6
  shell:
    '''
    fastaFromBed -fi {input.genome} -bed {input.GL_RE} -fo {output.GL_RE_fasta} -s
    fastaFromBed -fi {input.genome} -bed {input.nonGL_RE} -fo {output.nonGL_RE_fasta} -s
    meme-chip -oc {output.meme_results} -neg {output.nonGL_RE_fasta} -spamo-skip -fimo-skip -meme-p 6 -meme-nmotifs 6 -meme-minw 5 -meme-maxw 20 -dreme-m 0 {output.GL_RE_fasta}
    '''


rule cb_motif_RE_overlap:
  input:
    GL_RE = "RE_annotation/reg_elements_all.briggsae.gl_specific.bed",
    nonGL_RE = "RE_annotation/reg_elements_all.briggsae.not_gl_specific.bed",
    m1m2_pairs = "motif_enrichment/briggsae/briggsae.m1m2_clusters.bed"
  output:
    motif_re_summary = "RE_features/briggsae.motif_GL_overlap.summary",
    GL_RE_motif = "RE_annotation/reg_elements_all.briggsae.gl_specific.m1m2.bed",
    nonGL_RE_motif = "RE_annotation/reg_elements_all.briggsae.not_gl_specific.m1m2.bed"
  shell:
    '''
    echo -e "n_elements\twith_motif_pair\tconvergent\tdivergent\ttandem_m1m2\ttandem_m2m1\tn_promoters\twith_motif_pair\tconvergent\tdivergent\ttandem_m1m2\ttandem_m2m1" > {output.motif_re_summary}
    GL_n=$(wc -l < {input.GL_RE})
    notGL_n=$(wc -l < {input.nonGL_RE})
    intersectBed -a {input.GL_RE} -b {input.m1m2_pairs} -u > {output.GL_RE_motif}
    intersectBed -a {input.nonGL_RE} -b {input.m1m2_pairs} -u > {output.nonGL_RE_motif}
    GL_n_motif_all=$(wc -l < {output.GL_RE_motif})
    notGL_n_motif_all=$(wc -l < {output.nonGL_RE_motif})
    motif_GL_RE_type=$(intersectBed -a {input.m1m2_pairs} -b {input.GL_RE} -u | cut -f 5 | sort | uniq -c | awk '{{print $1}}' | paste -s -d '\t')
    motif_nonGL_RE_type=$(intersectBed -a {input.m1m2_pairs} -b {input.nonGL_RE} -u | cut -f 5 | sort | uniq -c | awk '{{print $1}}' | paste -s -d '\t')
    GL_prom_n=$(grep coding_promoter {input.GL_RE} | wc -l)
    notGL_prom_n=$(grep coding_promoter {input.nonGL_RE} | wc -l)
    GL_prom_motif_all=$(intersectBed -a {input.GL_RE} -b {input.m1m2_pairs} -u | grep coding_promoter | wc -l)
    notGL_prom_motif_all=$(intersectBed -a {input.nonGL_RE} -b {input.m1m2_pairs} -u | grep coding_promoter | wc -l)
    motif_GL_prom_type=$(grep coding_promoter {input.GL_RE} | intersectBed -a {input.m1m2_pairs} -b stdin -u | cut -f 5 | sort | uniq -c | awk '{{print $1}}' | paste -s -d '\t')
    motif_nonGL_prom_type=$(grep coding_promoter {input.nonGL_RE} | intersectBed -a {input.m1m2_pairs} -b stdin -u | cut -f 5 | sort | uniq -c | awk '{{print $1}}' | paste -s -d '\t')
    echo -e "GL_specific\t"$GL_n"\t"$GL_n_motif_all"\t"$motif_GL_RE_type"\t"$GL_prom_n"\t"$GL_prom_motif_all"\t"$motif_GL_prom_type >> {output.motif_re_summary}
    echo -e "not_GL_specific\t"$notGL_n"\t"$notGL_n_motif_all"\t"$motif_nonGL_RE_type"\t"$notGL_prom_n"\t"$notGL_prom_motif_all"\t"$motif_nonGL_prom_type >> {output.motif_re_summary}
    '''


rule cb_gl_m1m2_motif_enrichment:
  input:
    genome = "species/briggsae/genome/briggsae.fa",
    GL_RE = "RE_annotation/reg_elements_all.briggsae.gl_specific.bed",
    nonGL_RE = "RE_annotation/reg_elements_all.briggsae.not_gl_specific.bed",
    m1m2_pairs = "motif_enrichment/briggsae/briggsae.m1m2_clusters.bed"
  output:
    GL_RE_m1m2_fasta = temp("RE_annotation/reg_elements_all.briggsae.gl_specific.m1m2.fasta"),
    nonGL_RE_m1m2_fasta = temp("RE_annotation/reg_elements_all.briggsae.not_gl_specific.m1m2.fasta"),
    GL_RE_no_m1m2_fasta = temp("RE_annotation/reg_elements_all.briggsae.gl_specific.no_m1m2.fasta"),
    meme_results_1 = directory("meme/briggsae.GL_specific_m1m2_vs_nonGL_specific_m1m2"),
    meme_results_2 = directory("meme/briggsae.GL_specific_no_m1m2_vs_GL_specific_m1m2")
  resources:
    ntasks=6
  shell:
    '''
    intersectBed -a {input.GL_RE} -b {input.m1m2_pairs} -u | fastaFromBed -fi {input.genome} -bed stdin -fo {output.GL_RE_m1m2_fasta} -s
    intersectBed -a {input.nonGL_RE} -b {input.m1m2_pairs} -u | fastaFromBed -fi {input.genome} -bed stdin -fo {output.nonGL_RE_m1m2_fasta} -s
    intersectBed -a {input.GL_RE} -b {input.m1m2_pairs} -v | fastaFromBed -fi {input.genome} -bed stdin -fo {output.GL_RE_no_m1m2_fasta} -s
    meme -oc {output.meme_results_1} -objfun de -dna -revcomp -nmotifs 10 -p 6 -neg {output.nonGL_RE_m1m2_fasta} {output.GL_RE_m1m2_fasta}
    meme -oc {output.meme_results_2} -objfun de -dna -revcomp -nmotifs 10 -p 6 -neg {output.GL_RE_m1m2_fasta} {output.GL_RE_no_m1m2_fasta}
    '''


rule cb_RE_GC_content:
  input:
    GL_RE = "RE_annotation/reg_elements_all.briggsae.gl_specific.bed",
    nonGL_RE = "RE_annotation/reg_elements_all.briggsae.not_gl_specific.bed",
    m1m2_pairs = "motif_enrichment/briggsae/briggsae.m1m2_clusters.bed",
    genome = "species/briggsae/genome/briggsae.fa"
  output:
    GL_promoter_motif_nuc = "RE_features/reg_elements_all.briggsae.gl_specific.promoters.m1m2.nuc",
    nonGL_promoter_motif_nuc = "RE_features/reg_elements_all.briggsae.not_gl_specific.promoters.m1m2.nuc",
    GL_promoter_motif_fa = temp("RE_features/reg_elements_all.briggsae.gl_specific.promoters.m1m2.fa"),
    nonGL_promoter_motif_fa = temp("RE_features/reg_elements_all.briggsae.not_gl_specific.promoters.m1m2.fa"),
    GL_promoter_motif_cpg = "RE_features/reg_elements_all.briggsae.gl_specific.promoters.m1m2.cpg",
    nonGL_promoter_motif_cpg = "RE_features/reg_elements_all.briggsae.not_gl_specific.promoters.m1m2.cpg",
    GL_promoter_no_motif_nuc = "RE_features/reg_elements_all.briggsae.gl_specific.promoters.no_m1m2.nuc",
    nonGL_promoter_no_motif_nuc = "RE_features/reg_elements_all.briggsae.not_gl_specific.promoters.no_m1m2.nuc",
    GL_promoter_no_motif_fa = temp("RE_features/reg_elements_all.briggsae.gl_specific.promoters.no_m1m2.fa"),
    nonGL_promoter_no_motif_fa = temp("RE_features/reg_elements_all.briggsae.not_gl_specific.promoters.no_m1m2.fa"),
    GL_promoter_no_motif_cpg = "RE_features/reg_elements_all.briggsae.gl_specific.promoters.no_m1m2.cpg",
    nonGL_promoter_no_motif_cpg = "RE_features/reg_elements_all.briggsae.not_gl_specific.promoters.no_m1m2.cpg"
  shell:
    '''
    intersectBed -a {input.GL_RE} -b {input.m1m2_pairs} -u | grep coding_promoter | nucBed -fi {input.genome} -bed stdin -s > {output.GL_promoter_motif_nuc}
    intersectBed -a {input.nonGL_RE} -b {input.m1m2_pairs} -u | grep coding_promoter | nucBed -fi {input.genome} -bed stdin -s > {output.nonGL_promoter_motif_nuc}
    intersectBed -a {input.GL_RE} -b {input.m1m2_pairs} -u | grep coding_promoter | fastaFromBed -fi {input.genome} -bed stdin -s -fo {output.GL_promoter_motif_fa}
    intersectBed -a {input.nonGL_RE} -b {input.m1m2_pairs} -u | grep coding_promoter | fastaFromBed -fi {input.genome} -bed stdin -s -fo {output.nonGL_promoter_motif_fa}
    python scripts/CpG_dinucleotide_count.py {output.GL_promoter_motif_fa} {output.GL_promoter_motif_cpg}
    python scripts/CpG_dinucleotide_count.py {output.nonGL_promoter_motif_fa} {output.nonGL_promoter_motif_cpg}
    intersectBed -a {input.GL_RE} -b {input.m1m2_pairs} -v | grep coding_promoter | nucBed -fi {input.genome} -bed stdin -s > {output.GL_promoter_no_motif_nuc}
    intersectBed -a {input.nonGL_RE} -b {input.m1m2_pairs} -v | grep coding_promoter | nucBed -fi {input.genome} -bed stdin -s > {output.nonGL_promoter_no_motif_nuc}
    intersectBed -a {input.GL_RE} -b {input.m1m2_pairs} -v | grep coding_promoter | fastaFromBed -fi {input.genome} -bed stdin -s -fo {output.GL_promoter_no_motif_fa}
    intersectBed -a {input.nonGL_RE} -b {input.m1m2_pairs} -v | grep coding_promoter | fastaFromBed -fi {input.genome} -bed stdin -s -fo {output.nonGL_promoter_no_motif_fa}
    python scripts/CpG_dinucleotide_count.py {output.GL_promoter_no_motif_fa} {output.GL_promoter_no_motif_cpg}
    python scripts/CpG_dinucleotide_count.py {output.nonGL_promoter_no_motif_fa} {output.nonGL_promoter_no_motif_cpg}
    '''


rule cb_motif_associated_genes:
  input:
    annot_cb = "relmapping/annot_cb/reg_elements_cb.bed",
    re_all = "RE_annotation/reg_elements_all.briggsae.bed",
    GL_RE = "RE_annotation/reg_elements_all.briggsae.gl_specific.bed",
    nonGL_RE = "RE_annotation/reg_elements_all.briggsae.not_gl_specific.bed",
    m1m2_pairs = "motif_enrichment/briggsae/briggsae.m1m2_clusters.bed"
  output:
    promoters_genes = "promoter_annotation/promoters_all.briggsae.bed",
    m1m2_GL_genes = "promoter_annotation/promoters_gl_specific.briggsae.m1m2.any_promoter.genes",
    no_m1m2_GL_genes = "promoter_annotation/promoters_gl_specific.briggsae.no_m1m2.any_promoter.genes",
    m1m2_noGL_genes = "promoter_annotation/promoters_not_gl_specific.briggsae.m1m2.any_promoter.genes",
    no_m1m2_noGL_genes = "promoter_annotation/promoters_not_gl_specific.briggsae.no_m1m2.any_promoter.genes",
    unique_promoters_genes = temp("promoter_annotation/promoters_all.briggsae.unique.genes"),
    unique_promoters = "promoter_annotation/promoters_all.briggsae.unique.bed",
    m1m2_GL_genes_uniq = "promoter_annotation/promoters_gl_specific.briggsae.m1m2.unique_promoter.genes",
    no_m1m2_GL_genes_uniq = "promoter_annotation/promoters_gl_specific.briggsae.no_m1m2.unique_promoter.genes",
    m1m2_noGL_genes_uniq = "promoter_annotation/promoters_not_gl_specific.briggsae.m1m2.unique_promoter.genes",
    no_m1m2_noGL_genes_uniq = "promoter_annotation/promoters_not_gl_specific.briggsae.no_m1m2.unique_promoter.genes"
  shell:
    '''
    grep "annot=coding_promoter" {input.annot_cb} | awk 'BEGIN{{OFS="\t";FS="\t|=|;";}}{{if (length($5) == 35) print $1, $2, $3, substr($5, 0, 14) "," substr($5, 22, 35); else print $1, $2, $3, $5}}' | intersectBed -a stdin -b {input.re_all} -wa -wb | cut -f 1,2,3,4,9,10 > {output.promoters_genes}
    intersectBed -a {output.promoters_genes} -b {input.GL_RE} -u | intersectBed -a stdin -b {input.m1m2_pairs} -u | cut -f 4 | tr "," "\n" | sort > {output.m1m2_GL_genes}
    intersectBed -a {output.promoters_genes} -b {input.GL_RE} -u | intersectBed -a stdin -b {input.m1m2_pairs} -v | cut -f 4 | tr "," "\n" | sort | join -v 1 - {output.m1m2_GL_genes} > {output.no_m1m2_GL_genes}
    intersectBed -a {output.promoters_genes} -b {input.GL_RE} -v | intersectBed -a stdin -b {input.m1m2_pairs} -u | cut -f 4 | tr "," "\n" | sort > {output.m1m2_noGL_genes}
    intersectBed -a {output.promoters_genes} -b {input.GL_RE} -v | intersectBed -a stdin -b {input.m1m2_pairs} -v | cut -f 4 | tr "," "\n" | sort | join -v 1 - {output.m1m2_noGL_genes} > {output.no_m1m2_noGL_genes}
    awk 'BEGIN{{OFS="\t";}}{{if ($4 ~ ",") print $1, $2, $3, substr($4, 0, 14), $5, $6 "\\n" $1, $2, $3, substr($4, 16, 29), $5, $6; else print $0}}' {output.promoters_genes} | cut -f 4 | sort | uniq -c | awk '{{if ($1 == 1) print $2}}' > {output.unique_promoters_genes}
    awk 'BEGIN{{OFS="\t";}}{{if ($4 ~ ",") print $1, $2, $3, substr($4, 0, 14), $5, $6 "\\n" $1, $2, $3, substr($4, 16, 29), $5, $6; else print $0}}' {output.promoters_genes} | sort -k 4,4 | join -1 4 -2 1 - {output.unique_promoters_genes} | awk 'BEGIN{{OFS="\t";}}{{print $2, $3, $4, $1, $5, $6}}' > {output.unique_promoters}
    intersectBed -a {output.unique_promoters} -b {input.GL_RE} -u | intersectBed -a stdin -b {input.m1m2_pairs} -u | cut -f 4 | tr "," "\n" | sort > {output.m1m2_GL_genes_uniq}
    intersectBed -a {output.unique_promoters} -b {input.GL_RE} -u | intersectBed -a stdin -b {input.m1m2_pairs} -v | cut -f 4 | tr "," "\n" | sort > {output.no_m1m2_GL_genes_uniq}
    intersectBed -a {output.unique_promoters} -b {input.GL_RE} -v | intersectBed -a stdin -b {input.m1m2_pairs} -u | cut -f 4 | tr "," "\n" | sort > {output.m1m2_noGL_genes_uniq}
    intersectBed -a {output.unique_promoters} -b {input.GL_RE} -v | intersectBed -a stdin -b {input.m1m2_pairs} -v | cut -f 4 | tr "," "\n" | sort > {output.no_m1m2_noGL_genes_uniq}
    '''


rule GLRE_conservation:
  input:
    "promoter_annotation/promoters_all.elegans.bed",
    "promoter_annotation/promoters_all.briggsae.bed",
    "RE_annotation/reg_elements_all.elegans.gl_specific.m1m2.bed",
    "RE_annotation/reg_elements_all.briggsae.gl_specific.m1m2.bed",
    "species/elegans/orthologs/one2one.briggsae.orthologs.txt",
    "motif_enrichment/briggsae/briggsae.m1m2_clusters.bed",
    "motif_enrichment/briggsae/briggsae.m1.bed",
    "motif_enrichment/briggsae/briggsae.m2.bed",
    "motif_enrichment/elegans/elegans.m1m2_clusters.bed",
    "motif_enrichment/elegans/elegans.m1.bed",
    "motif_enrichment/elegans/elegans.m2.bed",
    "species/briggsae/gene_annotation/briggsae.gene_annotation.coding_gene_longest_isoform.bed",
    "species/briggsae/genome/briggsae.chrom.sizes.txt",
    "species/elegans/gene_annotation/elegans.gene_annotation.coding_gene_longest_isoform.bed",
    "species/elegans/genome/elegans.chrom.sizes.txt",
    "motif_enrichment/briggsae/briggsae.m1.permissive.bed",
    "motif_enrichment/briggsae/briggsae.m2.permissive.bed",
    "motif_enrichment/elegans/elegans.m1.permissive.bed",
    "motif_enrichment/elegans/elegans.m2.permissive.bed",
  output:
    "RE_conservation/promoters_all.elegans.associated_genes.bed",
    "RE_conservation/promoters_all.briggsae.associated_genes.bed",
    "RE_conservation/promoters_GL_m1m2.elegans.associated_genes.bed",
    "RE_conservation/promoters_GL_m1m2.briggsae.associated_genes.bed",
    "RE_conservation/elegans_to_briggsae/promoters_GL_m1m2.elegans.associated_genes.elegans_briggsae_orthologs",
    "RE_conservation/elegans_to_briggsae/briggsae.to_elegans_orthologs.promoters_all.bed",
    "RE_conservation/elegans_to_briggsae/briggsae.to_elegans_orthologs.promoters_m1m2.genes",
    "RE_conservation/elegans_to_briggsae/briggsae.to_elegans_orthologs.promoters_m1_dm2.genes",
    "RE_conservation/elegans_to_briggsae/briggsae.to_elegans_orthologs.promoters_m1_only.genes",
    "RE_conservation/elegans_to_briggsae/briggsae.to_elegans_orthologs.promoters_m2_dm1.genes",
    "RE_conservation/elegans_to_briggsae/briggsae.to_elegans_orthologs.promoters_m2_only.genes",
    "RE_conservation/elegans_to_briggsae/briggsae.to_elegans_orthologs.promoters_m1_m2_no_pair.genes",
    "RE_conservation/elegans_to_briggsae/briggsae.to_elegans_orthologs.promoters_no_motifs.genes",
    "RE_conservation/conservation_statistics.txt",
    "RE_conservation/briggsae_to_elegans/promoters_GL_m1m2.briggsae.associated_genes.briggsae_elegans_orthologs",
    "RE_conservation/briggsae_to_elegans/elegans.to_briggsae_orthologs.promoters_all.bed",
    "RE_conservation/briggsae_to_elegans/elegans.to_briggsae_orthologs.promoters_m1m2.genes",
    "RE_conservation/briggsae_to_elegans/elegans.to_briggsae_orthologs.promoters_m1_dm2.genes",
    "RE_conservation/briggsae_to_elegans/elegans.to_briggsae_orthologs.promoters_m1_only.genes",
    "RE_conservation/briggsae_to_elegans/elegans.to_briggsae_orthologs.promoters_m2_dm1.genes",
    "RE_conservation/briggsae_to_elegans/elegans.to_briggsae_orthologs.promoters_m2_only.genes",
    "RE_conservation/briggsae_to_elegans/elegans.to_briggsae_orthologs.promoters_m1_m2_no_pair.genes",
    "RE_conservation/briggsae_to_elegans/elegans.to_briggsae_orthologs.promoters_no_motifs.genes",
    "RE_conservation/elegans_to_briggsae/briggsae.to_elegans_orthologs.promoters_no_motifs.flanking_m1m2.genes",
    "RE_conservation/briggsae_to_elegans/elegans.to_briggsae_orthologs.promoters_no_motifs.flanking_m1m2.genes",
  shell:
    '''
    awk 'BEGIN{{OFS="\t";}}{{if ($4 ~ ",") print $1, $2, $3, substr($4, 1, 14), $5, $6 "\\n" $1, $2, $3, substr($4, 16, 29), $5, $6; else print $0}}' {input[0]} | sort -k 4,4 > {output[0]}
    awk 'BEGIN{{OFS="\t";}}{{if ($4 ~ ",") print $1, $2, $3, substr($4, 1, 14), $5, $6 "\\n" $1, $2, $3, substr($4, 16, 29), $5, $6; else print $0}}' {input[1]} | sort -k 4,4 > {output[1]}
    intersectBed -a {output[0]} -b {input[2]} -u > {output[2]}
    intersectBed -a {output[1]} -b {input[3]} -u > {output[3]}
    sort -k 1,1 {input[4]} | join -1 1 -2 4 - {output[2]} | awk '{{print $1 "\t" $2}}' | sort -k 2,2 > {output[4]}
    join -1 4 -2 2 {output[1]} {output[4]} | awk 'BEGIN{{OFS="\t";}}{{print $2, $3, $4, $1, $5, $6}}' > {output[5]}
    intersectBed -a {output[5]} -b {input[5]} -u | cut -f 4 | sort | uniq > {output[6]}
    intersectBed -a {output[5]} -b {input[6]} -u | intersectBed -a stdin -b {input[16]} -u | intersectBed -a stdin -b {input[7]} -v | cut -f 4 | sort | uniq | join -v 1 - {output[6]} > {output[7]}
    intersectBed -a {output[5]} -b {input[6]} -u | intersectBed -a stdin -b {input[16]} -v | cut -f 4 | sort | uniq | join -v 1 - {output[6]} > {output[8]}
    intersectBed -a {output[5]} -b {input[7]} -u | intersectBed -a stdin -b {input[15]} -u | intersectBed -a stdin -b {input[6]} -v | cut -f 4 | sort | uniq | join -v 1 - {output[6]} > {output[9]}
    intersectBed -a {output[5]} -b {input[7]} -u | intersectBed -a stdin -b {input[15]} -v | cut -f 4 | sort | uniq | join -v 1 - {output[6]} > {output[10]}
    intersectBed -a {output[5]} -b {input[6]} -u | intersectBed -a stdin -b {input[7]} -u | cut -f 4 | sort | uniq | join -v 1 - {output[6]} > {output[11]}
    intersectBed -a {output[5]} -b {input[6]} -v | intersectBed -a stdin -b {input[7]} -v | cut -f 4 | sort | uniq | join -v 1 - {output[6]} | join -v 1 - {output[7]} | join -v 1 - {output[8]} | join -v 1 - {output[9]} | join -v 1 - {output[10]} | join -v 1 - {output[11]} > {output[12]}
    sort -k 4,4 {input[11]} | join -1 1 -2 4 {output[12]} - | awk 'BEGIN{{OFS="\t";}}{{print $2, $3, $4, $1, $5, $6}}' | flankBed -l 1000 -r 0 -s -i stdin -g {input[12]} | slopBed -r 200 -l 0 -s -i stdin -g {input[12]} | intersectBed -a stdin -b {input[5]} -u | cut -f 4 | sort > {output[23]}
    echo -e "species1_GL_m1m2_genes\twith_orthologs\tspecies2_orthologs_with_promoters\tm1m2_orthologs\tm1_dm2_orthologs\tm1_only_orthologs\tm2_dm1_orthologs\tm2_only_orthologs\tm1m2_nopair_orthologs\tno_motif_orthologs_permissive\tno_motif_orthologs_flanking_m1m2" > {output[13]}
    species1_GL_m1m2_genes=$(cut -f 4 {output[2]} | sort | uniq | wc -l)
    with_orthologs=$(wc -l < {output[4]})
    species2_orthologs_with_promoters=$(cut -f 4 {output[5]} | sort | uniq | wc -l)
    m1m2_orthologs=$(wc -l < {output[6]})
    m1_dm2_orthologs=$(wc -l < {output[7]})
    m1_only_orthologs=$(wc -l < {output[8]})
    m2_dm1_orthologs=$(wc -l < {output[9]})
    m2_only_orthologs=$(wc -l < {output[10]})
    m1m2_nopair_orthologs=$(wc -l < {output[11]})
    no_motif_orthologs=$(wc -l < {output[12]})
    no_motif_orthologs_stringent=$(wc -l < {output[23]})
    echo -e "elegans.briggsae\t"$species1_GL_m1m2_genes"\t"$with_orthologs"\t"$species2_orthologs_with_promoters"\t"$m1m2_orthologs"\t"$m1_dm2_orthologs"\t"$m1_only_orthologs"\t"$m2_dm1_orthologs"\t"$m2_only_orthologs"\t"$m1m2_nopair_orthologs"\t"$no_motif_orthologs"\t"$no_motif_orthologs_stringent >> {output[13]}
    sort -k 2,2 {input[4]} | join -1 2 -2 4 - {output[3]} | awk '{{print $1 "\t" $2}}' | sort -k 2,2 > {output[14]}
    join -1 4 -2 2 {output[0]} {output[14]} | awk 'BEGIN{{OFS="\t";}}{{print $2, $3, $4, $1, $5, $6}}' > {output[15]}
    intersectBed -a {output[15]} -b {input[8]} -u | cut -f 4 | sort | uniq > {output[16]}
    intersectBed -a {output[15]} -b {input[9]} -u | intersectBed -a stdin -b {input[18]} -u | intersectBed -a stdin -b {input[10]} -v | cut -f 4 | sort | uniq | join -v 1 - {output[16]} > {output[17]}
    intersectBed -a {output[15]} -b {input[9]} -u | intersectBed -a stdin -b {input[18]} -v | cut -f 4 | sort | uniq | join -v 1 - {output[16]} > {output[18]}
    intersectBed -a {output[15]} -b {input[10]} -u | intersectBed -a stdin -b {input[17]} -u | intersectBed -a stdin -b {input[9]} -v | cut -f 4 | sort | uniq | join -v 1 - {output[16]} > {output[19]}
    intersectBed -a {output[15]} -b {input[10]} -u | intersectBed -a stdin -b {input[17]} -v | cut -f 4 | sort | uniq | join -v 1 - {output[16]} > {output[20]}
    intersectBed -a {output[15]} -b {input[9]} -u | intersectBed -a stdin -b {input[10]} -u | cut -f 4 | sort | uniq | join -v 1 - {output[16]} > {output[21]}    
    intersectBed -a {output[15]} -b {input[9]} -v | intersectBed -a stdin -b {input[10]} -v | cut -f 4 | sort | uniq | join -v 1 - {output[16]} | join -v 1 - {output[17]} | join -v 1 - {output[18]} | join -v 1 - {output[19]} | join -v 1 - {output[20]} | join -v 1 - {output[21]} > {output[22]}
    sort -k 4,4 {input[13]} | join -1 1 -2 4 {output[22]} - | awk 'BEGIN{{OFS="\t";}}{{print $2, $3, $4, $1, $5, $6}}' | flankBed -l 1000 -r 0 -s -i stdin -g {input[14]} | slopBed -r 200 -l 0 -s -i stdin -g {input[14]} | intersectBed -a stdin -b {input[8]} -u | cut -f 4 | sort > {output[24]}
    species1_GL_m1m2_genes=$(cut -f 4 {output[3]} | sort | uniq | wc -l)
    with_orthologs=$(wc -l < {output[14]})
    species2_orthologs_with_promoters=$(cut -f 4 {output[15]} | sort | uniq | wc -l)
    m1m2_orthologs=$(wc -l < {output[16]})
    m1_dm2_orthologs=$(wc -l < {output[17]})
    m1_only_orthologs=$(wc -l < {output[18]})
    m2_dm1_orthologs=$(wc -l < {output[19]})
    m2_only_orthologs=$(wc -l < {output[20]})
    m1m2_nopair_orthologs=$(wc -l < {output[21]})
    no_motif_orthologs=$(wc -l < {output[22]})
    no_motif_orthologs_stringent=$(wc -l < {output[24]})
    echo -e "briggsae.elegans\t"$species1_GL_m1m2_genes"\t"$with_orthologs"\t"$species2_orthologs_with_promoters"\t"$m1m2_orthologs"\t"$m1_dm2_orthologs"\t"$m1_only_orthologs"\t"$m2_dm1_orthologs"\t"$m2_only_orthologs"\t"$m1m2_nopair_orthologs"\t"$no_motif_orthologs"\t"$no_motif_orthologs_stringent >> {output[13]}
    '''

rule RE_conservation:
  input:
    "promoter_annotation/promoters_all.elegans.bed",
    "promoter_annotation/promoters_all.briggsae.bed",
    "RE_annotation/reg_elements_all.elegans.bed",
    "RE_annotation/reg_elements_all.briggsae.bed",
    "species/elegans/orthologs/one2one.briggsae.orthologs.txt",
    "motif_enrichment/briggsae/briggsae.m1m2_clusters.bed",
    "motif_enrichment/briggsae/briggsae.m1.bed",
    "motif_enrichment/briggsae/briggsae.m2.bed",
    "motif_enrichment/elegans/elegans.m1m2_clusters.bed",
    "motif_enrichment/elegans/elegans.m1.bed",
    "motif_enrichment/elegans/elegans.m2.bed",
    "species/briggsae/gene_annotation/briggsae.gene_annotation.coding_gene_longest_isoform.bed",
    "species/briggsae/genome/briggsae.chrom.sizes.txt",
    "species/elegans/gene_annotation/elegans.gene_annotation.coding_gene_longest_isoform.bed",
    "species/elegans/genome/elegans.chrom.sizes.txt",
    "motif_enrichment/briggsae/briggsae.m1.permissive.bed",
    "motif_enrichment/briggsae/briggsae.m2.permissive.bed",
    "motif_enrichment/elegans/elegans.m1.permissive.bed",
    "motif_enrichment/elegans/elegans.m2.permissive.bed",
  output:
    "RE_conservation_all/promoters_all.elegans.associated_genes.bed",
    "RE_conservation_all/promoters_all.briggsae.associated_genes.bed",
    "RE_conservation_all/promoters_m1m2.elegans.associated_genes.bed",
    "RE_conservation_all/promoters_m1m2.briggsae.associated_genes.bed",
    "RE_conservation_all/elegans_to_briggsae/promoters_m1m2.elegans.associated_genes.elegans_briggsae_orthologs",
    "RE_conservation_all/elegans_to_briggsae/briggsae.to_elegans_orthologs.promoters_all.bed",
    "RE_conservation_all/elegans_to_briggsae/briggsae.to_elegans_orthologs.promoters_m1m2.genes",
    "RE_conservation_all/elegans_to_briggsae/briggsae.to_elegans_orthologs.promoters_m1_dm2.genes",
    "RE_conservation_all/elegans_to_briggsae/briggsae.to_elegans_orthologs.promoters_m1_only.genes",
    "RE_conservation_all/elegans_to_briggsae/briggsae.to_elegans_orthologs.promoters_m2_dm1.genes",
    "RE_conservation_all/elegans_to_briggsae/briggsae.to_elegans_orthologs.promoters_m2_only.genes",
    "RE_conservation_all/elegans_to_briggsae/briggsae.to_elegans_orthologs.promoters_m1_m2_no_pair.genes",
    "RE_conservation_all/elegans_to_briggsae/briggsae.to_elegans_orthologs.promoters_no_motifs.genes",
    "RE_conservation_all/conservation_statistics.txt",
    "RE_conservation_all/briggsae_to_elegans/promoters_m1m2.briggsae.associated_genes.briggsae_elegans_orthologs",
    "RE_conservation_all/briggsae_to_elegans/elegans.to_briggsae_orthologs.promoters_all.bed",
    "RE_conservation_all/briggsae_to_elegans/elegans.to_briggsae_orthologs.promoters_m1m2.genes",
    "RE_conservation_all/briggsae_to_elegans/elegans.to_briggsae_orthologs.promoters_m1_dm2.genes",
    "RE_conservation_all/briggsae_to_elegans/elegans.to_briggsae_orthologs.promoters_m1_only.genes",
    "RE_conservation_all/briggsae_to_elegans/elegans.to_briggsae_orthologs.promoters_m2_dm1.genes",
    "RE_conservation_all/briggsae_to_elegans/elegans.to_briggsae_orthologs.promoters_m2_only.genes",
    "RE_conservation_all/briggsae_to_elegans/elegans.to_briggsae_orthologs.promoters_m1_m2_no_pair.genes",
    "RE_conservation_all/briggsae_to_elegans/elegans.to_briggsae_orthologs.promoters_no_motifs.genes",
    "RE_conservation_all/elegans_to_briggsae/briggsae.to_elegans_orthologs.promoters_no_motifs.flanking_m1m2.genes",
    "RE_conservation_all/briggsae_to_elegans/elegans.to_briggsae_orthologs.promoters_no_motifs.flanking_m1m2.genes",
    "RE_conservation_all/conserved_vs_species_specific_promoters.m1m2_statistics.txt",
    "plots/conserved_sp_specific_promoters.m1m2_overlap.pdf",
  shell:
    '''
    awk 'BEGIN{{OFS="\t";}}{{if ($4 ~ ",") print $1, $2, $3, substr($4, 1, 14), $5, $6 "\\n" $1, $2, $3, substr($4, 16, 29), $5, $6; else print $0}}' {input[0]} | sort -k 4,4 > {output[0]}
    awk 'BEGIN{{OFS="\t";}}{{if ($4 ~ ",") print $1, $2, $3, substr($4, 1, 14), $5, $6 "\\n" $1, $2, $3, substr($4, 16, 29), $5, $6; else print $0}}' {input[1]} | sort -k 4,4 > {output[1]}
    intersectBed -a {output[0]} -b {input[2]} -u | intersectBed -a stdin -b {input[8]} -u > {output[2]}
    intersectBed -a {output[1]} -b {input[3]} -u | intersectBed -a stdin -b {input[5]} -u > {output[3]}
    sort -k 1,1 {input[4]} | join -1 1 -2 4 - {output[2]} | awk '{{print $1 "\t" $2}}' | sort -k 2,2 > {output[4]}
    join -1 4 -2 2 {output[1]} {output[4]} | awk 'BEGIN{{OFS="\t";}}{{print $2, $3, $4, $1, $5, $6}}' > {output[5]}
    intersectBed -a {output[5]} -b {input[5]} -u | cut -f 4 | sort | uniq > {output[6]}
    intersectBed -a {output[5]} -b {input[6]} -u | intersectBed -a stdin -b {input[16]} -u | intersectBed -a stdin -b {input[7]} -v | cut -f 4 | sort | uniq | join -v 1 - {output[6]} > {output[7]}
    intersectBed -a {output[5]} -b {input[6]} -u | intersectBed -a stdin -b {input[16]} -v | cut -f 4 | sort | uniq | join -v 1 - {output[6]} > {output[8]}
    intersectBed -a {output[5]} -b {input[7]} -u | intersectBed -a stdin -b {input[15]} -u | intersectBed -a stdin -b {input[6]} -v | cut -f 4 | sort | uniq | join -v 1 - {output[6]} > {output[9]}
    intersectBed -a {output[5]} -b {input[7]} -u | intersectBed -a stdin -b {input[15]} -v | cut -f 4 | sort | uniq | join -v 1 - {output[6]} > {output[10]}
    intersectBed -a {output[5]} -b {input[6]} -u | intersectBed -a stdin -b {input[7]} -u | cut -f 4 | sort | uniq | join -v 1 - {output[6]} > {output[11]}
    intersectBed -a {output[5]} -b {input[6]} -v | intersectBed -a stdin -b {input[7]} -v | cut -f 4 | sort | uniq | join -v 1 - {output[6]} | join -v 1 - {output[7]} | join -v 1 - {output[8]} | join -v 1 - {output[9]} | join -v 1 - {output[10]} | join -v 1 - {output[11]} > {output[12]}
    sort -k 4,4 {input[11]} | join -1 1 -2 4 {output[12]} - | awk 'BEGIN{{OFS="\t";}}{{print $2, $3, $4, $1, $5, $6}}' | flankBed -l 1000 -r 0 -s -i stdin -g {input[12]} | slopBed -r 200 -l 0 -s -i stdin -g {input[12]} | intersectBed -a stdin -b {input[5]} -u | cut -f 4 | sort > {output[23]}
    echo -e "species1_GL_m1m2_genes\twith_orthologs\tspecies2_orthologs_with_promoters\tm1m2_orthologs\tm1_dm2_orthologs\tm1_only_orthologs\tm2_dm1_orthologs\tm2_only_orthologs\tm1m2_nopair_orthologs\tno_motif_orthologs_permissive\tno_motif_orthologs_flanking_m1m2" > {output[13]}
    species1_GL_m1m2_genes=$(cut -f 4 {output[2]} | sort | uniq | wc -l)
    with_orthologs=$(wc -l < {output[4]})
    species2_orthologs_with_promoters=$(cut -f 4 {output[5]} | sort | uniq | wc -l)
    m1m2_orthologs=$(wc -l < {output[6]})
    m1_dm2_orthologs=$(wc -l < {output[7]})
    m1_only_orthologs=$(wc -l < {output[8]})
    m2_dm1_orthologs=$(wc -l < {output[9]})
    m2_only_orthologs=$(wc -l < {output[10]})
    m1m2_nopair_orthologs=$(wc -l < {output[11]})
    no_motif_orthologs=$(wc -l < {output[12]})
    no_motif_orthologs_stringent=$(wc -l < {output[23]})
    echo -e "elegans.briggsae\t"$species1_GL_m1m2_genes"\t"$with_orthologs"\t"$species2_orthologs_with_promoters"\t"$m1m2_orthologs"\t"$m1_dm2_orthologs"\t"$m1_only_orthologs"\t"$m2_dm1_orthologs"\t"$m2_only_orthologs"\t"$m1m2_nopair_orthologs"\t"$no_motif_orthologs"\t"$no_motif_orthologs_stringent >> {output[13]}
    sort -k 2,2 {input[4]} | join -1 2 -2 4 - {output[3]} | awk '{{print $1 "\t" $2}}' | sort -k 2,2 > {output[14]}
    join -1 4 -2 2 {output[0]} {output[14]} | awk 'BEGIN{{OFS="\t";}}{{print $2, $3, $4, $1, $5, $6}}' > {output[15]}
    intersectBed -a {output[15]} -b {input[8]} -u | cut -f 4 | sort | uniq > {output[16]}
    intersectBed -a {output[15]} -b {input[9]} -u | intersectBed -a stdin -b {input[18]} -u | intersectBed -a stdin -b {input[10]} -v | cut -f 4 | sort | uniq | join -v 1 - {output[16]} > {output[17]}
    intersectBed -a {output[15]} -b {input[9]} -u | intersectBed -a stdin -b {input[18]} -v | cut -f 4 | sort | uniq | join -v 1 - {output[16]} > {output[18]}
    intersectBed -a {output[15]} -b {input[10]} -u | intersectBed -a stdin -b {input[17]} -u | intersectBed -a stdin -b {input[9]} -v | cut -f 4 | sort | uniq | join -v 1 - {output[16]} > {output[19]}
    intersectBed -a {output[15]} -b {input[10]} -u | intersectBed -a stdin -b {input[17]} -v | cut -f 4 | sort | uniq | join -v 1 - {output[16]} > {output[20]}
    intersectBed -a {output[15]} -b {input[9]} -u | intersectBed -a stdin -b {input[10]} -u | cut -f 4 | sort | uniq | join -v 1 - {output[16]} > {output[21]}
    intersectBed -a {output[15]} -b {input[9]} -v | intersectBed -a stdin -b {input[10]} -v | cut -f 4 | sort | uniq | join -v 1 - {output[16]} | join -v 1 - {output[17]} | join -v 1 - {output[18]} | join -v 1 - {output[19]} | join -v 1 - {output[20]} | join -v 1 - {output[21]} > {output[22]}
    sort -k 4,4 {input[13]} | join -1 1 -2 4 {output[22]} - | awk 'BEGIN{{OFS="\t";}}{{print $2, $3, $4, $1, $5, $6}}' | flankBed -l 1000 -r 0 -s -i stdin -g {input[14]} | slopBed -r 200 -l 0 -s -i stdin -g {input[14]} | intersectBed -a stdin -b {input[8]} -u | cut -f 4 | sort > {output[24]}
    species1_GL_m1m2_genes=$(cut -f 4 {output[3]} | sort | uniq | wc -l)
    with_orthologs=$(wc -l < {output[14]})
    species2_orthologs_with_promoters=$(cut -f 4 {output[15]} | sort | uniq | wc -l)
    m1m2_orthologs=$(wc -l < {output[16]})
    m1_dm2_orthologs=$(wc -l < {output[17]})
    m1_only_orthologs=$(wc -l < {output[18]})
    m2_dm1_orthologs=$(wc -l < {output[19]})
    m2_only_orthologs=$(wc -l < {output[20]})
    m1m2_nopair_orthologs=$(wc -l < {output[21]})
    no_motif_orthologs=$(wc -l < {output[22]})
    no_motif_orthologs_stringent=$(wc -l < {output[24]})
    echo -e "briggsae.elegans\t"$species1_GL_m1m2_genes"\t"$with_orthologs"\t"$species2_orthologs_with_promoters"\t"$m1m2_orthologs"\t"$m1_dm2_orthologs"\t"$m1_only_orthologs"\t"$m2_dm1_orthologs"\t"$m2_only_orthologs"\t"$m1m2_nopair_orthologs"\t"$no_motif_orthologs"\t"$no_motif_orthologs_stringent >> {output[13]}
    echo -e "convergent\tdivergent\ttandem_m1m2\ttandem_m2m1" > {output[25]}
    elegans_conserved=$(cat {output[6]} {output[7]} {output[8]} {output[9]} {output[10]} {output[11]} {output[23]} | sort | uniq | join -1 1 -2 2 - {output[4]} | awk '{{print $2}}' | sort | join -1 1 -2 4 - {output[2]} | awk 'BEGIN{{OFS="\t";}}{{print $2, $3, $4, $1, $5, $6}}' | intersectBed -a {input[8]} -b stdin -u | cut -f 5 | sort | uniq -c | awk '{{printf "\t" $1}}')
    elegans_specific=$(join -v 1 {output[12]} {output[23]} | join -1 1 -2 2  - {output[4]} | awk '{{print $2}}' | sort | join -1 1 -2 4 - {output[2]} | awk 'BEGIN{{OFS="\t";}}{{print $2, $3, $4, $1, $5, $6}}' | intersectBed -a {input[8]} -b stdin -u | cut -f 5 | sort | uniq -c | awk '{{printf "\t" $1}}')
    briggsae_conserved=$(cat {output[16]} {output[17]} {output[18]} {output[19]} {output[20]} {output[21]} {output[24]} | sort | uniq | join -1 1 -2 2 - {output[14]} | awk '{{print $2}}' | sort | join -1 1 -2 4 - {output[3]} | awk 'BEGIN{{OFS="\t";}}{{print $2, $3, $4, $1, $5, $6}}' | intersectBed -a {input[5]} -b stdin -u | cut -f 5 | sort | uniq -c | awk '{{printf "\t" $1}}')
    briggsae_specific=$(join -v 1 {output[22]} {output[24]} | join -1 1 -2 2 - {output[14]} | awk '{{print $2}}' | sort | join -1 1 -2 4 - {output[3]} | awk 'BEGIN{{OFS="\t";}}{{print $2, $3, $4, $1, $5, $6}}' | intersectBed -a {input[5]} -b stdin -u | cut -f 5 | sort | uniq -c | awk '{{printf "\t" $1}}')
    echo -e "elegans_conserved"$elegans_conserved"\nelegans_specific"$elegans_specific"\nbriggsae_conserved"$briggsae_conserved"\nbriggsae_specific"$briggsae_specific"\n" >> {output[25]}
    Rscript scripts/RE_conservation_m1m2_fraction.R
    '''


rule pgc_rnaseq_download:
  output:
    'data/external_data/rnaseq/fastq/sL1_pgc_rep1.r1.fq.gz',
    'data/external_data/rnaseq/fastq/sL1_pgc_rep2.r1.fq.gz',
    'data/external_data/rnaseq/fastq/sL1_pgc_rep3.r1.fq.gz',
    'data/external_data/rnaseq/fastq/sL1_pgc_rep4.r1.fq.gz',
    'data/external_data/rnaseq/fastq/fL1_pgc_rep1.r1.fq.gz',
    'data/external_data/rnaseq/fastq/fL1_pgc_rep2.r1.fq.gz',
    'data/external_data/rnaseq/fastq/fL1_pgc_rep3.r1.fq.gz',
  shell:
    '''
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR577/005/SRR5772155/SRR5772155.fastq.gz -O {output[0]}
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR577/006/SRR5772156/SRR5772156.fastq.gz -O {output[1]}
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR577/007/SRR5772157/SRR5772157.fastq.gz -O {output[2]}
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR577/008/SRR5772158/SRR5772158.fastq.gz -O {output[3]}
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR577/001/SRR5772161/SRR5772161.fastq.gz -O {output[4]}
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR577/002/SRR5772162/SRR5772162.fastq.gz -O {output[5]}
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR577/003/SRR5772163/SRR5772163.fastq.gz -O {output[6]}
    '''


rule whole_worm_rnaseq_download:
  output:
    'data/external_data/rnaseq/fastq/elegans_L1_rep1.r1.fq.gz',
    temp('data/external_data/rnaseq/fastq/elegans_L1_rep2.A.r1.fq.gz'),
    temp('data/external_data/rnaseq/fastq/elegans_L1_rep2.B.r1.fq.gz'),
    'data/external_data/rnaseq/fastq/elegans_L1_rep2.r1.fq.gz',
    temp('data/external_data/rnaseq/fastq/elegans_L1_rep3.A.r1.fq.gz'),
    temp('data/external_data/rnaseq/fastq/elegans_L1_rep3.B.r1.fq.gz'),
    'data/external_data/rnaseq/fastq/elegans_L1_rep3.r1.fq.gz',
    'data/external_data/rnaseq/fastq/elegans_L2_rep1.r1.fq.gz',
    'data/external_data/rnaseq/fastq/elegans_L3_rep1.r1.fq.gz',
    'data/external_data/rnaseq/fastq/elegans_L4_rep1.r1.fq.gz',
    'data/external_data/rnaseq/fastq/elegans_lateL4_rep1.r1.fq.gz',
    temp('data/external_data/rnaseq/fastq/elegans_YA_rep1.A.r1.fq.gz'),
    temp('data/external_data/rnaseq/fastq/elegans_YA_rep1.B.r1.fq.gz'),
    'data/external_data/rnaseq/fastq/elegans_YA_rep1.r1.fq.gz',
  shell:
    '''
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/001/SRR1050771/SRR1050771.fastq.gz -O {output[0]}
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/002/SRR1050772/SRR1050772_1.fastq.gz -O {output[1]}
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/002/SRR1050772/SRR1050772_2.fastq.gz -O {output[2]}
    zcat {output[1]} {output[2]} > {output[3]}
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/003/SRR1050773/SRR1050773_1.fastq.gz -O {output[4]}
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/003/SRR1050773/SRR1050773_2.fastq.gz -O {output[5]}
    zcat {output[4]} {output[5]} > {output[6]}
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/004/SRR1050774/SRR1050774.fastq.gz -O {output[7]}
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/005/SRR1050775/SRR1050775.fastq.gz -O {output[8]}
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/006/SRR1050776/SRR1050776.fastq.gz -O {output[9]}
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/007/SRR1050777/SRR1050777.fastq.gz -O {output[10]}
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/008/SRR1050778/SRR1050778_1.fastq.gz -O {output[11]}
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/008/SRR1050778/SRR1050778_2.fastq.gz -O {output[12]}
    zcat {output[11]} {output[12]} > {output[13]}
    '''


rule kallisto_whole_worm_elegans:
  input:
    "kallisto_index/elegans",
    "data/external_data/rnaseq/fastq/{sample}.r1.fq.gz",
    "species/elegans/gene_annotation/elegans.gene_annotation.gene_transcript_id.txt",
  output:
    "data/external_data/rnaseq/elegans/kallisto/{sample}/abundance.h5",
    "data/external_data/rnaseq/elegans/kallisto/{sample}/abundance.tsv",
    "data/external_data/rnaseq/elegans/kallisto/{sample}/run_info.json",
    "data/external_data/rnaseq/elegans/kallisto/{sample}/abundance.genes.txt",
  resources:
    cpus=10
  params:
    "data/external_data/rnaseq/elegans/kallisto/{sample}"
  shell:
    '''
    mkdir -p {params}
    kallisto quant -i {input[0]} --single -o {params} -l 200 -s 20 -t {resources.cpus} {input[1]}
    python scripts/sum_TPM_per_gene.py {input[2]} {output[1]} {output[3]}
    sort -k 1,1 {output[3]} > {output[3]}.sorted
    mv {output[3]}.sorted {output[3]}
    '''


rule elegans_gene_exp_table_generator_stages:
  input:
    kallisto_out_pgc=expand("data/external_data/rnaseq/elegans/kallisto/{sample}/abundance.genes.txt", sample=[item for item in config["elegans_external_rnaseq_samples"] if "elegans" not in item]),
    kallisto_out_stages=expand("data/external_data/rnaseq/elegans/kallisto/{sample}/abundance.genes.txt", sample=[item for item in config["elegans_external_rnaseq_samples"] if "pgc" not in item]),
  output:
    'data/external_data/rnaseq/elegans/gene_expression.pgc.txt',
    'data/external_data/rnaseq/elegans/gene_expression.stages.txt'
  params:
    filenames_pgc=[item for item in config["elegans_external_rnaseq_samples"] if "elegans" not in item],
    filenames_stages=[item for item in config["elegans_external_rnaseq_samples"] if "pgc" not in item]
  shell:
    '''
    echo {params.filenames_pgc} | tr " " "\t" > {output[0]}
    paste {input.kallisto_out_pgc} | awk '{{ for (i=3;i<=NF;i+=2) $i="" }} 1' | tr " " "\t" >> {output[0]}
    echo {params.filenames_stages} | tr " " "\t" > {output[1]}
    paste {input.kallisto_out_stages} | awk '{{ for (i=3;i<=NF;i+=2) $i="" }} 1' | tr " " "\t" >> {output[1]}
    '''


