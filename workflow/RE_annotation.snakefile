rule IGV_tracks_relmapping:
  input:
    'relmapping/atac/tg_se.bwa_se_ce11.rm_unmapped.rm_contigs_ce11.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all/annot_ce_atac_wt_ya_rep1.tg_se.bwa_se_ce11.rm_unmapped.rm_contigs_ce11.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all_treat_pileup.bw',
    'relmapping/atac/tg_se.bwa_se_ce11.rm_unmapped.rm_contigs_ce11.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all/annot_ce_atac_wt_ya_rep2.tg_se.bwa_se_ce11.rm_unmapped.rm_contigs_ce11.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all_treat_pileup.bw',
    'relmapping/atac/tg_se.bwa_se_ce11.rm_unmapped.rm_contigs_ce11.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all/annot_ce_atac_glp1_ya_rep1.tg_se.bwa_se_ce11.rm_unmapped.rm_contigs_ce11.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all_treat_pileup.bw',
    'relmapping/atac/tg_se.bwa_se_ce11.rm_unmapped.rm_contigs_ce11.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all/annot_ce_atac_glp1_ya_rep2.tg_se.bwa_se_ce11.rm_unmapped.rm_contigs_ce11.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all_treat_pileup.bw',
    'relmapping/atac/tg_se.bwa_se_cb4.rm_unmapped.rm_contigs_cb4.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all/annot_cb_atac_wt_ya_rep1.tg_se.bwa_se_cb4.rm_unmapped.rm_contigs_cb4.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all_treat_pileup.bw',
    'relmapping/atac/tg_se.bwa_se_cb4.rm_unmapped.rm_contigs_cb4.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all/annot_cb_atac_wt_ya_rep2.tg_se.bwa_se_cb4.rm_unmapped.rm_contigs_cb4.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all_treat_pileup.bw',
    'relmapping/atac/tg_se.bwa_se_cb4.rm_unmapped.rm_contigs_cb4.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all/annot_cb_atac_glp1_ya_rep1.tg_se.bwa_se_cb4.rm_unmapped.rm_contigs_cb4.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all_treat_pileup.bw',
    'relmapping/atac/tg_se.bwa_se_cb4.rm_unmapped.rm_contigs_cb4.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all/annot_cb_atac_glp1_ya_rep2.tg_se.bwa_se_cb4.rm_unmapped.rm_contigs_cb4.rm_q10.macs2_se_extsize150_shiftm75_keepdup_all_treat_pileup.bw',
  output:
    'relmapping/processed_tracks/elegans_atac_ya_wt.bw',
    'relmapping/processed_tracks/elegans_atac_ya_glp1.bw',
    'relmapping/processed_tracks/briggsae_atac_ya_wt.bw',
    'relmapping/processed_tracks/briggsae_atac_ya_glp1.bw',
  resources:
    cpus=10
  shell:
    '''
    bigwigCompare --bigwig1 {input[0]} --bigwig2 {input[1]} --operation mean -p {resources.cpus} -o {output[0]} -of bigwig -bs 1
    bigwigCompare --bigwig1 {input[2]} --bigwig2 {input[3]} --operation mean -p {resources.cpus} -o {output[1]} -of bigwig -bs 1
    bigwigCompare --bigwig1 {input[4]} --bigwig2 {input[5]} --operation mean -p {resources.cpus} -o {output[2]} -of bigwig -bs 1
    bigwigCompare --bigwig1 {input[6]} --bigwig2 {input[7]} --operation mean -p {resources.cpus} -o {output[3]} -of bigwig -bs 1
    '''
  

rule ce_gl_elements_annotation:
  input:
    "relmapping/annot_ce/reg_elements_ce.bed"
  output:
    parsed_RE = "RE_annotation/reg_elements_all.elegans.bed",
    GL_RE_diffbind = "RE_annotation/reg_elements_all.elegans.gl_specific.diffbind.bed",
    GL_RE = "RE_annotation/reg_elements_all.elegans.gl_specific.bed",
    nonGL_RE = "RE_annotation/reg_elements_all.elegans.not_gl_specific.bed"
  conda:
    "env/diffbind.yaml"
  shell:
    '''
    sed '1d' {input} | awk -F'[\t=;]' 'BEGIN{{OFS="\t";}}{{print $1, $2, $3, "ce.re" NR "." $7, $(NF-4), $(NF-3)}}' > {output.parsed_RE}
    Rscript scripts/diffbind_glre.R
    intersectBed -a {output.parsed_RE} -b {output.GL_RE_diffbind} -u > {output.GL_RE}
    intersectBed -a {output.parsed_RE} -b {output.GL_RE_diffbind} -v > {output.nonGL_RE}
    '''


rule ce_gl_motif_enrichment:
  input:
    genome = "species/elegans/genome/elegans.fa",
    GL_RE = "RE_annotation/reg_elements_all.elegans.gl_specific.bed",
    nonGL_RE = "RE_annotation/reg_elements_all.elegans.not_gl_specific.bed"
  output:
    GL_RE_fasta = temp("RE_annotation/reg_elements_all.elegans.gl_specific.fasta"),
    nonGL_RE_fasta = temp("RE_annotation/reg_elements_all.elegans.not_gl_specific.fasta"),
    meme_results_out = "meme/elegans.GL_specific_vs_nonGL_specific/meme_out/meme.txt",
    meme_results_tsv = "meme/elegans.GL_specific_vs_nonGL_specific/summary.tsv",
  params:
    directory("meme/elegans.GL_specific_vs_nonGL_specific"),
  resources:
    ntasks=6
  shell:
    '''
    fastaFromBed -fi {input.genome} -bed {input.GL_RE} -fo {output.GL_RE_fasta} -s
    fastaFromBed -fi {input.genome} -bed {input.nonGL_RE} -fo {output.nonGL_RE_fasta} -s
    meme-chip -oc {params} -neg {output.nonGL_RE_fasta} -spamo-skip -fimo-skip -meme-p {resources.ntasks} -meme-nmotifs 6 -meme-minw 5 -meme-maxw 20 -dreme-m 0 {output.GL_RE_fasta}
    '''


rule motif_mapping:
  input:
    meme_results = "meme/elegans.GL_specific_vs_nonGL_specific/meme_out/meme.txt",
    meme_summary = "meme/elegans.GL_specific_vs_nonGL_specific/summary.tsv",
    genome = 'species/{sample}/genome/{sample}.fa',
    chr_size = 'species/{sample}/genome/{sample}.chrom.sizes.txt'
  output:
    fimo_m1 = directory("motif_enrichment/{sample}/fimo_{sample}_m1"),
    fimo_m2 = directory("motif_enrichment/{sample}/fimo_{sample}_m2"),
    m1_bed = "motif_enrichment/{sample}/{sample}.m1.bed",
    m2_bed = "motif_enrichment/{sample}/{sample}.m2.bed",
    m1_permissive_bed = "motif_enrichment/{sample}/{sample}.m1.permissive.bed",
    m2_permissive_bed = "motif_enrichment/{sample}/{sample}.m2.permissive.bed",
    m1m2_dist_to_parse = "motif_enrichment/{sample}/{sample}.m1m2_clusters.distance.to_parse",
    m1m2_dist = "motif_enrichment/{sample}/{sample}.m1m2_clusters.distance",
    m1m2_orientation = "motif_enrichment/{sample}/{sample}.m1m2_clusters.motif_orientation.bed",
    m1m2_pairs = "motif_enrichment/{sample}/{sample}.m1m2_clusters.bed"
  shell:
    '''
    m1_consensus=$(sed -n 3p {input.meme_summary} | awk '{{print $3}}')
    m2_consensus=$(sed -n 4p {input.meme_summary} | awk '{{print $3}}')
    fimo --max-stored-scores 1000000 --motif $m1_consensus --oc {output.fimo_m1} --thresh 1e-3 {input.meme_results} {input.genome}
    fimo --max-stored-scores 1000000 --motif $m2_consensus --oc {output.fimo_m2} --thresh 1e-3 {input.meme_results} {input.genome}
    #sed '1d' {output.fimo_m1}/fimo.tsv | sed '/^$/d' | sed '/^#/d' | awk 'BEGIN{{OFS="\t";}}{{if ($7 > 1) print $3, $4, $5, $7, 0, $6}}' | sort -k 1,1 -k2,2n > {output.m1_bed}
    #sed '1d' {output.fimo_m2}/fimo.tsv | sed '/^$/d' | sed '/^#/d' | awk 'BEGIN{{OFS="\t";}}{{if ($7 > 5) print $3, $4, $5, $7, 0, $6}}' | sort -k 1,1 -k2,2n > {output.m2_bed}
    sed '1d' {output.fimo_m1}/fimo.tsv | sed '/^$/d' | sed '/^#/d' | awk 'BEGIN{{OFS="\t";}}{{if ($8 < 0.0005) print $3, $4, $5, $7, 0, $6}}' | sort -k 1,1 -k2,2n > {output.m1_bed}
    sed '1d' {output.fimo_m2}/fimo.tsv | sed '/^$/d' | sed '/^#/d' | awk 'BEGIN{{OFS="\t";}}{{if ($8 < 0.0005) print $3, $4, $5, $7, 0, $6}}' | sort -k 1,1 -k2,2n > {output.m2_bed}
    sed '1d' {output.fimo_m1}/fimo.tsv | sed '/^$/d' | grep -v "#" | awk 'BEGIN{{OFS="\t";}}{{print $3, $4, $5, $7, 0, $6}}' | sort -k 1,1 -k2,2n > {output.m1_permissive_bed}
    sed '1d' {output.fimo_m2}/fimo.tsv | sed '/^$/d' | grep -v "#" | awk 'BEGIN{{OFS="\t";}}{{print $3, $4, $5, $7, 0, $6}}' | sort -k 1,1 -k2,2n > {output.m2_permissive_bed}
    closestBed -a {output.m1_bed} -b {output.m2_bed} -d -io -k 3 -g {input.chr_size} | awk '{{if ($13 < 33 && $13 > 11) print}}' > {output.m1m2_dist_to_parse}
    python scripts/motif_cluster_definition.closest_motif.py {output.m1m2_dist_to_parse} {output.m1m2_dist}
    motif_1_len=$(awk '{{print $3-$2+1}}' {output.m1_bed} | sort | uniq)
    motif_2_len=$(awk '{{print $3-$2+1}}' {output.m2_bed} | sort | uniq)
    min_dist=$(echo $((motif_1_len+motif_2_len)))
    awk 'BEGIN{{OFS="\t";}}{{print $1, $2, $3, $6, $8, $9, $12}}' {output.m1m2_dist} | awk -v dist="$min_dist" 'BEGIN{{OFS="\t";}}{{if ($2 < $5) print $1, $2, $6, "cluster" NR "_" $6-$2+1-dist "bp_ele", "m1" $4 "_m2" $7, "+"; else print $1, $5, $3, "cluster" NR "_" $3-$5+1-dist "bp_ele", "m2" $7 "_m1" $4, "-"}}' > {output.m1m2_orientation}
    awk 'BEGIN{{OFS="\t";}}{{if ($5 == "m1-_m2-" || $5 == "m2+_m1+") print $1, $2, $3, $4, "tandem_m2m1", $6; else if ($5 == "m1+_m2-" || $5 == "m2+_m1-") print $1, $2, $3, $4, "convergent", $6; else if ($5 == "m1+_m2+" || $5 == "m2-_m1-") print $1, $2, $3, $4, "tandem_m1m2", $6; else if ($5 == "m1-_m2+" || $5 == "m2-_m1+") print $1, $2, $3, $4, "divergent", $6}}' {output.m1m2_orientation} > {output.m1m2_pairs}
    '''

rule motif_pair_GL_promoter_enrichment:
  input:
    'RE_annotation/reg_elements_all.elegans.gl_specific.bed',
    'species/elegans/genome/elegans.chrom.sizes.txt',
    'species/elegans/gene_annotation/elegans.gene_annotation.coding_gene_longest_isoform.bed',
    'motif_enrichment/elegans/elegans.m1m2_clusters.bed',
  output:
    'RE_shuffling/coding_prom_GL.bed',
    'RE_shuffling/coding_prom_GL.shuffle.bed',
    'RE_shuffling/coding_prom_GL.shuffling.results',
    'plots/m1m2_permutation_test.pdf',
  params:
    'RE_shuffling/coding_prom_GL.shuffle',
  shell:
    '''
    grep coding_promoter {input[0]} > {output[0]}
    for i in {{1..1000}}
    do
    shuffleBed -i {output[0]} -g {input[1]} -chrom -noOverlapping -excl {input[2]} > {output[1]}
    for arrangement in divergent convergent tandem_m1m2 tandem_m2m1
    do
    grep $arrangement {input[3]} | sort -k 1,1 -k2,2n | intersectBed -a {output[1]} -b stdin -u | wc -l >> {params}.$arrangement
    done
    done
    echo -e "hits\thigher_in_shuffled" > {output[2]}
    for arrangement in divergent convergent tandem_m1m2 tandem_m2m1
    do
    arrangement_hits=$(grep $arrangement {input[3]} | sort -k 1,1 -k2,2n | intersectBed -a {output[0]} -b stdin -u | wc -l)
    higher_in_shuffled=$(awk -v var="$arrangement_hits" '{{if ($1 >= var) print}}' {params}.$arrangement | wc -l)
    echo -e $arrangement"\t"$arrangement_hits"\t"$higher_in_shuffled >> {output[2]}
    done
    Rscript scripts/m1m2_permutation.R
    '''


rule motif_pair_stats:
  input:
    m1m2_pairs = "motif_enrichment/{sample}/{sample}.m1m2_clusters.bed",
    m1_bed = "motif_enrichment/{sample}/{sample}.m1.bed",
    m2_bed = "motif_enrichment/{sample}/{sample}.m2.bed"
  output:
    m1m2_arrangements = "motif_enrichment/{sample}/{sample}.m1m2_clusters.arrangements",
    m1m2_arrangements_n = "motif_enrichment/{sample}/{sample}.m1m2_clusters.arrangements_all.summary",
    m1m2_dist_path = directory("motif_enrichment/{sample}/cluster_distance"),
    m1m2_arrangements_by_dist = "motif_enrichment/{sample}/{sample}.m1m2_clusters.arrangements_by_distance.summary"
  shell:
    '''
    motif_1_len=$(awk '{{print $3-$2+1}}' {input.m1_bed} | sort | uniq)
    motif_2_len=$(awk '{{print $3-$2+1}}' {input.m2_bed} | sort | uniq)
    cut -f 5 {input.m1m2_pairs} | sort | uniq > {output.m1m2_arrangements}
    echo -e "motif_arrangement\tN_clusters" > {output.m1m2_arrangements_n}
    while read line
    do
    arrangement_n=$(grep $line {input.m1m2_pairs} | wc -l)
    echo -e $line"\t"$arrangement_n >> {output.m1m2_arrangements_n}
    done < {output.m1m2_arrangements}
    mkdir {output.m1m2_dist_path}
    header_summary="motif_distance\t"
    while read line
    do
    header_summary+=$line"_cluster\t"
    done < {output.m1m2_arrangements}
    echo -e ${{header_summary::-1}} > {output.m1m2_arrangements_by_dist}
    max_dist=$(cut -f 4 {input.m1m2_pairs} | awk -F"_" '{{print $2}}' | sort | uniq | sed 's/..$//' | sort -n | tail -n 1)
    min_dist=$(cut -f 4 {input.m1m2_pairs} | awk -F"_" '{{print $2}}' | sort | uniq | sed 's/..$//' | sort -n | head -n 1)
    for dist in $(seq $min_dist $max_dist)
    do
    distance=$((dist))bp
    distance_1=$((dist))bp_
    awk -v var_dist="$distance_1" '{{if ($4 ~ var_dist) print}}' {input.m1m2_pairs} > {output.m1m2_dist_path}/m1m2_clusters.$distance.bed
    distance_line=$distance"\t"
    while read line
    do
    arrangement_n=$(grep $line {output.m1m2_dist_path}/m1m2_clusters.$distance.bed | wc -l || true)
    echo $arrangement_n
    distance_line+=$arrangement_n"\t"
    done < {output.m1m2_arrangements}
    echo -e ${{distance_line::-1}} >> {output.m1m2_arrangements_by_dist}
    done
    '''

rule ce_motif_RE_overlap:
  input:
    GL_RE = "RE_annotation/reg_elements_all.elegans.gl_specific.bed",
    nonGL_RE = "RE_annotation/reg_elements_all.elegans.not_gl_specific.bed",
    m1m2_pairs = "motif_enrichment/elegans/elegans.m1m2_clusters.bed"
  output:
    motif_re_summary = "RE_features/elegans.motif_GL_overlap.summary",
    GL_RE_motif = "RE_annotation/reg_elements_all.elegans.gl_specific.m1m2.bed",
    nonGL_RE_motif = "RE_annotation/reg_elements_all.elegans.not_gl_specific.m1m2.bed",
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

rule serizay_RE_annotation_comparison:
  input:
    'data/external_data/serizay_RE_annotation.txt',
    'RE_annotation/reg_elements_all.elegans.gl_specific.bed',
    'motif_enrichment/elegans/elegans.m1m2_clusters.bed',
    'relmapping/processed_tracks/elegans_atac_ya_wt.bw',
    'relmapping/processed_tracks/elegans_atac_ya_glp1.bw',
    'RE_annotation/reg_elements_all.elegans.bed',
    'RE_annotation/reg_elements_all.elegans.not_gl_specific.m1m2.bed'
  output:
    'data/external_data/serizay_RE_annotation.GL_specific.bed',
    'RE_annotation_comparison/current_specific.GLRE.bed',
    'RE_annotation_comparison/serizay_specific.GLRE.bed',
    'RE_annotation_comparison/shared.GLRE.bed',
    'RE_annotation_comparison/current_vs_serizay_RE_annotation.txt',
    temp('plots/current_vs_serizay_RE_annotation.heatmap.mat.gz'),
    'plots/current_vs_serizay_RE_annotation.heatmap.pdf',
    'RE_annotation_comparison/current_specific.GLRE.annot_in_serizay.txt',
    'RE_annotation_comparison/annot_specific.GLRE.absent_in_other_annot.txt',
    'RE_annotation_comparison/non_GL_m1m2_RE.serizay_annotation.txt'
  resources:
    cpus=8
  shell:
    '''
    sed 's/^...//' {input[0]} | awk 'BEGIN{{OFS="\t";}}{{if ($32 == "Germline") print $1, $2, $3, $4}}' > {output[0]}
    intersectBed -a {input[1]} -b {output[0]} -v > {output[1]}
    intersectBed -b {input[1]} -a {output[0]} -v > {output[2]}
    intersectBed -a {input[1]} -b {output[0]} -u > {output[3]}
    echo -e 'GLRE\tGLRE_m1m2\tGL_promoters\tGL_promoters_m1m2' > {output[4]}
    current_RE=$(wc -l < {output[1]})
    serizay_RE=$(wc -l < {output[2]})
    shared_RE=$(wc -l < {output[3]})
    current_m1m2=$(intersectBed -a {output[1]} -b {input[2]} -u | wc -l)
    serizay_m1m2=$(intersectBed -a {output[2]} -b {input[2]} -u | wc -l)
    shared_m1m2=$(intersectBed -a {output[3]} -b {input[2]} -u | wc -l)
    current_P=$(grep "coding_promoter" {output[1]} | wc -l)
    serizay_P=$(grep "coding_promoter" {output[2]} | wc -l)
    shared_P=$(grep "coding_promoter" {output[3]} | wc -l)
    current_P_m1m2=$(grep "coding_promoter" {output[1]} | intersectBed -a stdin -b {input[2]} -u | wc -l)
    serizay_P_m1m2=$(grep "coding_promoter" {output[2]} | intersectBed -a stdin -b {input[2]} -u | wc -l)
    shared_P_m1m2=$(grep "coding_promoter" {output[3]} | intersectBed -a stdin -b {input[2]} -u | wc -l)
    echo -e "current_only\t"$current_RE"\t"$current_m1m2"\t"$current_P"\t"$current_P_m1m2 >> {output[4]}
    echo -e "serizay_only\t"$serizay_RE"\t"$serizay_m1m2"\t"$serizay_P"\t"$serizay_P_m1m2 >> {output[4]}
    echo -e "shared_only\t"$shared_RE"\t"$shared_m1m2"\t"$shared_P"\t"$shared_P_m1m2 >> {output[4]}
    computeMatrix scale-regions -R {output[1]} {output[2]} {output[3]} -S {input[3]} {input[4]} --beforeRegionStartLength 500 --afterRegionStartLength 500 --regionBodyLength 200 --startLabel start --endLabel end --binSize 10 --samplesLabel wt glp-1 -p {resources.cpus} -o {output[5]}
    plotHeatmap -m {output[5]} -out {output[6]} --averageTypeSummaryPlot mean --missingDataColor 0.5 --colorList 'white,purple' 'white,black' --plotTitle "elegans RE by annotation" -min 0 -max 8 8 --heatmapHeight 15 --heatmapWidth 6 --xAxisLabel "" --whatToShow 'heatmap and colorbar' --startLabel start --endLabel end
    sed 's/^...//' {input[0]} | intersectBed -a stdin -b {output[1]} -u | cut -f 23,24,25,26,27,32 > {output[7]}
    current_absent_in_serizay=$(sed 's/^...//' {input[0]} | cut -f 1,2,3 | intersectBed -a {output[1]} -b stdin -v | wc -l)
    serizay_absent_in_current=$(intersectBed -a {output[2]} -b {input[5]} -v | wc -l)
    echo -e "specific\tabsent_in_other" > {output[8]}
    echo -e "current\t"$current_RE"\t"$current_absent_in_serizay"\nserizay\t"$serizay_RE"\t"$serizay_absent_in_current >> {output[8]}
    sed 's/^...//' {input[0]} | intersectBed -a stdin -b {input[6]} -u | cut -f 32 | sort | uniq -c > {output[9]}
    Rscript scripts/RE_annotation_comparison.R
    '''


rule ce_motif_repeats_overlap:
  input:
    'species/elegans/repeats_dfam/elegans_dfam.repeats.id.bed',
    'motif_enrichment/elegans/elegans.m1m2_clusters.bed',
  output:
    temp('motif_enrichment/elegans/repeats_all.txt'),
    'motif_enrichment/elegans/repeat_m1m2_overlap.summary',
    temp('motif_enrichment/elegans/current_repeat'),
  shell:
    '''
    cut -f 4 {input[0]} | sort | uniq > {output[0]}
    header="repeat\ttotal_repeats\t"
    for arrang in convergent divergent tandem_m1m2 tandem_m2m1
    do
    header+=$arrang"\t"
    done
    echo -e ${{header::-1}} > {output[1]}
    while read line
    do
    newline=$line"\t"
    awk -v var="$line" '{{if ($4 == var) print}}' {input[0]} > {output[2]}
    newline+=$(wc -l < {output[2]})"\t"
    for arrang in convergent divergent tandem_m1m2 tandem_m2m1
    do
    newline+=$(grep $arrang {input[1]} | intersectBed -a {output[2]} -b stdin -u | wc -l)"\t"
    done
    echo -e ${{newline::-1}} >> {output[1]}
    done < {output[0]}
    '''


rule ce_GL_RE_repeats_overlap:
  input:
    'species/elegans/repeats_dfam/elegans_dfam.repeats.id.bed',
    'RE_annotation/reg_elements_all.elegans.gl_specific.bed',
    'RE_annotation/reg_elements_all.elegans.not_gl_specific.bed',
  output:
    temp('motif_enrichment/elegans/repeats_all.txt'),
    'motif_enrichment/elegans/repeat_GL_RE_overlap.summary',
    temp('motif_enrichment/elegans/current_repeat'),
    'motif_enrichment/elegans/repeat_GL_promoter_overlap.summary',
  shell:
    '''
    cut -f 4 {input[0]} | sort | uniq > {output[0]}
    echo -e "repeat_class\trepeat_number\trep_in_GL_specific_peaks\tGL_specific_peaks\trep_in_non_GL_specific_peaks\tnon_GL_specific_peaks" > {output[1]}
    echo -e "repeat_class\trepeat_number\trep_in_GL_specific_peaks\tGL_specific_peaks\trep_in_non_GL_specific_peaks\tnon_GL_specific_peaks" > {output[3]}
    while read line
    do
    newline=$line"\t"
    awk -v var="$line" '{{if ($4 == var) print}}' {input[0]} > {output[2]}
    newline+=$(wc -l < {output[2]})"\t"
    newline+=$(intersectBed -a {input[1]} -b {output[2]} -u | wc -l)"\t"
    newline+=$(wc -l < {input[1]})"\t"
    newline+=$(intersectBed -a {input[2]} -b {output[2]} -u | wc -l)"\t"
    newline+=$(wc -l < {input[2]})
    echo -e $newline >> {output[1]}
    newline=$line"\t"
    newline+=$(wc -l < {output[2]})"\t"
    newline+=$(grep coding_promoter {input[1]} | intersectBed -a stdin -b {output[2]} -u | wc -l)"\t"
    newline+=$(wc -l < {input[1]})"\t"
    newline+=$(grep coding_promoter {input[2]} | intersectBed -a stdin -b {output[2]} -u | wc -l)"\t"
    newline+=$(wc -l < {input[2]})
    echo -e $newline >> {output[3]}
    done < {output[0]}
    '''


#rule ce_RE_repeats_overlap:
#  input:
#    'species/elegans/repeats_dfam/elegans_dfam.repeats.id.bed',
#    'RE_annotation/reg_elements_all.elegans.bed',
#    'RE_annotation/reg_elements_all.elegans.gl_specific.bed',
#  output:
#    temp('motif_enrichment/elegans/repeats_all.txt'),
#    'motif_enrichment/elegans/repeat_RE_overlap.summary',
#    temp('motif_enrichment/elegans/current_repeat'),
#  shell:
#    '''
#    cut -f 4 {input[0]} | sort | uniq > {output[0]}
#    echo -e "repeat_class\ttotal_REs\ttotal_repeats\tGL_specific_REs\tGL_specific_RE_repeat_overlap" > {output[1]}
#    while read line
#    do
#    newline=$line"\t"
#    awk -v var="$line" '{{if ($4 == var) print}}' {input[0]} > {output[2]}
#    newline+=$(wc -l < {input[1]})"\t"
#    newline+=$(wc -l < {output[2]})"\t"
#    newline+=$(wc -l < {input[2]})"\t"
#    newline+=$(intersectBed -a {input[2]} -b {output[2]} -u | wc -l)
#    echo -e $newline >> {output[1]}
#    done < {output[0]}
#    '''


rule ce_gl_m1m2_motif_enrichment:
  input:
    genome = "species/elegans/genome/elegans.fa",
    GL_RE = "RE_annotation/reg_elements_all.elegans.gl_specific.bed",
    nonGL_RE = "RE_annotation/reg_elements_all.elegans.not_gl_specific.bed",
    m1m2_pairs = "motif_enrichment/elegans/elegans.m1m2_clusters.bed"
  output:
    GL_RE_m1m2_fasta = temp("RE_annotation/reg_elements_all.elegans.gl_specific.m1m2.fasta"),
    nonGL_RE_m1m2_fasta = temp("RE_annotation/reg_elements_all.elegans.not_gl_specific.m1m2.fasta"),
    GL_RE_no_m1m2_fasta = temp("RE_annotation/reg_elements_all.elegans.gl_specific.no_m1m2.fasta"),
    meme_results_1 = directory("meme/elegans.GL_specific_m1m2_vs_nonGL_specific_m1m2"),
    meme_results_2 = directory("meme/elegans.GL_specific_no_m1m2_vs_GL_specific_m1m2")
  resources:
    ntasks=6
  shell:
    '''
    intersectBed -a {input.GL_RE} -b {input.m1m2_pairs} -u | fastaFromBed -fi {input.genome} -bed stdin -fo {output.GL_RE_m1m2_fasta} -s
    intersectBed -a {input.nonGL_RE} -b {input.m1m2_pairs} -u | fastaFromBed -fi {input.genome} -bed stdin -fo {output.nonGL_RE_m1m2_fasta} -s
    intersectBed -a {input.GL_RE} -b {input.m1m2_pairs} -v | fastaFromBed -fi {input.genome} -bed stdin -fo {output.GL_RE_no_m1m2_fasta} -s
    meme -oc {output.meme_results_1} -objfun de -dna -revcomp -nmotifs 10 -p {resources.ntasks} -neg {output.nonGL_RE_m1m2_fasta} {output.GL_RE_m1m2_fasta}
    meme -oc {output.meme_results_2} -objfun de -dna -revcomp -nmotifs 10 -p {resources.ntasks} -neg {output.GL_RE_m1m2_fasta} {output.GL_RE_no_m1m2_fasta}
    '''


rule ce_RE_GC_content:
  input:
    GL_RE = "RE_annotation/reg_elements_all.elegans.gl_specific.bed",
    nonGL_RE = "RE_annotation/reg_elements_all.elegans.not_gl_specific.bed",
    m1m2_pairs = "motif_enrichment/elegans/elegans.m1m2_clusters.bed",
    genome = "species/elegans/genome/elegans.fa"
  output:
    GL_promoter_motif_nuc = "RE_features/reg_elements_all.elegans.gl_specific.promoters.m1m2.nuc",
    nonGL_promoter_motif_nuc = "RE_features/reg_elements_all.elegans.not_gl_specific.promoters.m1m2.nuc",
    GL_promoter_motif_fa = temp("RE_features/reg_elements_all.elegans.gl_specific.promoters.m1m2.fa"),
    nonGL_promoter_motif_fa = temp("RE_features/reg_elements_all.elegans.not_gl_specific.promoters.m1m2.fa"),
    GL_promoter_motif_cpg = "RE_features/reg_elements_all.elegans.gl_specific.promoters.m1m2.cpg",
    nonGL_promoter_motif_cpg = "RE_features/reg_elements_all.elegans.not_gl_specific.promoters.m1m2.cpg",
    GL_promoter_no_motif_nuc = "RE_features/reg_elements_all.elegans.gl_specific.promoters.no_m1m2.nuc",
    nonGL_promoter_no_motif_nuc = "RE_features/reg_elements_all.elegans.not_gl_specific.promoters.no_m1m2.nuc",
    GL_promoter_no_motif_fa = temp("RE_features/reg_elements_all.elegans.gl_specific.promoters.no_m1m2.fa"),
    nonGL_promoter_no_motif_fa = temp("RE_features/reg_elements_all.elegans.not_gl_specific.promoters.no_m1m2.fa"),
    GL_promoter_no_motif_cpg = "RE_features/reg_elements_all.elegans.gl_specific.promoters.no_m1m2.cpg",
    nonGL_promoter_no_motif_cpg = "RE_features/reg_elements_all.elegans.not_gl_specific.promoters.no_m1m2.cpg"
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

rule motif_3_elegans:
  input:
    meme_results = "meme/elegans.GL_specific_vs_nonGL_specific/meme_out/meme.txt",
    meme_summary = "meme/elegans.GL_specific_vs_nonGL_specific/summary.tsv",
    genome = "species/elegans/genome/elegans.fa",
    GL_RE = "RE_annotation/reg_elements_all.elegans.gl_specific.bed",
    nonGL_RE = "RE_annotation/reg_elements_all.elegans.not_gl_specific.bed",
    m1m2_pairs = "motif_enrichment/elegans/elegans.m1m2_clusters.bed"
  output:
    fimo_m3 = directory("motif_enrichment/elegans/fimo_elegans_m3"),
    m3_bed = "motif_enrichment/elegans/elegans.m3.bed",
    m3_stats = "motif_enrichment/elegans/elegans.motif3_association.txt"
  shell:
    '''
    m3_consensus=$(sed -n 2p {input.meme_summary} | awk '{{print $3}}')
    fimo --max-stored-scores 1000000 --motif $m3_consensus --oc {output.fimo_m3} --thresh 1e-4 {input.meme_results} {input.genome}
    sed '1d' {output.fimo_m3}/fimo.tsv | sed '/^$/d' | sed '/^#/d' | awk 'BEGIN{{OFS="\t";}}{{print $3, $4, $5, $7, 0, $6}}' | sort -k 1,1 -k2,2n > {output.m3_bed}
    echo -e "GL_with_m3\tGL_without_m3\tnonGL_with_m3\tnonGL_without_m3" > {output.m3_stats}
    GL_with_m3=$(intersectBed -a {input.GL_RE} -b {output.m3_bed} -u | intersectBed -a {input.m1m2_pairs} -b stdin -u | cut -f 5 | sort | uniq -c | grep -E 'divergent|tandem_m2m1')
    GL_without_m3=$(intersectBed -a {input.GL_RE} -b {output.m3_bed} -v | intersectBed -a {input.m1m2_pairs} -b stdin -u | cut -f 5 | sort | uniq -c | grep -E 'divergent|tandem_m2m1')
    nonGL_with_m3=$(intersectBed -a {input.nonGL_RE} -b {output.m3_bed} -u | intersectBed -a {input.m1m2_pairs} -b stdin -u | cut -f 5 | sort | uniq -c | grep -E 'divergent|tandem_m2m1')
    nonGL_without_m3=$(intersectBed -a {input.nonGL_RE} -b {output.m3_bed} -v | intersectBed -a {input.m1m2_pairs} -b stdin -u | cut -f 5 | sort | uniq -c | grep -E 'divergent|tandem_m2m1')
    echo $GL_with_m3 $GL_without_m3 $nonGL_with_m3 $nonGL_without_m3 | awk 'BEGIN{{OFS="\t";}}{{print $2, $1, $5, $9, $13}}' >> {output.m3_stats}
    echo $GL_with_m3 $GL_without_m3 $nonGL_with_m3 $nonGL_without_m3 | awk 'BEGIN{{OFS="\t";}}{{print $4, $3, $7, $11, $15}}' >> {output.m3_stats}
    '''


rule motif_3_other_species:
  input:
    meme_results = "meme/elegans.GL_specific_vs_nonGL_specific/meme_out/meme.txt",
    meme_summary = "meme/elegans.GL_specific_vs_nonGL_specific/summary.tsv",
    genome = "species/{sample}/genome/{sample}.fa",
  output:
    fimo_m3 = directory("motif_enrichment/{sample}/fimo_{sample}_m3"),
    m3_bed = "motif_enrichment/{sample}/{sample}.m3.bed",
  shell:
    '''
    m3_consensus=$(sed -n 2p {input.meme_summary} | awk '{{print $3}}')
    fimo --max-stored-scores 1000000 --motif $m3_consensus --oc {output.fimo_m3} --thresh 1e-4 {input.meme_results} {input.genome}
    sed '1d' {output.fimo_m3}/fimo.tsv | sed '/^$/d' | sed '/^#/d' | awk 'BEGIN{{OFS="\t";}}{{print $3, $4, $5, $7, 0, $6}}' | sort -k 1,1 -k2,2n > {output.m3_bed}
    '''


rule ce_motif_associated_genes:
  input:
    annot_ce = "relmapping/annot_ce/reg_elements_ce.bed",
    re_all = "RE_annotation/reg_elements_all.elegans.bed",
    GL_RE = "RE_annotation/reg_elements_all.elegans.gl_specific.bed",
    nonGL_RE = "RE_annotation/reg_elements_all.elegans.not_gl_specific.bed",
    m1m2_pairs = "motif_enrichment/elegans/elegans.m1m2_clusters.bed"
  output:
    promoters_genes = "promoter_annotation/promoters_all.elegans.bed",
    m1m2_GL_genes = "promoter_annotation/promoters_gl_specific.elegans.m1m2.any_promoter.genes",
    no_m1m2_GL_genes = "promoter_annotation/promoters_gl_specific.elegans.no_m1m2.any_promoter.genes",
    m1m2_noGL_genes = "promoter_annotation/promoters_not_gl_specific.elegans.m1m2.any_promoter.genes",
    no_m1m2_noGL_genes = "promoter_annotation/promoters_not_gl_specific.elegans.no_m1m2.any_promoter.genes",
    unique_promoters_genes = temp("promoter_annotation/promoters_all.elegans.unique.genes"),
    unique_promoters = "promoter_annotation/promoters_all.elegans.unique.bed",
    m1m2_GL_genes_uniq = "promoter_annotation/promoters_gl_specific.elegans.m1m2.unique_promoter.genes",
    no_m1m2_GL_genes_uniq = "promoter_annotation/promoters_gl_specific.elegans.no_m1m2.unique_promoter.genes",
    m1m2_noGL_genes_uniq = "promoter_annotation/promoters_not_gl_specific.elegans.m1m2.unique_promoter.genes",
    no_m1m2_noGL_genes_uniq = "promoter_annotation/promoters_not_gl_specific.elegans.no_m1m2.unique_promoter.genes"
  shell:
    '''
    grep "annot=coding_promoter" {input.annot_ce} | awk 'BEGIN{{OFS="\t";FS="\t|=|;";}}{{if (length($5) == 35) print $1, $2, $3, substr($5, 0, 14) "," substr($5, 22, 35); else print $1, $2, $3, $5}}' | intersectBed -a stdin -b {input.re_all} -wa -wb | cut -f 1,2,3,4,9,10 > {output.promoters_genes}
    intersectBed -a {output.promoters_genes} -b {input.GL_RE} -u | intersectBed -a stdin -b {input.m1m2_pairs} -u | cut -f 4 | tr "," "\n" | sort > {output.m1m2_GL_genes}
    intersectBed -a {output.promoters_genes} -b	{input.GL_RE} -u | intersectBed	-a stdin -b {input.m1m2_pairs} -v | cut	-f 4 | tr "," "\n" | sort | join -v 1 - {output.m1m2_GL_genes} > {output.no_m1m2_GL_genes}
    intersectBed -a {output.promoters_genes} -b	{input.GL_RE} -v | intersectBed	-a stdin -b {input.m1m2_pairs} -u | cut	-f 4 | tr "," "\n" | sort > {output.m1m2_noGL_genes}
    intersectBed -a {output.promoters_genes} -b	{input.GL_RE} -v | intersectBed	-a stdin -b {input.m1m2_pairs} -v | cut	-f 4 | tr "," "\n" | sort | join -v 1 - {output.m1m2_noGL_genes} > {output.no_m1m2_noGL_genes}
    awk 'BEGIN{{OFS="\t";}}{{if ($4 ~ ",") print $1, $2, $3, substr($4, 0, 14), $5, $6 "\\n" $1, $2, $3, substr($4, 16, 29), $5, $6; else print $0}}' {output.promoters_genes} | cut -f 4 | sort | uniq -c | awk '{{if ($1 == 1) print $2}}' > {output.unique_promoters_genes}
    awk	'BEGIN{{OFS="\t";}}{{if ($4 ~ ",") print $1, $2, $3, substr($4, 0, 14), $5, $6 "\\n" $1, $2, $3, substr($4, 16, 29), $5, $6; else print $0}}' {output.promoters_genes} | sort -k 4,4 | join -1 4 -2 1 - {output.unique_promoters_genes} | awk 'BEGIN{{OFS="\t";}}{{print $2, $3, $4, $1, $5, $6}}' > {output.unique_promoters}
    intersectBed -a {output.unique_promoters} -b {input.GL_RE} -u | intersectBed -a stdin -b {input.m1m2_pairs} -u | cut -f 4 | tr "," "\n" | sort > {output.m1m2_GL_genes_uniq}
    intersectBed -a {output.unique_promoters} -b {input.GL_RE} -u | intersectBed -a stdin -b {input.m1m2_pairs} -v | cut -f 4 | tr "," "\n" | sort > {output.no_m1m2_GL_genes_uniq}
    intersectBed -a {output.unique_promoters} -b {input.GL_RE} -v | intersectBed -a stdin -b {input.m1m2_pairs} -u | cut -f 4 | tr "," "\n" | sort > {output.m1m2_noGL_genes_uniq}
    intersectBed -a {output.unique_promoters} -b {input.GL_RE} -v | intersectBed -a stdin -b {input.m1m2_pairs} -v | cut -f 4 | tr "," "\n" | sort > {output.no_m1m2_noGL_genes_uniq}
    '''

rule boeck_expression_parser:
  output:
    gene_ids = "species/elegans/gene_annotation/elegans.gene_annotation.id",
    boeck_data = "data/external_data/boeck_exp_data.txt",
    boeck_data_emb1 = "data/external_data/boeck_exp_data.embryo_series_1.txt",
    boeck_data_emb2 = "data/external_data/boeck_exp_data.embryo_series_2.txt",
    boeck_data_emb3 = "data/external_data/boeck_exp_data.embryo_series_3.txt",
    boeck_data_emb4 = "data/external_data/boeck_exp_data.embryo_series_4.txt",
    boeck_data_postemb = "data/external_data/boeck_exp_data.postembryonic.txt",
  params:
    "data/external_data/boeck_exp_data"
  shell:
    '''
    wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS275/species/c_elegans/PRJNA13758/annotation/c_elegans.PRJNA13758.WS275.geneIDs.txt.gz -O {output.gene_ids}.gz; gunzip {output.gene_ids}.gz
    wget https://genome.cshlp.org/content/suppl/2016/09/20/gr.202663.115.DC1/Supplemental_Table_S2.gz -O {output.boeck_data}.gz; gunzip {output.boeck_data}.gz
    python scripts/boeck_data_parser.py {params} {output.gene_ids}
    '''

rule boeck_unified_exp_parser:
  input:
    "species/elegans/gene_annotation/elegans.gene_annotation.id"
  output:
    temp("data/external_data/boeck_unified_exp_data.to_parse.txt"),
    "data/external_data/boeck_unified_exp_data.txt",
  shell:
    '''
    wget https://genome.cshlp.org/content/suppl/2016/09/20/gr.202663.115.DC1/Supplemental_Table_S13.gz -O {output[0]}.gz; gunzip {output[0]}.gz
    python scripts/boeck_unified_data_parser.py {input} {output[0]} {output[1]}
    '''


rule boeck_gl_expression:
  input:
    m1m2_GL_genes_uniq = "promoter_annotation/promoters_gl_specific.elegans.m1m2.unique_promoter.genes",
    no_m1m2_GL_genes_uniq = "promoter_annotation/promoters_gl_specific.elegans.no_m1m2.unique_promoter.genes",
    m1m2_noGL_genes_uniq = "promoter_annotation/promoters_not_gl_specific.elegans.m1m2.unique_promoter.genes",
    no_m1m2_noGL_genes_uniq = "promoter_annotation/promoters_not_gl_specific.elegans.no_m1m2.unique_promoter.genes",
    boeck_exp = "data/external_data/boeck_exp_data.{sample}.txt"
  output:
    m1m2_GL_genes_uniq_exp = "GL_genes_exp/promoters_gl_specific.elegans.m1m2.unique_promoter.{sample}.exp",
    no_m1m2_GL_genes_uniq_exp = "GL_genes_exp/promoters_gl_specific.elegans.no_m1m2.unique_promoter.{sample}.exp",
    m1m2_noGL_genes_uniq_exp = "GL_genes_exp/promoters_not_gl_specific.elegans.m1m2.unique_promoter.{sample}.exp",
    no_m1m2_noGL_genes_uniq_exp = "GL_genes_exp/promoters_not_gl_specific.elegans.no_m1m2.unique_promoter.{sample}.exp"
  shell:
    '''
    head -n 1 {input.boeck_exp} > {output.m1m2_GL_genes_uniq_exp}
    head -n 1 {input.boeck_exp} > {output.no_m1m2_GL_genes_uniq_exp}
    head -n 1 {input.boeck_exp} > {output.m1m2_noGL_genes_uniq_exp}
    head -n 1 {input.boeck_exp} > {output.no_m1m2_noGL_genes_uniq_exp}
    sed '1d' {input.boeck_exp} | sort -k 1,1 | join -1 1 -2 1 - {input.m1m2_GL_genes_uniq} | tr " " "\t" >> {output.m1m2_GL_genes_uniq_exp}
    sed	'1d' {input.boeck_exp} | sort -k 1,1 | join -1 1 -2 1 -	{input.no_m1m2_GL_genes_uniq} | tr " " "\t" >> {output.no_m1m2_GL_genes_uniq_exp}
    sed	'1d' {input.boeck_exp} | sort -k 1,1 | join -1 1 -2 1 -	{input.m1m2_noGL_genes_uniq} | tr " " "\t" >> {output.m1m2_noGL_genes_uniq_exp}
    sed	'1d' {input.boeck_exp} | sort -k 1,1 | join -1 1 -2 1 -	{input.no_m1m2_noGL_genes_uniq} | tr " " "\t" >> {output.no_m1m2_noGL_genes_uniq_exp}
    '''


rule packer_expression_parser:
  output:
    packer_data_all = "data/external_data/packer_exp_data.all.zip",
    packer_data = "data/external_data/packer_exp_data.txt",
    packer_data_parsed = "data/external_data/packer_exp_data.germline.txt"
  params:
    "data/external_data"
  shell:
    '''
    wget https://science.sciencemag.org/highwire/filestream/731368/field_highwire_adjunct_files/3/aax1971_Tables_S7_S8_S10_S11_S14.zip -O {output.packer_data_all}; unzip -d {params} {output.packer_data_all}
    gunzip {params}/aax1971_Table_S7.gz; mv {params}/aax1971_Table_S7 {output.packer_data}
    python scripts/packer_data_parser.py {output.packer_data} {output.packer_data_parsed}
    '''

rule packer_gl_expression:
  input:
    m1m2_GL_genes_uniq = "promoter_annotation/promoters_gl_specific.elegans.m1m2.unique_promoter.genes",
    no_m1m2_GL_genes_uniq = "promoter_annotation/promoters_gl_specific.elegans.no_m1m2.unique_promoter.genes",
    m1m2_noGL_genes_uniq = "promoter_annotation/promoters_not_gl_specific.elegans.m1m2.unique_promoter.genes",
    no_m1m2_noGL_genes_uniq = "promoter_annotation/promoters_not_gl_specific.elegans.no_m1m2.unique_promoter.genes",
    packer_exp = "data/external_data/packer_exp_data.germline.txt"
  output:
    m1m2_GL_genes_uniq_exp_gl = "GL_genes_exp/promoters_gl_specific.elegans.m1m2.unique_promoter.embryo_germline.exp",
    no_m1m2_GL_genes_uniq_exp_gl = "GL_genes_exp/promoters_gl_specific.elegans.no_m1m2.unique_promoter.embryo_germline.exp",
    m1m2_noGL_genes_uniq_exp_gl = "GL_genes_exp/promoters_not_gl_specific.elegans.m1m2.unique_promoter.embryo_germline.exp",
    no_m1m2_noGL_genes_uniq_exp_gl = "GL_genes_exp/promoters_not_gl_specific.elegans.no_m1m2.unique_promoter.embryo_germline.exp"
  shell:
    '''
    head -n 1 {input.packer_exp} > {output.m1m2_GL_genes_uniq_exp_gl}
    head -n 1 {input.packer_exp} > {output.no_m1m2_GL_genes_uniq_exp_gl}
    head -n 1 {input.packer_exp} > {output.m1m2_noGL_genes_uniq_exp_gl}
    head -n 1 {input.packer_exp} > {output.no_m1m2_noGL_genes_uniq_exp_gl}
    sed '1d' {input.packer_exp} | sort -k 1,1 | join -1 1 -2 1 - {input.m1m2_GL_genes_uniq} | tr " " "\t" >> {output.m1m2_GL_genes_uniq_exp_gl}
    sed '1d' {input.packer_exp} | sort -k 1,1 | join -1 1 -2 1 - {input.no_m1m2_GL_genes_uniq} | tr " " "\t" >> {output.no_m1m2_GL_genes_uniq_exp_gl}
    sed '1d' {input.packer_exp} | sort -k 1,1 | join -1 1 -2 1 - {input.m1m2_noGL_genes_uniq} | tr " " "\t" >> {output.m1m2_noGL_genes_uniq_exp_gl}
    sed '1d' {input.packer_exp} | sort -k 1,1 | join -1 1 -2 1 - {input.no_m1m2_noGL_genes_uniq} | tr " " "\t" >> {output.no_m1m2_noGL_genes_uniq_exp_gl}
    '''

rule motif_bw:
  input:
    'motif_enrichment/{sample}/{sample}.m1m2_clusters.bed',
    'species/{sample}/genome/{sample}.chrom.sizes.txt'
  output:
    temp('motif_enrichment/{sample}/{sample}.m1m2_clusters.bg'),
    'motif_enrichment/{sample}/{sample}.m1m2_clusters.bw'
  shell:
    '''
    genomeCoverageBed -i {input[0]} -g {input[1]} -bga > {output[0]}.unsorted
    LC_COLLATE=C sort -k 1,1 -k2,2n {output[0]}.unsorted > {output[0]}
    rm {output[0]}.unsorted
    bedGraphToBigWig {output[0]} {input[1]} {output[1]}
    '''

#rule GL_gene_annotation_serizay:
#  input:
#    'serizay_data/LCAP_final-annotations_uniqueTPMthreshold.gff3',
#    'promoter_annotation/promoters_gl_specific.elegans.m1m2.any_promoter.genes',
#    'promoter_annotation/promoters_gl_specific.elegans.m1m2.unique_promoter.genes',
#    'promoter_annotation/promoters_gl_specific.elegans.no_m1m2.any_promoter.genes',
#    'promoter_annotation/promoters_not_gl_specific.elegans.m1m2.any_promoter.genes',
#    'promoter_annotation/promoters_not_gl_specific.elegans.m1m2.unique_promoter.genes',
#    'promoter_annotation/promoters_not_gl_specific.elegans.no_m1m2.any_promoter.genes',
#  output:
#    'gene_annotation/elegans.germline_specific.serizay.bed',
#    'gene_annotation/elegans.germline_expressed.serizay.bed',
#    'gene_annotation/elegans.germline_specific.m1m2.any.serizay.bed',
#    'gene_annotation/elegans.germline_specific.m1m2.unique.serizay.bed',
#    'gene_annotation/elegans.germline_specific.no_m1m2.serizay.bed',
#    'gene_annotation/elegans.germline_expressed.m1m2.any.serizay.bed',
#    'gene_annotation/elegans.germline_expressed.m1m2.unique.serizay.bed',
#    'gene_annotation/elegans.germline_expressed.no_m1m2.serizay.bed',
#  shell:
#    '''
#    sed '1,3d' {input[0]} | awk 'BEGIN{{FS="\t|=|;"; OFS="\t";}}{{if ($(NF-2) == "Germline") print $1, $4, $5, $10, $6, $7}}' | sort -k 4,4 > {output[0]}
#    sed	'1,3d' {input[0]} | awk	'BEGIN{{FS="\t|=|;"; OFS="\t";}}{{if ($(NF-2) ~ "Germline" || $(NF-2) ~ "Ubiq") print $1, $4, $5, $10, $6, $7}}' | sort -k 4,4 > {output[1]}
#    '''

rule DESeq_GL_gene_annotation_ce:
  input:
    expand("relmapping/lcap/trim20.bwa_pe_ce11.rm_unmapped_pe.rm_contigs_ce11.rm_q10/{sample}.trim20.bwa_pe_ce11.rm_unmapped_pe.rm_contigs_ce11.rm_q10.bam", sample=config["lcap_samples_ce"]),
    "species/elegans/gene_annotation/elegans.annotations.coding_genes.gtf"
  output:
    "gene_annotation/elegans.genes_DESeq_glp1_vs_wt.txt",
    "gene_annotation/elegans.genes_downregulated_glp1_vs_wt.txt",
    "gene_annotation/elegans.genes_upregulated_glp1_vs_wt.txt",
  conda:
    "env/diffbind.yaml"
  shell:
    '''
    Rscript scripts/DESeq_lcap_ce.R
    '''

rule GL_gene_annotation:
  input:
    'promoter_annotation/promoters_gl_specific.elegans.m1m2.unique_promoter.genes',
    'promoter_annotation/promoters_gl_specific.elegans.no_m1m2.unique_promoter.genes',
    'gene_annotation/elegans.genes_downregulated_glp1_vs_wt.txt',
    'species/elegans/gene_annotation/elegans.gene_annotation.coding_gene_longest_isoform.bed'
  output:
    'gene_annotation/elegans.GL_specific_genes.m1m2_unique.genes',
    'gene_annotation/elegans.GL_specific_genes.no_m1m2_unique.genes',
    'gene_annotation/elegans.GL_specific_genes.m1m2_unique.bed',
    'gene_annotation/elegans.GL_specific_genes.no_m1m2_unique.bed',
    'gene_annotation/elegans.GL_specific_genes.bed',
  shell:
    '''
    join {input[0]} {input[2]} | cut -f 1 | sort > {output[0]}
    join {input[1]} {input[2]} | cut -f 1 | sort > {output[1]}
    sort -k 4,4 {input[3]} | join -1 4 -2 1 - {output[0]} | awk 'BEGIN{{OFS="\t";}}{{print $2, $3, $4, $1, $5, $6, $7, $8, $9, $10, $11, $12}}' > {output[2]}
    sort -k 4,4	{input[3]} | join -1 4 -2 1 - {output[1]} | awk	'BEGIN{{OFS="\t";}}{{print $2, $3, $4, $1, $5, $6, $7, $8, $9, $10, $11, $12}}' > {output[3]}
    sort -k 4,4 {input[3]} | join -1 4 -2 1 - {input[2]} | awk 'BEGIN{{OFS="\t";}}{{print $2, $3, $4, $1, $5, $6, $7, $8, $9, $10, $11, $12}}' > {output[4]}
    '''

rule DESeq_GL_gene_annotation_cb:
  input:
    expand("relmapping/lcap/trim20.bwa_pe_cb4.rm_unmapped_pe.rm_contigs_cb4.rm_q10/{sample}.trim20.bwa_pe_cb4.rm_unmapped_pe.rm_contigs_cb4.rm_q10.bam", sample=config["lcap_samples_cb"]),
    "species/briggsae/gene_annotation/briggsae.annotations.coding_genes.gtf"
  output:
    "gene_annotation/briggsae.genes_DESeq_glp1_vs_wt.txt",
    "gene_annotation/briggsae.genes_downregulated_glp1_vs_wt.txt",
    "gene_annotation/briggsae.genes_upregulated_glp1_vs_wt.txt",
  conda:
    "env/diffbind.yaml"
  shell:
    '''
    Rscript scripts/DESeq_lcap_cb.R
    '''

rule external_atac_download:
  output:
    'data/external_data/atac/fastq/elegans.pgc.rep1.fq.gz',
    'data/external_data/atac/fastq/elegans.pgc.rep2.fq.gz',
    'data/external_data/atac/fastq/elegans.adult_gl.rep1.fq.gz',
    'data/external_data/atac/fastq/elegans.adult_gl.rep2.fq.gz',
  shell:
    '''
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR577/003/SRR5772133/SRR5772133.fastq.gz -O {output[0]}
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR577/004/SRR5772134/SRR5772134.fastq.gz -O {output[1]}
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/033/SRR10564633/SRR10564633.fastq.gz -O {output[2]}
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/057/SRR10564657/SRR10564657.fastq.gz -O {output[3]}
    '''

rule bwa_idx:
  input:
    "species/elegans/genome/elegans.fa"
  output:
    "bwa_idx/elegans.amb",
    "bwa_idx/elegans.ann",
    "bwa_idx/elegans.bwt",
    "bwa_idx/elegans.pac",
    "bwa_idx/elegans.sa"
  params:
    "bwa_idx/elegans"
  shell:
    '''
    bwa index {input} -p {params}
    '''

rule external_atac_aln:
  input:
    'bwa_idx/elegans.amb',
    'bwa_idx/elegans.ann',
    'bwa_idx/elegans.bwt',
    'bwa_idx/elegans.pac',
    'bwa_idx/elegans.sa',
    'species/elegans/genome/elegans.fa',
    'data/external_data/atac/fastq/elegans.{sample}.fq.gz',
  output:
    'data/external_data/atac/aligned/elegans.{sample}.bam'
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

rule external_atac_macs2:
  input:
    'data/external_data/atac/aligned/elegans.{sample}.bam',
    'species/elegans/genome/elegans.chrom.sizes.txt'
  output:
    'data/external_data/atac/peaks/elegans.{sample}/elegans.{sample}_peaks.narrowPeak',
    'data/external_data/atac/peaks/elegans.{sample}/elegans.{sample}_peaks.xls',
    'data/external_data/atac/peaks/elegans.{sample}/elegans.{sample}_summits.bed',
    'data/external_data/atac/peaks/elegans.{sample}/elegans.{sample}_treat_pileup.bdg',
    temp('data/external_data/atac/peaks/elegans.{sample}/elegans.{sample}_treat_pileup.sorted.bdg'),
    'data/external_data/atac/peaks/elegans.{sample}/elegans.{sample}_treat_pileup.sorted.bw'
  params:
    'data/external_data/atac/peaks/elegans.{sample}'
  shell:
    '''
    macs2 callpeak --bdg --SPMR --gsize ce --nolambda --extsize 150 --shift -75 --nomodel -n elegans.{wildcards.sample} --outdir {params} -t {input[0]}
    sort -k 1,1 -k2,2n {output[3]} > {output[4]}
    bedGraphToBigWig {output[4]} {input[1]} {output[5]}
    '''

rule external_atac_yapc:
  input:
    'data/external_data/atac/peaks/elegans.{sample}.rep1/elegans.{sample}.rep1_treat_pileup.sorted.bw',
    'data/external_data/atac/peaks/elegans.{sample}.rep2/elegans.{sample}.rep2_treat_pileup.sorted.bw',
  output:
    'data/external_data/atac/yapc/elegans.{sample}.smooth_100_yapc_0.001.bed',
    'data/external_data/atac/yapc/elegans.{sample}.smooth_100_yapc_0.005.bed',
    'data/external_data/atac/yapc/elegans.{sample}.smooth_100_yapc_0.01.bed',
    'data/external_data/atac/yapc/elegans.{sample}.smooth_100_yapc_0.05.bed',
    'data/external_data/atac/yapc/elegans.{sample}.smooth_100_yapc_coverage.bw',
    'data/external_data/atac/yapc/elegans.{sample}.smooth_100_yapc_d2smooth.bw',
    'data/external_data/atac/yapc/elegans.{sample}.smooth_100_yapc_peaksall.bed',
    'data/external_data/atac/yapc/elegans.{sample}.smooth_100_yapc_peaksall.tsv',
    'data/external_data/atac/yapc/elegans.{sample}.smooth_100_yapc.tsv'
  params:
    'data/external_data/atac/yapc/elegans.{sample}.smooth_100_yapc'
  shell:
    '''
    yapc --smoothing-window-width 100 {params} elegans.{wildcards.sample} {input[0]} {input[1]}
    '''

rule lcap_coverage_m1m2_genes:
  input:
    'relmapping/lcap/trim20.bwa_pe_ce11.rm_unmapped_pe.rm_contigs_ce11.rm_q10/annot_ce_lcap_wt_ya_rep1.trim20.bwa_pe_ce11.rm_unmapped_pe.rm_contigs_ce11.rm_q10.bam',
    'relmapping/lcap/trim20.bwa_pe_ce11.rm_unmapped_pe.rm_contigs_ce11.rm_q10/annot_ce_lcap_wt_ya_rep2.trim20.bwa_pe_ce11.rm_unmapped_pe.rm_contigs_ce11.rm_q10.bam',
    'species/elegans/gene_annotation/elegans.gene_annotation.coding_gene_longest_isoform.bed',
    'promoter_annotation/promoters_gl_specific.elegans.m1m2.unique_promoter.genes',
    'promoter_annotation/promoters_gl_specific.elegans.no_m1m2.unique_promoter.genes',
    'species/elegans/genome/elegans.chrom.sizes.txt',
  output:
    temp('lcap_cov/elegans_lcap_wt_rep1.bed'),
    temp('lcap_cov/elegans_lcap_wt_rep2.bed'),
    'lcap_cov/elegans_lcap_wt.bed',
    'lcap_cov/elegans.gene_annotation.first_1000_exon.bed',
    'lcap_cov/promoters_gl_specific.elegans.m1m2.unique_promoter.cov.stranded.bed',
    'lcap_cov/promoters_gl_specific.elegans.no_m1m2.unique_promoter.cov.stranded.bed',
  resources:
    cpus=10
  shell:
    '''
    samtools view -b -f 80 -@ {resources.cpus} {input[0]} | bedtools bamtobed -i stdin > {output[0]}
    samtools view -b -f 64 -F 16 -@ {resources.cpus} {input[0]} | bedtools bamtobed -i stdin >> {output[0]}
    samtools view -b -f 80 -@ {resources.cpus} {input[1]} | bedtools bamtobed -i stdin > {output[1]}
    samtools view -b -f 64 -F 16 -@ {resources.cpus} {input[1]} | bedtools bamtobed -i stdin >> {output[1]}
    cat {output[0]} {output[1]} | sort -k 1,1 -k2,2n > {output[2]}
    python scripts/gene_exon_extract.py {input[2]} {output[3]} 1000
    sort -k 4,4 {output[3]} | join -1 1 -2 4 {input[3]} - | awk 'BEGIN{{OFS="\t";}}{{print $2, $3, $4, $1, $5, $6, $7, $8, $9, $10, $11, $12}}' | sort -k 1,1 -k2,2n | coverageBed -a stdin -b {output[2]} -g {input[5]} -sorted -split -S > {output[4]}
    sort -k 4,4 {output[3]} | join -1 1 -2 4 {input[4]} - | awk 'BEGIN{{OFS="\t";}}{{print $2, $3, $4, $1, $5, $6, $7, $8, $9, $10, $11, $12}}' | sort -k 1,1 -k2,2n | coverageBed -a stdin -b {output[2]} -g {input[5]} -sorted -split -S > {output[5]}
    Rscript 
    '''


rule ATAC_coverage_m1m2_GLRE:
  input:
    'RE_annotation/reg_elements_all.elegans.gl_specific.bed',
    'motif_enrichment/elegans/elegans.m1m2_clusters.bed',
    'relmapping/atac/tg_se.bwa_se_ce11.rm_unmapped.rm_contigs_ce11.rm_q10/annot_ce_atac_wt_ya_rep1.tg_se.bwa_se_ce11.rm_unmapped.rm_contigs_ce11.rm_q10.bam',
    'relmapping/atac/tg_se.bwa_se_ce11.rm_unmapped.rm_contigs_ce11.rm_q10/annot_ce_atac_wt_ya_rep2.tg_se.bwa_se_ce11.rm_unmapped.rm_contigs_ce11.rm_q10.bam',
    'species/elegans/genome/elegans.chrom.sizes.txt',
    'relmapping/atac/tg_se.bwa_se_ce11.rm_unmapped.rm_contigs_ce11.rm_q10/annot_ce_atac_glp1_ya_rep1.tg_se.bwa_se_ce11.rm_unmapped.rm_contigs_ce11.rm_q10.bam',
    'relmapping/atac/tg_se.bwa_se_ce11.rm_unmapped.rm_contigs_ce11.rm_q10/annot_ce_atac_glp1_ya_rep2.tg_se.bwa_se_ce11.rm_unmapped.rm_contigs_ce11.rm_q10.bam',
    'data/external_data/atac/aligned/elegans.pgc.rep1.bam',
    'data/external_data/atac/aligned/elegans.pgc.rep2.bam',
    'data/external_data/atac/aligned/elegans.adult_gl.rep1.bam',
    'data/external_data/atac/aligned/elegans.adult_gl.rep2.bam',
  output:
    'atac_cov/ele_m1m2_RE_atac_adult.txt',
    'atac_cov/ele_no_m1m2_RE_atac_adult.txt',
    'atac_cov/ele_m1m2_RE_atac_glp1.txt',
    'atac_cov/ele_no_m1m2_RE_atac_glp1.txt',
    'atac_cov/ele_m1m2_RE_atac_pgc.txt',
    'atac_cov/ele_no_m1m2_RE_atac_pgc.txt',
    'atac_cov/ele_m1m2_RE_atac_adult_gl.txt',
    'atac_cov/ele_no_m1m2_RE_atac_adult_gl.txt',
    'plots/ATAC_coverage_GLRE_elegans.pdf'
  shell:
    '''
    intersectBed -a {input[0]} -b {input[1]} -u | coverageBed -a stdin -b {input[2]} {input[3]} -g {input[4]} -sorted > {output[0]}
    intersectBed -a {input[0]} -b {input[1]} -v | coverageBed -a stdin -b {input[2]} {input[3]} -g {input[4]} -sorted >	{output[1]}
    intersectBed -a {input[0]} -b {input[1]} -u | coverageBed -a stdin -b {input[5]} {input[6]} -g {input[4]} -sorted > {output[2]}
    intersectBed -a {input[0]} -b {input[1]} -v | coverageBed -a stdin -b {input[5]} {input[6]} -g {input[4]} -sorted > {output[3]}
    intersectBed -a {input[0]} -b {input[1]} -u | coverageBed -a stdin -b {input[7]} {input[8]} -g {input[4]} -sorted > {output[4]}
    intersectBed -a {input[0]} -b {input[1]} -v | coverageBed -a stdin -b {input[7]} {input[8]} -g {input[4]} -sorted > {output[5]}
    intersectBed -a {input[0]} -b {input[1]} -u | coverageBed -a stdin -b {input[9]} {input[10]} -g {input[4]} -sorted > {output[6]}
    intersectBed -a {input[0]} -b {input[1]} -v | coverageBed -a stdin -b {input[9]} {input[10]} -g {input[4]} -sorted > {output[7]}
    Rscript scripts/ATAC_coverage_GLRE.R
    '''


