rule CERP2_CELE2_hmm_fetch:
  input:
    'software/dfamscan/caenorhabditis_elegans_dfam.hmm',
  output:
    'data/external_data/CERP2.hmm',
    'data/external_data/CELE2.hmm',
  shell:
    '''
    hmmfetch {input[0]} CERP2 > {output[0]}
    hmmfetch {input[0]} CELE2 > {output[1]}
    '''


rule CERP2_CELE2_re_annotation_elegans_hmm:
  input:
    'species/elegans/repeats_dfam/elegans_dfam.repeats.id.bed',
    'species/elegans/genome/elegans.fa',
    'data/external_data/CERP2.hmm',
    'data/external_data/CELE2.hmm',
  output:
    'CELE2_CERP2_other_species/elegans.CERP2.dfam.fa',
    'CELE2_CERP2_other_species/elegans.CELE2.dfam.fa',
    'CELE2_CERP2_other_species/elegans.CERP2.dfam.hmmalign',
    'CELE2_CERP2_other_species/elegans.CELE2.dfam.hmmalign',
    'CELE2_CERP2_other_species/elegans.CERP2.dfam.hmmalign.hmm',
    'CELE2_CERP2_other_species/elegans.CELE2.dfam.hmmalign.hmm',
  shell:
    '''
    grep CERP2 {input[0]} | fastaFromBed -fi {input[1]} -bed stdin -fo {output[0]}
    grep CELE2 {input[0]} | fastaFromBed -fi {input[1]} -bed stdin -fo {output[0]}
    hmmalign --trim -o {output[2]} {input[2]} {output[0]}
    hmmalign --trim -o {output[3]} {input[3]} {output[1]}
    hmmbuild {output[4]} {output[2]}
    hmmbuild {output[5]} {output[3]}
    '''


rule CERP2_CELE2_re_annotation_all_species:
  input:
    'CELE2_CERP2_other_species/elegans.CERP2.dfam.hmmalign.hmm',
    'CELE2_CERP2_other_species/elegans.CELE2.dfam.hmmalign.hmm',
    'species/{sample}/genome/{sample}.fa',
    'motif_enrichment/{sample}/{sample}.m1m2_clusters.bed',
    'motif_enrichment/{sample}/{sample}.m1.bed',
    'motif_enrichment/{sample}/{sample}.m2.bed',
  output:
    'CELE2_CERP2_other_species/{sample}.CERP2.hmmalign.out',
    'CELE2_CERP2_other_species/{sample}.CELE2.hmmalign.out',
    'CELE2_CERP2_other_species/{sample}.CERP2.hmmalign.bed',
    'CELE2_CERP2_other_species/{sample}.CELE2.hmmalign.bed',
    'CELE2_CERP2_other_species/{sample}.CERP2.hmmalign.nooverlap.bed',
    'CELE2_CERP2_other_species/{sample}.CELE2.hmmalign.nooverlap.bed',
    'CELE2_CERP2_other_species/{sample}.summary',
  resources:
    cpus=10
  shell:
    '''
    nhmmer --tblout {output[0]} --noali --cpu {resources.cpus} {input[0]} {input[2]}
    nhmmer --tblout {output[1]} --noali --cpu {resources.cpus} {input[1]} {input[2]}
    sed '1,3d' {output[0]} | grep -v "#" | awk 'BEGIN{{OFS="\t";}}{{if ($13 < 0.001 && $7 < $8) print $1, $7, $8, $3, $13, $12; else if ($13 < 0.001 && $7 > $8) print $1, $8, $7, $3, $13, $12}}' > {output[2]}
    sed '1,3d' {output[1]} | grep -v "#" | awk 'BEGIN{{OFS="\t";}}{{if ($13 < 0.001 && $7 < $8) print $1, $7, $8, $3, $13, $12; else if ($13 < 0.001 && $7 > $8) print $1, $8, $7, $3, $13, $12}}' > {output[3]}
    intersectBed -a {output[2]} -b {output[3]} -wa -wb | awk 'BEGIN{{OFS="\t";}}{{if ($5 > $11) print $1, $2, $3, $4, $5}}' | intersectBed -a {output[2]} -b stdin -v | sort -k 1,1 -k2,2n | mergeBed -i stdin > {output[4]}
    intersectBed -a {output[3]}	-b {output[2]} -wa -wb | awk 'BEGIN{{OFS="\t";}}{{if ($5 > $11) print $1, $2, $3, $4, $5}}' | intersectBed -a {output[3]} -b stdin -v | sort -k 1,1 -k2,2n | mergeBed -i stdin > {output[5]}
    CERP2_n=$(wc -l < {output[4]})
    CELE2_n=$(wc -l < {output[5]})
    echo -e "CERP2\t"$CERP2_n > {output[6]}
    for orientation in convergent divergent tandem_m1m2 tandem_m2m1
    do
    m1m2_peak=$(grep $orientation {input[3]} | intersectBed -a {output[4]} -b stdin -u | wc -l)
    echo -e $orientation"\t"$m1m2_peak >> {output[6]}
    done
    m1_m2_peak=$(intersectBed -a {output[4]} -b {input[4]} {input[5]} -u | wc -l)
    echo -e "any_m1_m2\t"$m1_m2_peak >> {output[6]}
    echo -e "CELE2\t"$CELE2_n >> {output[6]}
    for	orientation in convergent divergent tandem_m1m2	tandem_m2m1
    do
    m1m2_peak=$(grep $orientation {input[3]} | intersectBed -a {output[5]} -b stdin -u | wc -l)
    echo -e $orientation"\t"$m1m2_peak >> {output[6]}
    done
    m1_m2_peak=$(intersectBed -a {output[5]} -b {input[4]} {input[5]} -u | wc -l)
    echo -e "any_m1_m2\t"$m1_m2_peak >>	{output[6]}
    '''


rule dfam_intact_CERP2_CELE2:
  input:
    'species/elegans/repeats_dfam/elegans_dfam.repeats.id.bed',
    'motif_enrichment/elegans/elegans.m1m2_clusters.bed',
    'species/elegans/genome/elegans.fa',
    'software/dfamscan/caenorhabditis_elegans_dfam.hmm',
  output:
    temp('intact_repeats/merged_m1m2_elegans.bed'),
    'intact_repeats/intact_CERP2_elegans.bed',
    'intact_repeats/intact_CERP2_elegans.fa',
    'intact_repeats/intact_CERP2_elegans.fa.CERP2_aln',
    'intact_repeats/intact_CERP2_elegans.hmm',
    'software/dfamscan/caenorhabditis_elegans_dfam.CERP2.hmm',
    'intact_repeats/intact_CELE2_elegans.bed',
    'intact_repeats/intact_CELE2_elegans.fa',
    'intact_repeats/intact_CELE2_elegans.fa.CELE2_aln',
    'intact_repeats/intact_CELE2_elegans.hmm',
    'software/dfamscan/caenorhabditis_elegans_dfam.CELE2.hmm',
  shell:
    '''
    grep divergent {input[1]} | sort -k 1,1 -k2,2n | mergeBed -i stdin -c 4,5,6 -o distinct,distinct,distinct > {output[0]}
    grep CERP2 {input[0]} | intersectBed -a stdin -b {output[0]} -wa -wb | sort -k 1,1 -k2,2n | mergeBed -i stdin -c 12 -o distinct | awk '{{if ($4 ~ ",") print}}' > {output[1]}
    fastaFromBed -bed {output[1]} -fi {input[2]} -fo {output[2]}
    hmmfetch -o {output[5]} {input[3]} CERP2
    hmmalign -o {output[3]} {output[5]} {output[2]}
    hmmbuild {output[4]} {output[3]}
    grep tandem_m2m1 {input[1]} | sort -k 1,1 -k2,2n | mergeBed -i stdin -c 4,5,6 -o distinct,distinct,distinct > {output[0]}
    grep CELE2 {input[0]} | intersectBed -a stdin -b {output[0]} -wa -wb | sort -k 1,1 -k2,2n | mergeBed -i stdin -c 12 -o distinct | awk '{{if ($4 ~ ",") print}}' > {output[6]}
    fastaFromBed -bed {output[6]} -fi {input[2]} -fo {output[7]}
    hmmfetch -o {output[10]} {input[3]} CELE2
    hmmalign -o {output[8]} {output[10]} {output[7]}
    hmmbuild {output[9]} {output[8]}
    '''


rule dfam_active_CERP2:
  input:
    'species/elegans/repeats_dfam/elegans_dfam.repeats.id.bed',
    'RE_annotation/reg_elements_all.elegans.gl_specific.bed',
    'motif_enrichment/elegans/elegans.m1m2_clusters.bed',
    'species/elegans/genome/elegans.fa',
    'RE_annotation/reg_elements_all.elegans.bed',
    'intact_repeats/intact_CERP2_elegans.hmm',
  output:
    'intact_repeats/GL_active_CERP2_elegans.fa',
    'intact_repeats/inactive_CERP2_elegans.fa',
    'intact_repeats/GL_active_CERP2_elegans.fa.CERP2_aln',
    'intact_repeats/GL_active_CERP2_elegans.hmm',
    'intact_repeats/inactive_CERP2_elegans.fa.CERP2_aln',
    'intact_repeats/inactive_CERP2_elegans.hmm',
  shell:
    '''
    grep CERP2 {input[0]} | intersectBed -a stdin -b {input[1]} -u | intersectBed -a stdin -b {input[2]} -wa -wb | awk 'BEGIN{{OFS="\t";}}{{print $1, $2, $3, $4, $5, $12}}' | sort -k 1,1 -k2,2n | fastaFromBed -bed stdin -fi {input[3]} -fo {output[0]} -s
    grep CERP2 {input[0]} | intersectBed -a stdin -b {input[1]} -v | intersectBed -a stdin -b {input[2]} -u | sort -k 1,1 -k2,2n | mergeBed -i stdin | fastaFromBed -bed stdin -fi {input[3]} -fo {output[1]}
    hmmalign --trim -o {output[2]} {input[5]} {output[0]}
    hmmbuild {output[3]} {output[2]}
    hmmalign --trim -o {output[4]} {input[5]} {output[1]}
    hmmbuild {output[5]} {output[4]}
    '''


rule other_species_CERP2:
  input:
    'motif_enrichment/{sample}/{sample}.m1m2_clusters.bed',
    'species/{sample}/genome/{sample}.chrom.sizes.txt',
    'species/{sample}/genome/{sample}.fa',
    'intact_repeats/intact_CERP2_elegans.hmm',
  output:
    'intact_repeats/{sample}/{sample}_CERP2_all.fa',
    'intact_repeats/{sample}/{sample}_CERP2_all.clw',
    'intact_repeats/{sample}/{sample}_CERP2_all.clw.hmm',
  shell:
    '''
    grep divergent {input[0]} | awk 'BEGIN{{OFS="\t";}}{{if ($4 ~ "14bp" || $4 ~ "15bp" || $4 ~ "16bp" || $4 ~ "17bp" || $4 ~ "18bp" || $4 ~ "19bp") print}}' | slopBed -b 100 -s -i stdin -g {input[1]} | fastaFromBed -fi {input[2]} -bed stdin -fo {output[0]} -s
    muscle -in {output[0]} -clw -out {output[1]}
    hmmbuild --fragthresh 0 {output[2]} {output[1]}
    '''


rule elegans_specific_CERP2:
  input:
    'motif_enrichment/elegans/elegans.m1m2_clusters.bed',
    'species/elegans/genome/elegans.chrom.sizes.txt',
    'species/elegans/genome/elegans.fa',
    'RE_annotation/reg_elements_all.elegans.gl_specific.bed',
    'RE_annotation/reg_elements_all.elegans.bed',
    'intact_repeats/intact_CERP2_elegans.bed',
  output:
    'intact_repeats/elegans/elegans_CERP2.GLRE.fa',
    'intact_repeats/elegans/elegans_CERP2.inactive.fa',
    'intact_repeats/elegans/elegans_CERP2.GLRE.clw',
    'intact_repeats/elegans/elegans_CERP2.GLRE.clw.hmm',
    'intact_repeats/elegans/elegans_CERP2.inactive.clw',
    'intact_repeats/elegans/elegans_CERP2.inactive.clw.hmm',
    'intact_repeats/elegans/elegans_CERP2.intact.fa',
    'intact_repeats/elegans/elegans_CERP2.intact.clw',
    'intact_repeats/elegans/elegans_CERP2.intact.clw.hmm',
  shell:
    '''
    grep divergent {input[0]} | awk 'BEGIN{{OFS="\t";}}{{if ($4 ~ "14bp" || $4 ~ "15bp" || $4 ~ "16bp" || $4 ~ "17bp" || $4 ~ "18bp" || $4 ~ "19bp") print}}' | intersectBed -a stdin -b {input[3]} -u | slopBed -b 100 -s -i stdin -g {input[1]} | fastaFromBed -fi {input[2]} -bed stdin -fo {output[0]} -s
    grep divergent {input[0]} | awk 'BEGIN{{OFS="\t";}}{{if ($4 ~ "14bp" || $4 ~ "15bp" || $4 ~ "16bp" || $4 ~ "17bp" || $4 ~ "18bp" || $4 ~ "19bp") print}}' | intersectBed -a stdin -b {input[4]} -v | slopBed -b 100 -s -i stdin -g {input[1]} | fastaFromBed -fi {input[2]} -bed stdin -fo {output[1]} -s
    muscle -in {output[0]} -clw -out {output[2]}
    hmmbuild --fragthresh 0 {output[3]} {output[2]}
    muscle -in {output[1]} -clw -out {output[4]}
    hmmbuild --fragthresh 0 {output[5]} {output[4]}
    grep divergent {input[0]} | awk 'BEGIN{{OFS="\t";}}{{if ($4 ~ "14bp" || $4 ~ "15bp" || $4 ~ "16bp" || $4 ~ "17bp" || $4 ~ "18bp" || $4 ~ "19bp") print}}' | intersectBed -a stdin -b {input[5]} -u | slopBed -b 100 -s -i stdin -g {input[1]} | fastaFromBed -fi {input[2]} -bed stdin -fo {output[6]} -s
    muscle -in {output[6]} -clw -out {output[7]}
    hmmbuild --fragthresh 0 {output[8]} {output[7]}
    '''


rule HIM17_bound_motifs:
  input:
    'motif_enrichment/elegans/elegans.m1m2_clusters.bed',
    'species/elegans/genome/elegans.chrom.sizes.txt',
    'species/elegans/genome/elegans.fa',
    'data/elegans/chipseq/yapc/elegans_HIM17.smooth_100_yapc_0.00001.bed',
  output:
    'intact_repeats/elegans/elegans_CERP2.HIM17.fa',
    'intact_repeats/elegans/elegans_CERP2.no_HIM17.fa',
    'intact_repeats/elegans/elegans_CERP2.HIM17.clw',
    'intact_repeats/elegans/elegans_CERP2.HIM17.clw.hmm',
    'intact_repeats/elegans/elegans_CERP2.no_HIM17.clw',
    'intact_repeats/elegans/elegans_CERP2.no_HIM17.clw.hmm',
  shell:
    '''
    grep divergent {input[0]} | awk 'BEGIN{{OFS="\t";}}{{if ($4 ~ "14bp" || $4 ~ "15bp" || $4 ~ "16bp" || $4 ~ "17bp" || $4 ~ "18bp" || $4 ~ "19bp") print}}' | intersectBed -a stdin -b {input[3]} -u | slopBed -b 100 -s -i stdin -g {input[1]} | fastaFromBed -fi {input[2]} -bed stdin -fo {output[0]} -s
    grep divergent {input[0]} | awk 'BEGIN{{OFS="\t";}}{{if ($4 ~ "14bp" || $4 ~ "15bp" || $4 ~ "16bp" || $4 ~ "17bp" || $4 ~ "18bp" || $4 ~ "19bp") print}}' | intersectBed -a stdin -b {input[3]} -v | slopBed -b 100 -s -i stdin -g {input[1]} | fastaFromBed -fi {input[2]} -bed stdin -fo {output[1]} -s
    muscle -in {output[0]} -clw -out {output[2]}
    hmmbuild --fragthresh 0 {output[3]} {output[2]}
    muscle -in {output[1]} -clw -out {output[4]}
    hmmbuild --fragthresh 0 {output[5]} {output[4]}
    '''


rule elegans_specific_CELE2:
  input:
    'motif_enrichment/elegans/elegans.m1m2_clusters.bed',
    'species/elegans/genome/elegans.chrom.sizes.txt',
    'species/elegans/genome/elegans.fa',
    'RE_annotation/reg_elements_all.elegans.gl_specific.bed',
    'RE_annotation/reg_elements_all.elegans.bed',
    'intact_repeats/intact_CELE2_elegans.bed',
  output:
    'intact_repeats/elegans/elegans_CELE2.GLRE.fa',
    'intact_repeats/elegans/elegans_CELE2.inactive.fa',
    'intact_repeats/elegans/elegans_CELE2.GLRE.clw',
    'intact_repeats/elegans/elegans_CELE2.GLRE.clw.hmm',
    'intact_repeats/elegans/elegans_CELE2.inactive.clw',
    'intact_repeats/elegans/elegans_CELE2.inactive.clw.hmm',
    'intact_repeats/elegans/elegans_CELE2.intact.fa',
    'intact_repeats/elegans/elegans_CELE2.intact.clw',
    'intact_repeats/elegans/elegans_CELE2.intact.clw.hmm',
  shell:
    '''
    grep tandem_m2m1 {input[0]} | awk 'BEGIN{{OFS="\t";}}{{if ($4 ~ "24bp" || $4 ~ "25bp" || $4 ~ "26bp" || $4 ~ "27bp" || $4 ~ "28bp" || $4 ~ "29bp") print}}' | intersectBed -a stdin -b {input[3]} -u | slopBed -b 100 -s -i stdin -g {input[1]} | fastaFromBed -fi {input[2]} -bed stdin -fo {output[0]} -s
    grep tandem_m2m1 {input[0]} | awk 'BEGIN{{OFS="\t";}}{{if ($4 ~ "24bp" || $4 ~ "25bp" || $4 ~ "26bp" || $4 ~ "27bp" || $4 ~ "28bp" || $4 ~ "29bp") print}}' | intersectBed -a stdin -b {input[4]} -v | slopBed -b 100 -s -i stdin -g {input[1]} | fastaFromBed -fi {input[2]} -bed stdin -fo {output[1]} -s
    muscle -in {output[0]} -clw -out {output[2]}
    hmmbuild --fragthresh 0 {output[3]} {output[2]}
    muscle -in {output[1]} -clw -out {output[4]}
    hmmbuild --fragthresh 0 {output[5]} {output[4]}
    grep tandem_m2m1 {input[0]} | awk 'BEGIN{{OFS="\t";}}{{if ($4 ~ "24bp" || $4 ~ "25bp" || $4 ~ "26bp" || $4 ~ "27bp" || $4 ~ "28bp" || $4 ~ "29bp") print}}' | intersectBed -a stdin -b {input[5]} -u | slopBed -b 100 -s -i stdin -g {input[1]} | fastaFromBed -fi {input[2]} -bed stdin -fo {output[6]} -s
    muscle -in {output[6]} -clw -out {output[7]}
    hmmbuild --fragthresh 0 {output[8]} {output[7]}
    '''


rule TT_periodicity_elegans_repeats:
  input:
    'motif_enrichment/elegans/elegans.m1m2_clusters.bed',
    'RE_annotation/reg_elements_all.elegans.gl_specific.bed',
    'species/elegans/genome/elegans.chrom.sizes.txt',
    'species/elegans/genome/elegans.fa',
    'RE_annotation/reg_elements_all.elegans.bed',
    'intact_repeats/intact_CERP2_elegans.bed',
    'species/elegans/repeats_dfam/elegans_dfam.repeats.id.bed',
    'intact_repeats/intact_CELE2_elegans.bed',
  output:
    temp('intact_repeats/elegans/elegans_divergent_m1m2.bed'),
    'intact_repeats/elegans/elegans_divergent_m1m2.GLRE.for_clustering.fa',
    'intact_repeats/elegans/elegans_divergent_m1m2.GLRE.for_clustering.TT_freq',
    'intact_repeats/elegans/elegans_divergent_m1m2.inactive.for_clustering.fa',
    'intact_repeats/elegans/elegans_divergent_m1m2.inactive.for_clustering.TT_freq',
    'intact_repeats/elegans/elegans_divergent_m1m2.intact.for_clustering.fa',
    'intact_repeats/elegans/elegans_divergent_m1m2.intact.for_clustering.TT_freq',
    'intact_repeats/elegans/elegans_CERP2.GLRE.for_clustering.fa',
    'intact_repeats/elegans/elegans_CERP2.GLRE.for_clustering.TT_freq',
    'intact_repeats/elegans/elegans_CERP2.inactive.for_clustering.fa',
    'intact_repeats/elegans/elegans_CERP2.inactive.for_clustering.TT_freq',
    temp('intact_repeats/elegans/elegans_tandem_m2m1.bed'),
    'intact_repeats/elegans/elegans_tandem_m2m1.GLRE.for_clustering.fa',
    'intact_repeats/elegans/elegans_tandem_m2m1.GLRE.for_clustering.TT_freq',
    'intact_repeats/elegans/elegans_tandem_m2m1.inactive.for_clustering.fa',
    'intact_repeats/elegans/elegans_tandem_m2m1.inactive.for_clustering.TT_freq',
    'intact_repeats/elegans/elegans_tandem_m2m1.intact.for_clustering.fa',
    'intact_repeats/elegans/elegans_tandem_m2m1.intact.for_clustering.TT_freq',
    'intact_repeats/elegans/elegans_CELE2.GLRE.for_clustering.fa',
    'intact_repeats/elegans/elegans_CELE2.GLRE.for_clustering.TT_freq',
    'intact_repeats/elegans/elegans_CELE2.inactive.for_clustering.fa',
    'intact_repeats/elegans/elegans_CELE2.inactive.for_clustering.TT_freq',
  shell:
    '''
    grep divergent {input[0]} | awk 'BEGIN{{OFS="\t";}}{{if ($4 ~ "14bp" || $4 ~ "15bp" || $4 ~ "16bp" || $4 ~ "17bp" || $4 ~ "18bp" || $4 ~ "19bp") print $1, $2, $3, $4, 0, $6}}' > {output[0]}
    intersectBed -a {output[0]} -b {input[1]} -u | slopBed -l 0 -r 200 -s -i stdin -g {input[2]} | fastaFromBed -fi {input[3]} -bed stdin -fo {output[1]} -s
    python scripts/dinuleotide_count.py {output[1]} {output[2]} TT
    intersectBed -a {output[0]} -b {input[4]} -v | slopBed -l 0 -r 200 -s -i stdin -g {input[2]} | fastaFromBed -fi {input[3]} -bed stdin -fo {output[3]} -s
    python scripts/dinuleotide_count.py {output[3]} {output[4]} TT
    intersectBed -a {output[0]} -b {input[5]} -u | slopBed -l 0 -r 200 -s -i stdin -g {input[2]} | fastaFromBed -fi {input[3]} -bed stdin -fo {output[5]} -s
    python scripts/dinuleotide_count.py {output[5]} {output[6]} TT
    grep CERP2 {input[6]} | intersectBed -a {output[0]} -b stdin -u | intersectBed -a stdin -b {input[1]} -u | slopBed -l 0 -r 200 -s -i stdin -g {input[2]} | fastaFromBed -fi {input[3]} -bed stdin -fo {output[7]} -s
    python scripts/dinuleotide_count.py {output[7]} {output[8]} TT
    grep CERP2 {input[6]} | intersectBed -a {output[0]} -b stdin -u | intersectBed -a stdin -b {input[4]} -v | slopBed -l 0 -r 200 -s -i stdin -g {input[2]} | fastaFromBed -fi {input[3]} -bed stdin -fo {output[9]} -s
    python scripts/dinuleotide_count.py {output[9]} {output[10]} TT
    grep tandem_m2m1 {input[0]} | awk 'BEGIN{{OFS="\t";}}{{if ($4 ~ "24bp" || $4 ~ "25bp" || $4 ~ "26bp" || $4 ~ "27bp" || $4 ~ "28bp" || $4 ~ "29bp") print $1, $2, $3, $4, 0, $6}}' > {output[11]}
    intersectBed -a {output[11]} -b {input[1]} -u | slopBed -l 0 -r 200 -s -i stdin -g {input[2]} | fastaFromBed -fi {input[3]} -bed stdin -fo {output[12]} -s
    python scripts/dinuleotide_count.py {output[12]} {output[13]} TT
    intersectBed -a {output[11]} -b {input[4]} -v | slopBed -l 0 -r 200 -s -i stdin -g {input[2]} | fastaFromBed -fi {input[3]} -bed stdin -fo {output[14]} -s
    python scripts/dinuleotide_count.py {output[14]} {output[15]} TT
    intersectBed -a {output[11]} -b {input[7]} -u | slopBed -l 0 -r 200 -s -i stdin -g {input[2]} | fastaFromBed -fi {input[3]} -bed stdin -fo {output[16]} -s
    python scripts/dinuleotide_count.py {output[16]} {output[17]} TT
    grep CELE2 {input[6]} | intersectBed -a {output[11]} -b stdin -u | intersectBed -a stdin -b {input[1]} -u | slopBed -l 0 -r 200 -s -i stdin -g {input[2]} | fastaFromBed -fi {input[3]} -bed stdin -fo {output[18]} -s
    python scripts/dinuleotide_count.py {output[18]} {output[19]} TT
    grep CELE2 {input[6]} | intersectBed -a {output[11]} -b stdin -u | intersectBed -a stdin -b {input[4]} -v | slopBed -l 0 -r 200 -s -i stdin -g {input[2]} | fastaFromBed -fi {input[3]} -bed stdin -fo {output[20]} -s
    python scripts/dinuleotide_count.py {output[20]} {output[21]} TT
    '''


rule TT_plots:
  input:
    'intact_repeats/elegans/elegans_divergent_m1m2.GLRE.for_clustering.TT_freq',
    'intact_repeats/elegans/elegans_divergent_m1m2.inactive.for_clustering.TT_freq',
    'intact_repeats/elegans/elegans_divergent_m1m2.intact.for_clustering.TT_freq',
    'intact_repeats/elegans/elegans_CERP2.GLRE.for_clustering.TT_freq',
    'intact_repeats/elegans/elegans_CERP2.inactive.for_clustering.TT_freq',
    'intact_repeats/elegans/elegans_tandem_m2m1.GLRE.for_clustering.TT_freq',
    'intact_repeats/elegans/elegans_tandem_m2m1.inactive.for_clustering.TT_freq',
    'intact_repeats/elegans/elegans_tandem_m2m1.intact.for_clustering.TT_freq',
    'intact_repeats/elegans/elegans_CELE2.GLRE.for_clustering.TT_freq',
    'intact_repeats/elegans/elegans_CELE2.inactive.for_clustering.TT_freq',
  output:
    'plots/TT_periodicity_CERP2_elements.pdf',
    'plots/TT_periodicity_divergent_m1m2_elements.pdf',
    'plots/TT_periodicity_CELE2_elements.pdf',
    'plots/TT_periodicity_tandem_m2m1_elements.pdf',
  shell:
    '''
    Rscript scripts/nucleosomal_TT.R
    '''


rule dfam_repeat_consensus:
  input:
    'software/dfamscan/caenorhabditis_elegans_dfam.hmm',
  output:
    'data/external_data/dfam_repeats.fa',
    'data/external_data/CERP2.fa',
    'data/external_data/CELE2.fa',
  shell:
    '''
    hmmemit -c {input} > {output[0]}
    samtools faidx {output[0]} CERP2-consensus | awk -F"-" '{{print $1}}' > {output[1]}
    samtools faidx {output[0]} CELE2-consensus | awk -F"-" '{{print $1}}' > {output[2]}
    '''
