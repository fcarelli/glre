rule genome_download_main_species:
  output:
    elegans='species/elegans/genome/elegans.fa',
    inopinata='species/inopinata/genome/inopinata.fa',
    briggsae='species/briggsae/genome/briggsae.fa',
    remanei='species/remanei/genome/remanei.fa',
    nigoni='species/nigoni/genome/nigoni.fa',
    contortus='species/contortus/genome/contortus.fa',
    tipulae='species/tipulae/genome/tipulae.fa',
    bacteriophora='species/bacteriophora/genome/bacteriophora.fa',
    ceylanicum='species/ceylanicum/genome/ceylanicum.fa',
    pacificus='species/pacificus/genome/pacificus.fa',
    redivivus='species/redivivus/genome/redivivus.fa',
    bovis='species/bovis/genome/bovis.fa',
    monodelphis='species/monodelphis/genome/monodelphis.fa',
    becei='species/becei/genome/becei.fa',
    quiockensis='species/quiockensis/genome/quiockensis.fa',
    plicata='species/plicata/genome/plicata.fa',
    uteleia='species/uteleia/genome/uteleia.fa',
    exspectatus='species/exspectatus/genome/exspectatus.fa',
  shell:
    '''
    wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS275/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS275.genomic.fa.gz -O {output.elegans}.gz; gunzip {output.elegans}.gz
    wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS275/species/c_inopinata/PRJDB5687/c_inopinata.PRJDB5687.WS275.genomic.fa.gz -O {output.inopinata}.gz; gunzip {output.inopinata}.gz
    wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS275/species/c_briggsae/PRJNA10731/c_briggsae.PRJNA10731.WS275.genomic.fa.gz -O {output.briggsae}.gz; gunzip {output.briggsae}.gz
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/010/183/535/GCA_010183535.1_CRPX506/GCA_010183535.1_CRPX506_genomic.fna.gz -O {output.remanei}.gz; gunzip {output.remanei}.gz
    wget ftp://ftp.wormbase.org/pub/wormbase/releases/WS275/species/c_nigoni/PRJNA384657/c_nigoni.PRJNA384657.WS275.genomic.fa.gz -O {output.nigoni}.gz; gunzip {output.nigoni}.gz
    wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species/haemonchus_contortus/PRJEB506/haemonchus_contortus.PRJEB506.WBPS14.genomic.fa.gz -O {output.contortus}.gz; gunzip {output.contortus}.gz
    wget ftp://ftp.wormbase.org/pub/wormbase/species/o_tipulae/sequence/genomic/o_tipulae.PRJEB15512.WS275.genomic.fa.gz -O {output.tipulae}.gz; gunzip {output.tipulae}.gz
    wget ftp://ftp.wormbase.org/pub/wormbase/species/h_bacteriophora/sequence/genomic/h_bacteriophora.PRJNA13977.WS248.genomic.fa.gz -O {output.bacteriophora}.gz; gunzip {output.bacteriophora}.gz
    wget ftp://ftp.wormbase.org/pub/wormbase/species/a_ceylanicum/PRJNA231479/sequence/genomic/a_ceylanicum.PRJNA231479.WS248.genomic.fa.gz -O {output.ceylanicum}.gz; gunzip {output.ceylanicum}.gz
    wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species/pristionchus_pacificus/PRJNA12644/pristionchus_pacificus.PRJNA12644.WBPS14.genomic.fa.gz -O {output.pacificus}.gz; gunzip {output.pacificus}.gz
    wget ftp://ftp.wormbase.org/pub/wormbase/species/p_redivivus/sequence/genomic/p_redivivus.PRJNA186477.WS275.genomic.fa.gz -O {output.redivivus}.gz; gunzip {output.redivivus}.gz
    wget http://download.caenorhabditis.org/v1/sequence/Caenorhabditis_bovis_LS_v1.scaffolds.fa.gz -O {output.bovis}.gz; gunzip {output.bovis}.gz
    wget http://download.caenorhabditis.org/v1/sequence/Caenorhabditis_monodelphis_JU1667_v1.scaffolds.fa.gz -O {output.monodelphis}.gz; gunzip {output.monodelphis}.gz
    wget http://download.caenorhabditis.org/v1/sequence/Caenorhabditis_sp29_QG2083_v1.scaffolds.fa.gz -O {output.becei}.gz; gunzip {output.becei}.gz
    wget http://download.caenorhabditis.org/v1/sequence/Caenorhabditis_sp38_JU2809_v1.scaffolds.fa.gz -O {output.quiockensis}.gz; gunzip {output.quiockensis}.gz
    wget http://download.caenorhabditis.org/v1/sequence/Caenorhabditis_plicata_SB355_v1.scaffolds.fa.gz -O {output.plicata}.gz; gunzip {output.plicata}.gz
    wget http://download.caenorhabditis.org/v1/sequence/Caenorhabditis_sp31_JU2585_v1.scaffolds.fa.gz -O {output.uteleia}.gz; gunzip {output.uteleia}.gz
    wget http://pristionchus.org/variome/exspectatus_genome.fa.gz -O {output.exspectatus}.gz; gunzip {output.exspectatus}.gz
    '''


rule genome_idx:
  input:
    'species/{sample}/genome/{sample}.fa'
  output:
    fai_idx = 'species/{sample}/genome/{sample}.fa.fai',
    chr_len = 'species/{sample}/genome/{sample}.chrom.sizes.txt'
  shell:
    '''
    samtools faidx {input}
    cut -f 1,2 {output.fai_idx} | sort -k 1,1 > {output.chr_len}
    '''


rule gene_annotation_download_all:
  output:
    elegans='species/elegans/gene_annotation/elegans.annotations.all.gff3.gz',
    inopinata='species/inopinata/gene_annotation/inopinata.annotations.all.gff3.gz',
    briggsae='species/briggsae/gene_annotation/briggsae.annotations.all.gff3.gz',
    remanei='species/remanei/gene_annotation/remanei.annotations.all.gff3.gz',
    nigoni='species/nigoni/gene_annotation/nigoni.annotations.all.gff3.gz',
    contortus='species/contortus/gene_annotation/contortus.annotations.all.gff3.gz',
    tipulae='species/tipulae/gene_annotation/tipulae.annotations.all.gff3.gz',
    bacteriophora='species/bacteriophora/gene_annotation/bacteriophora.annotations.all.gff3.gz',
    ceylanicum='species/ceylanicum/gene_annotation/ceylanicum.annotations.all.gff3.gz',
    pacificus='species/pacificus/gene_annotation/pacificus.annotations.all.gff3.gz',
    redivivus='species/redivivus/gene_annotation/redivivus.annotations.all.gff3.gz',
  shell:
    '''
    wget ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/gff/c_elegans.PRJNA13758.WS275.annotations.gff3.gz -O {output.elegans}
    wget ftp://ftp.wormbase.org/pub/wormbase/species/c_inopinata/gff/c_inopinata.PRJDB5687.WS275.annotations.gff3.gz -O {output.inopinata}
    wget ftp://ftp.wormbase.org/pub/wormbase/species/c_briggsae/gff/c_briggsae.PRJNA10731.WS275.annotations.gff3.gz -O {output.briggsae}
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/010/183/535/GCA_010183535.1_CRPX506/GCA_010183535.1_CRPX506_genomic.gff.gz -O {output.remanei}
    wget ftp://ftp.wormbase.org/pub/wormbase/species/c_nigoni/gff/c_nigoni.PRJNA384657.WS275.annotations.gff3.gz -O {output.nigoni}
    wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species/haemonchus_contortus/PRJEB506/haemonchus_contortus.PRJEB506.WBPS14.annotations.gff3.gz -O {output.contortus}
    wget ftp://ftp.wormbase.org/pub/wormbase/species/o_tipulae/gff/o_tipulae.PRJEB15512.WS275.annotations.gff3.gz -O {output.tipulae}
    wget ftp://ftp.wormbase.org/pub/wormbase/species/h_bacteriophora/gff/h_bacteriophora.PRJNA13977.WS248.annotations.gff3.gz -O {output.bacteriophora}
    wget ftp://ftp.wormbase.org/pub/wormbase/species/a_ceylanicum/PRJNA231479/gff/a_ceylanicum.PRJNA231479.WS248.annotations.gff3.gz -O {output.ceylanicum}
    wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species/pristionchus_pacificus/PRJNA12644/pristionchus_pacificus.PRJNA12644.WBPS14.annotations.gff3.gz -O {output.pacificus}
    wget ftp://ftp.wormbase.org/pub/wormbase/species/p_redivivus/gff/p_redivivus.PRJNA186477.WS275.annotations.gff3.gz -O {output.redivivus}
    '''


rule gene_annotation_processing:
  input:
    gff = "species/{sample}/gene_annotation/{sample}.annotations.all.gff3.gz",
    genome = "species/{sample}/genome/{sample}.fa"
  output:
    all_gff = temp("species/{sample}/gene_annotation/{sample}.annotations.all.gff3"),
    genes_gff = temp("species/{sample}/gene_annotation/{sample}.annotations.genes.gff3"),
    genes_gtf = "species/{sample}/gene_annotation/{sample}.annotations.genes.gtf",
    coding_genes_gtf = "species/{sample}/gene_annotation/{sample}.annotations.coding_genes.gtf",
    genes_annot = "species/{sample}/gene_annotation/{sample}.gene_annotation.bed",
    genes_transcripts = "species/{sample}/gene_annotation/{sample}.gene_annotation.gene_transcript_id.txt",
    gene_type = "species/{sample}/gene_annotation/{sample}.gene_annotation.gene_type.txt",
    coding_genes_bed = "species/{sample}/gene_annotation/{sample}.gene_annotation.coding_gene.bed"
  shell:
    """
    gunzip -k {input.gff}
    echo -e "##gff-version 3" > {output.genes_gff}
    sed '/^##/d' {output.all_gff} | awk '{{if ($2 == "WormBase" || $2 == "WormBase_imported" || $2 == "Genbank" && $3 != "assembly_component" && $3 != "region") print}}' >> {output.genes_gff}
    gffread {output.genes_gff} -g {input.genome} -T -o {output.genes_gtf}
    gffread {output.genes_gff} -g {input.genome} -T -C -o {output.coding_genes_gtf}
    gtfToGenePred {output.genes_gtf} {output.genes_annot}.genePred
    genePredToBed {output.genes_annot}.genePred {output.genes_annot}.temp
    cat {output.genes_annot}.temp | awk 'BEGIN{{FS=":|\t";OFS="\t";}}{{print $1, $2, $3, $5, $6, $7, $8, $9, $10, $11, $12, $13}}' > {output.genes_annot}
    rm {output.genes_annot}.genePred
    rm {output.genes_annot}.temp
    python scripts/gene_transcript_id_extractor.py {output.genes_gff} {output.genes_transcripts}
    python scripts/gene_biotype_def.py {output.genes_gff} {output.genes_annot} {output.genes_transcripts} {output.gene_type} {output.coding_genes_bed}
    """


rule longest_isoform:
  input:
    coding_genes = "species/{sample}/gene_annotation/{sample}.gene_annotation.coding_gene.bed",
    genes_transcripts = "species/{sample}/gene_annotation/{sample}.gene_annotation.gene_transcript_id.txt",
    genome = "species/{sample}/genome/{sample}.fa"
  output:
    longest_isoform_bed = "species/{sample}/gene_annotation/{sample}.gene_annotation.coding_gene_longest_isoform.bed",
    longest_isoform_fa = "species/{sample}/gene_annotation/{sample}.gene_annotation.coding_gene_longest_isoform.fa",
    longest_isoform = "species/{sample}/gene_annotation/{sample}.gene_annotation.coding_gene_longest_isoform.longest_isoform",
  params:
    "species/{sample}/gene_annotation/{sample}.gene_annotation.coding_gene_longest_isoform"
  shell:
    """
    python scripts/gene_longest_isoform_annotation.py {input.coding_genes} {input.genes_transcripts} {params} {wildcards.sample}
    fastaFromBed -fi {input.genome} -bed {output.longest_isoform_bed} -fo {output.longest_isoform_fa} -name -split -s
    awk 'BEGIN{{FS="(";}}{{print $1}}' {output.longest_isoform_fa} > {output.longest_isoform_fa}.temp 
    mv {output.longest_isoform_fa}.temp {output.longest_isoform_fa}
    """


rule protein_annotation_download_all:
  output:
    elegans='species/elegans/protein_sequences/elegans.protein.to_parse.fa',
    inopinata='species/inopinata/protein_sequences/inopinata.protein.to_parse.fa',
    briggsae='species/briggsae/protein_sequences/briggsae.protein.to_parse.fa',
    remanei='species/remanei/protein_sequences/remanei.protein.to_parse.fa',
    nigoni='species/nigoni/protein_sequences/nigoni.protein.to_parse.fa',
    contortus='species/contortus/protein_sequences/contortus.protein.to_parse.fa',
    tipulae='species/tipulae/protein_sequences/tipulae.protein.to_parse.fa',
    bacteriophora='species/bacteriophora/protein_sequences/bacteriophora.protein.to_parse.fa',
    ceylanicum='species/ceylanicum/protein_sequences/ceylanicum.protein.to_parse.fa',
    pacificus='species/pacificus/protein_sequences/pacificus.protein.to_parse.fa',
    redivivus='species/redivivus/protein_sequences/redivivus.protein.to_parse.fa'
    becei='species/becei/protein_sequences/becei.protein.longest_isoform.fa',
    bovis='species/bovis/protein_sequences/bovis.protein.longest_isoform.fa',
    monodelphis='species/monodelphis/protein_sequences/monodelphis.protein.longest_isoform.fa',
    plicata='species/plicata/protein_sequences/plicata.protein.longest_isoform.fa',
    quiockensis='species/quiockensis/protein_sequences/quiockensis.protein.longest_isoform.fa',
    uteleia='species/uteleia/protein_sequences/uteleia.protein.longest_isoform.fa',
    exspectatus='species/exspectatus/protein_sequences/exspectatus.protein.longest_isoform.fa',
  shell:
    '''
    wget ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/sequence/protein/c_elegans.PRJNA13758.WS275.protein.fa.gz -O {output.elegans}.gz; gunzip -c {output.elegans}.gz > {output.elegans}; rm {output.elegans}.gz
    wget ftp://ftp.wormbase.org/pub/wormbase/species/c_inopinata/sequence/protein/c_inopinata.PRJDB5687.WS275.protein.fa.gz -O {output.inopinata}.gz; gunzip -c {output.inopinata}.gz > {output.inopinata}; rm {output.inopinata}.gz
    wget ftp://ftp.wormbase.org/pub/wormbase/species/c_briggsae/sequence/protein/c_briggsae.PRJNA10731.WS275.protein.fa.gz -O {output.briggsae}.gz; gunzip -c {output.briggsae}.gz > {output.briggsae}; rm {output.briggsae}.gz
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/010/183/535/GCA_010183535.1_CRPX506/GCA_010183535.1_CRPX506_protein.faa.gz -O {output.remanei}.gz; gunzip -c {output.remanei}.gz > {output.remanei}; rm {output.remanei}.gz
    wget ftp://ftp.wormbase.org/pub/wormbase/species/c_nigoni/sequence/protein/c_nigoni.PRJNA384657.WS275.protein.fa.gz -O {output.nigoni}.gz; gunzip -c {output.nigoni}.gz > {output.nigoni}; rm {output.nigoni}.gz
    wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species/haemonchus_contortus/PRJEB506/haemonchus_contortus.PRJEB506.WBPS14.protein.fa.gz -O {output.contortus}.gz; gunzip -c {output.contortus}.gz > {output.contortus}; rm {output.contortus}.gz
    wget ftp://ftp.wormbase.org/pub/wormbase/species/o_tipulae/sequence/protein/o_tipulae.PRJEB15512.WS275.protein.fa.gz -O {output.tipulae}.gz; gunzip -c {output.tipulae}.gz > {output.tipulae}; rm {output.tipulae}.gz
    wget ftp://ftp.wormbase.org/pub/wormbase/species/h_bacteriophora/sequence/protein/h_bacteriophora.PRJNA13977.WS248.protein.fa.gz -O {output.bacteriophora}.gz; gunzip -c {output.bacteriophora}.gz > {output.bacteriophora}; rm {output.bacteriophora}.gz
    wget ftp://ftp.wormbase.org/pub/wormbase/species/a_ceylanicum/PRJNA231479/sequence/protein/a_ceylanicum.PRJNA231479.WS248.protein.fa.gz -O {output.ceylanicum}.gz; gunzip -c {output.ceylanicum}.gz > {output.ceylanicum}; rm {output.ceylanicum}.gz
    wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS14/species/pristionchus_pacificus/PRJNA12644/pristionchus_pacificus.PRJNA12644.WBPS14.protein.fa.gz -O {output.pacificus}.gz; gunzip -c {output.pacificus}.gz > {output.pacificus}; rm {output.pacificus}.gz
    wget ftp://ftp.wormbase.org/pub/wormbase/species/p_redivivus/sequence/protein/p_redivivus.PRJNA186477.WS275.protein.fa.gz -O {output.redivivus}.gz; gunzip -c {output.redivivus}.gz > {output.redivivus}; rm {output.redivivus}.gz
    wget http://download.caenorhabditis.org/v1/sequence/Caenorhabditis_sp29_QG2083_v1.proteins.fa.gz -O {output.becei}.to_clean.gz; gunzip {output.becei}.to_clean.gz; awk 'BEGIN{{FS=" ";}}{{if ($1 ~ ">") print $1; else print $0}}' {output.becei}.to_clean | fold -s -w60 > {output.becei}.clean; rm {output.becei}.to_clean; python scripts/other_species_protein_fix.py {output.becei}.clean {output.becei}; rm {output.becei}.clean
    wget http://download.caenorhabditis.org/v1/sequence/Caenorhabditis_bovis_LS_v1.proteins.fa.gz -O {output.bovis}.to_clean.gz; gunzip {output.bovis}.to_clean.gz; awk 'BEGIN{{FS=" ";}}{{if ($1 ~ ">") print $1; else print $0}}' {output.bovis}.to_clean | fold -s -w60 > {output.bovis}.clean; rm {output.bovis}.to_clean; python scripts/other_species_protein_fix.py {output.bovis}.clean {output.bovis}; rm {output.bovis}.clean
    wget http://download.caenorhabditis.org/v1/sequence/Caenorhabditis_monodelphis_JU1667_v1.proteins.fa.gz -O {output.monodelphis}.to_clean.gz; gunzip {output.monodelphis}.to_clean.gz; awk 'BEGIN{{FS=" ";}}{{if ($1 ~ ">") print $1; else print $0}}' {output.monodelphis}.to_clean | fold -s -w60 > {output.monodelphis}.clean; rm {output.monodelphis}.to_clean; python scripts/other_species_protein_fix.py {output.monodelphis}.clean {output.monodelphis}; rm {output.monodelphis}.clean
    wget http://download.caenorhabditis.org/v1/sequence/Caenorhabditis_plicata_SB355_v1.proteins.fa.gz -O {output.plicata}.to_clean.gz; gunzip {output.plicata}.to_clean.gz; awk 'BEGIN{{FS=" ";}}{{if ($1 ~ ">") print $1; else print $0}}' {output.plicata}.to_clean | fold -s -w60 > {output.plicata}.clean; rm {output.plicata}.to_clean; python scripts/other_species_protein_fix.py {output.plicata}.clean {output.plicata}; rm {output.plicata}.clean
    wget http://download.caenorhabditis.org/v1/sequence/Caenorhabditis_sp38_JU2809_v1.proteins.fa.gz -O {output.quiockensis}.to_clean.gz; gunzip {output.quiockensis}.to_clean.gz; awk 'BEGIN{{FS=" ";}}{{if ($1 ~ ">") print $1; else print $0}}' {output.quiockensis}.to_clean | fold -s -w60 > {output.quiockensis}.clean; rm {output.quiockensis}.to_clean; python scripts/other_species_protein_fix.py {output.quiockensis}.clean {output.quiockensis}; rm {output.quiockensis}.clean
    wget http://download.caenorhabditis.org/v1/sequence/Caenorhabditis_sp31_JU2585_v1.proteins.fa.gz -O {output.uteleia}.to_clean.gz; gunzip {output.uteleia}.to_clean.gz; awk 'BEGIN{{FS=" ";}}{{if ($1 ~ ">") print $1; else print $0}}' {output.uteleia}.to_clean | fold -s -w60 > {output.uteleia}.clean; rm {output.uteleia}.to_clean; python scripts/other_species_protein_fix.py {output.uteleia}.clean {output.uteleia}; rm {output.uteleia}.clean
    wget http://pristionchus.org/variome/SNAP2012_exspectatus_proteins.fa.gz -O {output.exspectatus}.to_clean.gz; gunzip {output.exspectatus}.to_clean.gz; awk 'BEGIN{{FS=" ";}}{{if ($1 ~ ">") print $1; else print $0}}' {output.exspectatus}.to_clean | fold -s -w60 > {output.exspectatus}.clean; rm {output.exspectatus}.to_clean; python scripts/other_species_protein_fix.py {output.exspectatus}.clean {output.exspectatus}; rm {output.exspectatus}.clean
    '''


rule protein_annotation_processing:
  input:
    fa_to_parse = 'species/{sample}/protein_sequences/{sample}.protein.to_parse.fa',
    longest_isoform_annot = 'species/{sample}/gene_annotation/{sample}.gene_annotation.coding_gene_longest_isoform.longest_isoform',
    all_gff = 'species/{sample}/gene_annotation/{sample}.annotations.all.gff3.gz'
  output:
    fa_all = 'species/{sample}/protein_sequences/{sample}.protein.fa',
    longest_isoform = 'species/{sample}/protein_sequences/{sample}.protein.longest_isoform.fa',
    all_gff = temp('species/{sample}/gene_annotation/{sample}.annotations.all.gff3'),
  shell:
    '''
    awk '{{if ($1 ~ ">" && $3 ~ "=") print $1 ":" substr($3, 6); else if ($1 ~ ">" && $3 !~ "=" && $1 ~ $3) print $1 ":" $3; else if ($1 ~ ">" && $3 !~ "=" && $1 !~ $3) print $1 ":" substr($1,2); else print $0}}' {input.fa_to_parse} > {output.fa_all}
    gunzip -k {input.all_gff}
    python scripts/longest_protein_filtering.py {output.fa_all} {input.longest_isoform_annot} {output.longest_isoform} {output.all_gff} {wildcards.sample}
    '''


rule dfam_repeats:
  input:
    genome = "species/{sample}/genome/{sample}.fa",
    hmm_file = "software/dfamscan/caenorhabditis_elegans_dfam.hmm"
  output:
    dfam_out = "species/{sample}/repeats_dfam/{sample}_dfam.out",
    dfam_rep_id = "species/{sample}/repeats_dfam/{sample}_dfam.repeats.id.bed"
  resources:
    cpus=6
  shell:
    '''
    dfamscan.pl --cut_ga --species "Caenorhabditis elegans" --fastafile {input.genome} --hmmfile {input.hmm_file} --dfam_outfile {output.dfam_out} --cpu {resources.cpus}
    sed '1,5d' {output.dfam_out} | awk 'BEGIN{{OFS="\t";}}{{if ($9 == "+") print $3, $10, $11, $1, $2, $9; else print $3, $11, $10, $1, $2, "-"}}' > {output.dfam_rep_id}
    '''

