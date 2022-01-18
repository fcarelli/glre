rule extra_species_THAP_analysis:
  output:
    placei='species/placei/protein_sequences/placei.protein.longest_isoform.fa',
    caninum='species/caninum/protein_sequences/caninum.protein.longest_isoform.fa',
    polygyrus='species/polygyrus/protein_sequences/polygyrus.protein.longest_isoform.fa',
    brasiliensis='species/brasiliensis/protein_sequences/brasiliensis.protein.longest_isoform.fa',
    viviparus='species/viviparus/protein_sequences/viviparus.protein.longest_isoform.fa',
    rhodensis='species/rhodensis/protein_sequences/rhodensis.protein.longest_isoform.fa',
    sinica='species/sinica/protein_sequences/sinica.protein.longest_isoform.fa',
    afra='species/afra/protein_sequences/afra.protein.longest_isoform.fa',
    panamensis='species/panamensis/protein_sequences/panamensis.protein.longest_isoform.fa',
    virilis='species/virilis/protein_sequences/virilis.protein.longest_isoform.fa',
  resources:
    download_streams=1
  shell:
    '''
    wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS15/species/haemonchus_placei/PRJEB509/haemonchus_placei.PRJEB509.WBPS15.protein.fa.gz -O {output.placei}.to_clean.gz; gunzip {output.placei}.to_clean.gz; awk 'BEGIN{{FS=" ";}}{{if ($1 ~ ">") print $1; else print $0}}' {output.placei}.to_clean | fold -s -w60 > {output.placei}.clean; rm {output.placei}.to_clean; python scripts/other_species_protein_fix.py {output.placei}.clean {output.placei}; rm {output.placei}.clean
    wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS15/species/ancylostoma_caninum/PRJNA72585/ancylostoma_caninum.PRJNA72585.WBPS15.protein.fa.gz -O {output.caninum}.to_clean.gz; gunzip {output.caninum}.to_clean.gz; awk 'BEGIN{{FS=" ";}}{{if ($1 ~ ">") print $1; else print $0}}' {output.caninum}.to_clean | fold -s -w60 > {output.caninum}.clean; rm {output.caninum}.to_clean; python scripts/other_species_protein_fix.py {output.caninum}.clean {output.caninum}; rm {output.caninum}.clean
    wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS15/species/heligmosomoides_polygyrus/PRJEB15396/heligmosomoides_polygyrus.PRJEB15396.WBPS15.protein.fa.gz -O {output.polygyrus}.to_clean.gz; gunzip {output.polygyrus}.to_clean.gz; awk 'BEGIN{{FS=" ";}}{{if ($1 ~ ">") print $1; else print $0}}' {output.polygyrus}.to_clean | fold -s -w60 > {output.polygyrus}.clean; rm {output.polygyrus}.to_clean; python scripts/other_species_protein_fix.py {output.polygyrus}.clean {output.polygyrus}; rm {output.polygyrus}.clean
    wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS15/species/nippostrongylus_brasiliensis/PRJEB511/nippostrongylus_brasiliensis.PRJEB511.WBPS15.protein.fa.gz -O {output.brasiliensis}.to_clean.gz; gunzip {output.brasiliensis}.to_clean.gz; awk 'BEGIN{{FS=" ";}}{{if ($1 ~ ">") print $1; else print $0}}' {output.brasiliensis}.to_clean | fold -s -w60 > {output.brasiliensis}.clean; rm {output.brasiliensis}.to_clean; python scripts/other_species_protein_fix.py {output.brasiliensis}.clean {output.brasiliensis}; rm {output.brasiliensis}.clean
    wget ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS15/species/dictyocaulus_viviparus/PRJNA72587/dictyocaulus_viviparus.PRJNA72587.WBPS15.protein.fa.gz -O {output.viviparus}.to_clean.gz; gunzip {output.viviparus}.to_clean.gz; awk 'BEGIN{{FS=" ";}}{{if ($1 ~ ">") print $1; else print $0}}' {output.viviparus}.to_clean | fold -s -w60 > {output.viviparus}.clean; rm {output.viviparus}.to_clean; python scripts/other_species_protein_fix.py {output.viviparus}.clean {output.viviparus}; rm {output.viviparus}.clean
    wget https://figshare.com/ndownloader/files/24649328?private_link=7787e173bba6f44d0aa5 -O {output.rhodensis}.to_clean; fold -s -w60 {output.rhodensis}.to_clean > {output.rhodensis}.clean; rm {output.rhodensis}.to_clean; python scripts/other_species_protein_fix.py {output.rhodensis}.clean {output.rhodensis}; rm {output.rhodensis}.clean
    wget http://download.caenorhabditis.org/v1/sequence/Caenorhabditis_sinica_JU800_v1.proteins.fa.gz -O {output.sinica}.to_clean.gz; gunzip {output.sinica}.to_clean.gz; awk 'BEGIN{{FS=" ";}}{{if ($1 ~ ">") print $1; else print $0}}' {output.sinica}.to_clean | fold -s -w60 > {output.sinica}.clean; rm {output.sinica}.to_clean; python scripts/other_species_protein_fix.py {output.sinica}.clean {output.sinica}; rm {output.sinica}.clean
    wget http://download.caenorhabditis.org/v1/sequence/Caenorhabditis_afra_JU1286_v1.proteins.fa.gz -O {output.afra}.to_clean.gz; gunzip {output.afra}.to_clean.gz; awk 'BEGIN{{FS=" ";}}{{if ($1 ~ ">") print $1; else print $0}}' {output.afra}.to_clean | fold -s -w60 > {output.afra}.clean; rm {output.afra}.to_clean; python scripts/other_species_protein_fix.py {output.afra}.clean {output.afra}; rm {output.afra}.clean
    wget http://download.caenorhabditis.org/v1/sequence/Caenorhabditis_sp28_QG2080_v1.proteins.fa.gz -O {output.panamensis}.to_clean.gz; gunzip {output.panamensis}.to_clean.gz; awk 'BEGIN{{FS=" ";}}{{if ($1 ~ ">") print $1; else print $0}}' {output.panamensis}.to_clean | fold -s -w60 > {output.panamensis}.clean; rm {output.panamensis}.to_clean; python scripts/other_species_protein_fix.py {output.panamensis}.clean {output.panamensis}; rm {output.panamensis}.clean
    wget http://download.caenorhabditis.org/v1/sequence/Caenorhabditis_virilis_JU1968_v1.proteins.fa.gz -O {output.virilis}.to_clean.gz; gunzip {output.virilis}.to_clean.gz; awk 'BEGIN{{FS=" ";}}{{if ($1 ~ ">") print $1; else print $0}}' {output.virilis}.to_clean | fold -s -w60 > {output.virilis}.clean; rm {output.virilis}.to_clean; python scripts/other_species_protein_fix.py {output.virilis}.clean {output.virilis}; rm {output.virilis}.clean
    '''
    

rule protein_blastdb:
  input:
    'species/{sample}/protein_sequences/{sample}.protein.longest_isoform.fa'
  output:
    'species/{sample}/protein_sequences/blastdb/{sample}.protein.longest_isoform.fa.phr',
    'species/{sample}/protein_sequences/blastdb/{sample}.protein.longest_isoform.fa.pin',
    'species/{sample}/protein_sequences/blastdb/{sample}.protein.longest_isoform.fa.psq',
  shell:
    '''
    makeblastdb -in {input} -dbtype prot -out species/{wildcards.sample}/protein_sequences/blastdb/{wildcards.sample}.protein.longest_isoform.fa
    '''


rule him17_elegans:
  input:
    'species/elegans/protein_sequences/elegans.protein.longest_isoform.fa',
  output:
    'TF_evolution/HIM17.prot.fa',
  shell:
    '''
    samtools faidx {input[0]} WBGene00001874 > {output[0]}
    '''


rule him17_orthologs:
  input:
    'TF_evolution/HIM17.prot.fa',
    'species/{sample}/protein_sequences/{sample}.protein.longest_isoform.fa',
    'species/{sample}/protein_sequences/blastdb/{sample}.protein.longest_isoform.fa.phr',
    'species/{sample}/protein_sequences/blastdb/{sample}.protein.longest_isoform.fa.pin',
    'species/{sample}/protein_sequences/blastdb/{sample}.protein.longest_isoform.fa.psq',
  output:
    'TF_evolution/orthologs/HIM17.{sample}.aln',
    'TF_evolution/orthologs/HIM17.{sample}.fa',
    temp('species/{sample}/protein_sequences/{sample}.protein.longest_isoform.fa.clean'),
    temp('species/{sample}/protein_sequences/{sample}.protein.longest_isoform.fa.clean.fai'),
  resources:
    cpus=8
  shell:
    '''
    short_species=$(echo {wildcards.sample} | cut -c1-3)
    blastp -query {input[0]} -db species/{wildcards.sample}/protein_sequences/blastdb/{wildcards.sample}.protein.longest_isoform.fa -outfmt 6 > {output[0]}
    touch {output[1]}
    if [ -s {output[0]} ]
    then
    evalue=$(head -n 1 {output[0]} | awk '{{print $11}}' | sed 's/e/*10^/g;s/ /*/')
    if (( $(echo "$evalue < 0.00001" |bc -l) ))
    then
    HIM17_ortho=$(head -n 1 {output[0]} | awk '{{print $2}}')
    echo -e ">"$short_species"_HIM17" > {output[1]}
    python scripts/fasta_seq_len_fix.py {input[1]} {output[2]}
    samtools faidx {output[2]} $HIM17_ortho | sed '1d' >> {output[1]}
    fi
    fi
    touch {output[2]}
    touch {output[3]}
    '''


rule him17_alignment:
  input:
    HIM17_all=expand('TF_evolution/orthologs/HIM17.{sample}.fa', sample=THAP_ANALYSIS_EURHABDITIS),
  output:
    'TF_evolution/orthologs/HIM17.all.fa',
    'TF_evolution/orthologs/HIM17.all.dist',
    'TF_evolution/orthologs/HIM17.all.aln',
  shell:
    '''
    cat {input.HIM17_all} >> {output[0]}
    clustalo --infile={output[0]} --distmat-out={output[1]}.temp --outfile={output[2]} --full-iter --full --percent-id --output-order=tree-order
    sed '1d' {output[1]}.temp > {output[1]}; rm {output[1]}.temp
    '''


rule elegans_other_species_ortholog_id:
  input:
    'data/external_data/elegans.other_caenorhabditis.orthologs.id.gz',
    'data/external_data/elegans.other_nematodes.orthologs.id.gz',
  output:
    temp('data/external_data/elegans.other_caenorhabditis.orthologs.id'),
    temp('data/external_data/elegans.other_nematodes.orthologs.id'),
    'data/external_data/elegans.inopinata.1to1.txt',
    'data/external_data/elegans.briggsae.1to1.txt',
    'data/external_data/elegans.nigoni.1to1.txt',
    'data/external_data/elegans.remanei.1to1.txt',
    'data/external_data/elegans.ceylanicum.1to1.txt',
    'data/external_data/elegans.contortus.1to1.txt',
    'data/external_data/elegans.bacteriophora.1to1.txt',
    'data/external_data/elegans.tipulae.1to1.txt',
  shell:
    '''
    gunzip -c {input[0]} > {output[0]}
    gunzip -c {input[1]} > {output[1]}
    python scripts/biomart_elegans_orthologs_parser.caenorhabditis.py {output[0]} {output[2]} {output[3]} {output[4]} {output[5]}
    python scripts/biomart_elegans_orthologs_parser.other_nematodes.py {output[1]} {output[6]} {output[7]} {output[8]} {output[9]}
    '''


rule briggsae_other_species_ortholog_id:
  input:
    'data/external_data/briggsae.other_caenorhabditis.orthologs.id.gz'
  output:
    temp('data/external_data/briggsae.other_caenorhabditis.orthologs.id'),
    'data/external_data/briggsae.nigoni.1to1.txt',
    'data/external_data/briggsae.remanei.1to1.txt',
    'data/external_data/briggsae.elegans.1to1.txt',
    'data/external_data/briggsae.inopinata.1to1.txt',
  shell:
    '''
    gunzip -c {input[0]} > {output[0]}
    python scripts/biomart_briggsae_orthologs_parser.caenorhabditis.py {output[0]} {output[1]} {output[2]} {output[3]} {output[4]}
    '''


rule THAP_profile_download:
  output:
    'data/external_data/THAP.hmm'
  shell:
    '''
    wget http://pfam.xfam.org/family/PF05485/hmm -O {output[0]}
    '''


rule THAP_alignment:
  input:
    'data/external_data/THAP.hmm',
    'species/{sample}/protein_sequences/{sample}.protein.longest_isoform.fa',
    'TF_evolution/orthologs/HIM17.{sample}.aln',
  output:
    'TF_evolution/THAP_conservation/{sample}.THAP.out',
    'TF_evolution/THAP_conservation/{sample}.THAP.parsed',
    temp('species/{sample}/protein_sequences/{sample}.protein.longest_isoform.fa.clean'),
    temp('species/{sample}/protein_sequences/{sample}.protein.longest_isoform.fa.clean.fai'),
  shell:
    '''
    hmmsearch --pfamtblout {output[0]} {input[0]} {input[1]}
    python scripts/hmmer_THAP_parser.py {output[0]} {output[1]}.temp {input[2]}
    python scripts/fasta_seq_len_fix.py {input[1]} {output[2]}
    if [ -s {output[1]}.temp ]
    then
    ortholog=$(head -n 1 {output[1]}.temp | cut -f 1)
    python scripts/fasta_seq_len_fix.py {input[1]} {output[2]}
    ortholog_len=$(samtools faidx {output[2]} $ortholog | bioawk -c fastx '{{ print length($seq) }}')
    awk -v var="$ortholog_len" '{{print $0 "\t" var}}' {output[1]}.temp > {output[1]}
    rm {output[1]}.temp
    else
    rm {output[1]}.temp
    touch {output[1]}
    touch {output[3]}
    fi
    '''


rule orthologs_plots_R:
  input:
    'TF_evolution/orthologs/HIM17.all.dist',
    'data/external_data/elegans.inopinata.1to1.txt',
    'data/external_data/elegans.briggsae.1to1.txt',
    'data/external_data/elegans.nigoni.1to1.txt',
    'data/external_data/elegans.remanei.1to1.txt',
    'data/external_data/briggsae.nigoni.1to1.txt',
    'data/external_data/briggsae.remanei.1to1.txt',
    'data/external_data/briggsae.elegans.1to1.txt',
    'data/external_data/briggsae.inopinata.1to1.txt',
    expand('TF_evolution/THAP_conservation/{sample}.THAP.parsed', sample = THAP_ANALYSIS_SPECIES),
  output:
    'plots/HIM17_identity.heat.pdf',
    'plots/HIM17_orthologs_identity.pdf',
    'plots/THAP_domains_HIM17.pdf'
  shell:
    '''
    Rscript scripts/him17_orthologs.R
    '''

rule positive_selection_test:
  input:
    'TF_evolution/positive_selection/HIM17_prot.all.fa',
    'TF_evolution/positive_selection/HIM17_cdna.all.fa',
    'TF_evolution/positive_selection/HIM17.all.tree',
    'TF_evolution/positive_selection/branch_site_null/codeml_him17_null.ctl',
    'TF_evolution/positive_selection/branch_site_alt/codeml_him17_alt.ctl',
  output:
    'TF_evolution/positive_selection/HIM17_prot.all.mafft.fa',
    'TF_evolution/positive_selection/HIM17.all.linsi.pal2nal.codon.phy',
    'TF_evolution/positive_selection/branch_site_null/HIM17.all.modelA.null',
    'TF_evolution/positive_selection/branch_site_alt/HIM17.all.modelA.alt',
  shell:
    '''
    mafft --localpair --maxiterate 1000 {input[0]} > {output[0]}
    pal2nal.pl {outpu[0]} {input[1]} -output paml -nogap > {output[1]}
    cd TF_evolution/positive_selection/branch_site_null; codeml codeml_him17_null.ctl
    cd ../branch_site_alt; codeml codeml_him17_alt.ctl
    '''

