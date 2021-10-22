df_annot_ce_data = pd.DataFrame()

for k, v in config['annot_ce']['atac_samples'].items():
    sample_rep1 = 'annot_ce_atac_%(k)s_rep1' % locals()
    sample_rep2 = 'annot_ce_atac_%(k)s_rep2' % locals()
    l_fp_rep1 = v['rep1']
    l_fp_rep2 = v['rep2']

    df_annot_ce_data = df_annot_ce_data.append(pd.Series({
        'sample': sample_rep1,
        'condition': k,
        'assay': 'atac',
        'r1_raw': l_fp_rep1,
        'r1_pooled': 'samples/%(sample_rep1)s.r1.fq.gz' % locals(),
    }), ignore_index=True)
    df_annot_ce_data = df_annot_ce_data.append(pd.Series({
        'sample': sample_rep2,
        'condition': k,
        'assay': 'atac',
        'r1_raw': l_fp_rep2,
        'r1_pooled': 'samples/%(sample_rep2)s.r1.fq.gz' % locals(),
    }), ignore_index=True)


for k, v in config['annot_ce']['lcap_samples'].items():
    sample_rep1 = 'annot_ce_lcap_%(k)s_rep1' % locals()
    sample_rep2 = 'annot_ce_lcap_%(k)s_rep2' % locals()
    l_fp_rep1_read1 = v['rep1_read1']
    l_fp_rep1_read2 = v['rep1_read2']
    l_fp_rep2_read1 = v['rep2_read1']
    l_fp_rep2_read2 = v['rep2_read2']

    df_annot_ce_data = df_annot_ce_data.append(pd.Series({
        'sample': sample_rep1,
        'condition': k,
        'assay': 'lcap',
        'r1_raw': l_fp_rep1_read1,
        'r2_raw': l_fp_rep1_read2,
        'r1_pooled': 'samples/%(sample_rep1)s.r1.fq.gz' % locals(),
        'r2_pooled': 'samples/%(sample_rep1)s.r2.fq.gz' % locals(),
    }), ignore_index=True)
    df_annot_ce_data = df_annot_ce_data.append(pd.Series({
        'sample': sample_rep2,
        'condition': k,
        'assay': 'lcap',
        'r1_raw': l_fp_rep2_read1,
        'r2_raw': l_fp_rep2_read2,
        'r1_pooled': 'samples/%(sample_rep2)s.r1.fq.gz' % locals(),
        'r2_pooled': 'samples/%(sample_rep2)s.r2.fq.gz' % locals(),
    }), ignore_index=True)

df_annot_ce_data = df_annot_ce_data[['sample', 'assay', 'condition', 'r1_pooled', 'r2_pooled', 'r1_raw', 'r2_raw']]
df_annot_ce_data = df_annot_ce_data.set_index('sample')

def samples_ce_r1_input(wildcards):
    return df_annot_ce_data.query('index == "%s"' % (wildcards.sample,))['r1_raw'][0]

def samples_ce_r2_input(wildcards):
    return df_annot_ce_data.query('index == "%s"' % (wildcards.sample,))['r2_raw'][0]

rule samples_ce_r1:
    input:
        samples_ce_r1_input,
    output:
        'samples/{sample}.r1.fq.gz'
    run:
        arg_input = ' '.join(input)
        shell('cat %(arg_input)s > %(output)s' % locals())

rule samples_ce_r2:
    input:
        samples_ce_r2_input,
    output:
        'samples/{sample}.r2.fq.gz'
    run:
        arg_input = ' '.join(input)
        shell('cat %(arg_input)s > %(output)s' % locals())

rule rm_contigs_ce11:
    input:
        pf('{bid}', '{step}', '.bam', '{prefix}'),
        pf('{bid}', '{step}', '.bam.bai', '{prefix}')
    output:
        pf('{bid}', '{step}.rm_contigs_ce11', '.bam', '{prefix}'),
    shell:
        'samtools view -b {input[0]} I II III IV V X > {output}'

rule lcap_mean_by_stage:
    input:
        pf('{cond}_rep1', '{step}', '.bw', 'lcap'),
        pf('{cond}_rep2', '{step}', '.bw', 'lcap'),
    output:
        pf('{cond}', '{step}.mean_by_stage', '.bw', 'lcap'),
    shell: '''
        scripts/bigWiggleTools.ipy write_bg {output[0]} mean {input[0]} {input[1]}
    '''

rule annot_ce_atac:
    input:
        expand(pf('{sample}', config['annot_ce']['atac_step'], config['annot_ce']['atac_suffix'], 'atac'), sample=[*df_annot_ce_data.query('assay == "atac"').index]),
        expand(pf('{sample}', config['annot_ce']['atac_step_prp'], config['annot_ce']['atac_suffix'], 'atac'), sample=[*df_annot_ce_data.query('assay == "atac"').index]),
    output:
        'annot_ce/accessible_sites_ce.tsv',
        'annot_ce/notebooks/annot_ce_atac.html',
    shell: '''
        mkdir annot_ce/metrics_atac
        jupyter nbconvert --execute annot_ce/notebooks/annot_ce_atac.ipynb --ExecutePreprocessor.timeout=-1
    '''

rule annot_ce_canonical_geneset:
    input:
    output:
        'annot_ce/notebooks/annot_ce_canonical_geneset.html',
    shell: '''
        mkdir annot_ce/canonical_geneset
        jupyter nbconvert --execute annot_ce/notebooks/annot_ce_canonical_geneset.ipynb --ExecutePreprocessor.timeout=-1
    '''

rule annot_ce_exon:
    input:
        'annot_ce/accessible_sites_ce.tsv',
        'annot_ce/notebooks/annot_ce_canonical_geneset.html',
    output:
        'annot_ce/notebooks/annot_ce_exon.html',
    shell: '''
        mkdir annot_ce/metrics_exon
        jupyter nbconvert --execute annot_ce/notebooks/annot_ce_exon.ipynb --ExecutePreprocessor.timeout=-1
    '''

rule annot_ce_lcap:
    input:
        'annot_ce/accessible_sites_ce.tsv',
        expand(pf('{sample}', config['annot_ce']['lcap_firstbp_fwd'], '.bw', 'lcap'), sample=[*df_annot_ce_data.query('assay == "lcap"').index]),
        expand(pf('{sample}', config['annot_ce']['lcap_firstbp_rev'], '.bw', 'lcap'), sample=[*df_annot_ce_data.query('assay == "lcap"').index]),
    output:
        'annot_ce/notebooks/annot_ce_lcap.html',
    shell: '''
        mkdir annot_ce/metrics_lcap
        jupyter nbconvert --execute annot_ce/notebooks/annot_ce_lcap.ipynb --ExecutePreprocessor.timeout=-1
    '''

rule annot_ce_maxgap:
    input:
        'annot_ce/accessible_sites_ce.tsv',
        'annot_ce/notebooks/annot_ce_canonical_geneset.html',
        'annot_ce/notebooks/annot_ce_exon.html',
        expand(pf('annot_ce_lcap_{cond}', config['annot_ce']['lcap_filled_fwd'] + '.mean_by_stage', '.bw', 'lcap'), cond=[* config['annot_ce']['lcap_samples'].keys() ]),
        expand(pf('annot_ce_lcap_{cond}', config['annot_ce']['lcap_filled_rev'] + '.mean_by_stage', '.bw', 'lcap'), cond=[* config['annot_ce']['lcap_samples'].keys() ]),
    output:
        'annot_ce/notebooks/annot_ce_maxgap.html',
    shell: '''
	mkdir annot_ce/metrics_maxgap
        jupyter nbconvert --execute annot_ce/notebooks/annot_ce_maxgap.ipynb --ExecutePreprocessor.timeout=-1
    '''

rule annot_ce_type:
    input:
        'annot_ce/accessible_sites_ce.tsv',
        'annot_ce/notebooks/annot_ce_atac.html',
        'annot_ce/notebooks/annot_ce_canonical_geneset.html',
        'annot_ce/notebooks/annot_ce_exon.html',
        'annot_ce/notebooks/annot_ce_lcap.html',
        'annot_ce/notebooks/annot_ce_maxgap.html',
    output:
        'annot_ce/notebooks/annot_ce_type.html',
        'annot_ce/reg_elements_ce.bed',
    shell: '''
        mkdir annot_ce/metrics_type
        jupyter nbconvert --execute annot_ce/notebooks/annot_ce_type.ipynb --ExecutePreprocessor.timeout=-1
    '''

rule annot_ce:
    input:
        'annot_ce/reg_elements_ce.bed',

rule annot_ce_tests:
    input:
        expand('samples/{sample}.r1.fq.gz', sample=[*df_annot_ce_data.query('assay == "atac"').index]),
        expand('samples/{sample}.r1.fq.gz', sample=[*df_annot_ce_data.query('assay == "lcap"').index]),
        expand('samples/{sample}.r2.fq.gz', sample=[*df_annot_ce_data.query('assay == "lcap"').index]),
        #expand(pf('{sample}', config['annot_ce']['lcap_filled_fwd'], '.bw', 'lcap'), sample=[*df_annot_ce_data.query('assay == "lcap"').index]),
        #expand(pf('{sample}', config['annot_ce']['lcap_filled_rev'], '.bw', 'lcap'), sample=[*df_annot_ce_data.query('assay == "lcap"').index]),
