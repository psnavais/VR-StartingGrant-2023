

rule format_sumstats_LDAK:
	'Format proxy-SNP sumstats for LDAK.'
	input:
		'results/GWAS/sumstats_non_rel_lm/GWAS-{pheno}/allchr-{sample}.txt.gz',
	output:
		'results/LDAK/sumstats/temp/{pheno}/allchr-{sample}-LDAK.txt'
	run:
		d= pd.read_csv(input[0], sep= ' ', header= 0)
		d= d.loc[(d.A1FREQ< 0.995) & (d.A1FREQ> 0.005), :]
		d['Predictor']= d.CHROM.astype(str) + ':' + d.GENPOS.astype(str)
		d['Z']= d.BETA / d.SE
		if wildcards.sample== 'fets':
			d.N= d.N / 4
		d= d.loc[:, ['Predictor', 'ALLELE1', 'ALLELE0', 'N', 'Z']]
		d.drop_duplicates('Predictor', inplace= True)
		d.columns= ['Predictor', 'A1', 'A2', 'n', 'Z']
		d.to_csv(output[0], sep= '\t', header= True, index= False)

rule LDAK:
	''
	input:
		'results/LDAK/sumstats/temp/{pheno}/allchr-{sample}-LDAK.txt',
		'/home/pol.sole.navais/soft/ldak/bld.ldak.hapmap.gbr.tagging'
	output:
		'results/LDAK/sumstats/heritability/{pheno}/{sample}-LDAK.hers'
	params:
		'results/LDAK/sumstats/heritability/{pheno}/{sample}-LDAK'
	shell:
		'/home/pol.sole.navais/soft/ldak/ldak5.2.linux --tagfile {input[1]} --summary {input[0]} --check-sums NO --sum-hers {params[0]}'

rule check_results_LDAK:
        ''
        input:
                expand('results/LDAK/sumstats/heritability/{pheno}/{sample}-LDAK.hers', pheno= pheno_file['phenotypes'], sample= fam_ids['fam_id'])
        output:
                'results/LDAK/checks/h2_performed.txt'
        shell:
                'touch {output[0]}'

