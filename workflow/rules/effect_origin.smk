

rule list_vcf_ids:
	'Obtain list of IID present in each chromosome.'
	input:
		'/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/vcf/{CHR}.vcf.gz'
	output:
		temp('results/effect_origin/aux/vcf_ids/temp/{CHR}-ids.txt')
	shell:
		'bcftools query -l {input[0]} > {output[0]}'

rule merge_vcf_ids:
	'Keep only IIDs present in all chromosomes.'
	input:
		expand('results/effect_origin/aux/vcf_ids/temp/{CHR}-ids.txt', CHR= CHROM)
	output:
		'results/effect_origin/aux/vcf_ids/allchr-ids.txt'
	run:
		df_list= list()
		for infile in input:
			d= pd.read_csv(infile, header= None, names= ['IID'])
			df_list.append(d)
		d= reduce(lambda x, y: pd.merge(x, y, on = 'IID', how = 'inner'), df_list)
		d.to_csv(output[0], sep= '\t', header= True, index= False)
	

rule list_trio_ids:
        'Make a list of trio IDs with genotype data.'
        input:
                '/mnt/work/pol/MoBaGenetics-1.0/delivery/sample-ids.txt',
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/aux/flaglist-merged/mobagen-flaglist-n99259.txt',
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/aux/pca/ethnic-core-samples',
                'results/effect_origin/aux/vcf_ids/allchr-ids.txt',
		'/mnt/work/pol/MoBaGenetics-1.0/delivery/trio-ids.txt'
        output:
                'results/effect_origin/aux/ids/fets_toextract.txt',
                'results/effect_origin/aux/ids/moms_toextract.txt',
                'results/effect_origin/aux/ids/dads_toextract.txt',
                'results/effect_origin/aux/ids/parent_offspring_trios.txt'
        run:
                d= pd.read_csv(input[0], sep= '\t', header= 0)
                flag= pd.read_csv(input[1], sep= '\t', header= 0)
                flag= flag.loc[flag.genotypesOK== True, :]
                flag= flag.loc[flag.phenoOK== True, :]
                pcs= [line.strip() for line in open(input[2], 'r')]
                x= pd.read_csv(input[3], sep= '\t', header= 0)
		fets= format_trios(d, 'Child', flag, pcs, x)
		print(fets.shape)
		moms= format_trios(d, 'Mother', flag, pcs, x)
		dads= format_trios(d, 'Father', flag, pcs, x)
                fets.to_csv(output[0], columns= ['IID'], sep= '\t', header= False, index= False)
                moms.to_csv(output[1], columns= ['IID'], sep= '\t', header= False, index= False)
                dads.to_csv(output[2], columns= ['IID'], sep= '\t', header= False, index= False)
		df= pd.read_csv(input[4], sep= '\t', header= 0)
		df= df.loc[df.Child.isin(pcs) & (df.Mother.isin(pcs)) & (df.Father.isin(pcs)), :]
		df= df.loc[(df.Child.isin(flag.IID.values)) & (df.Mother.isin(flag.IID.values)) & (df.Father.isin(flag.IID.values)), :]
		df= df.loc[df.BATCH != 'TED', ]
		df.drop_duplicates('PREG_ID_1724', inplace= True, keep= 'first')
		df.drop_duplicates('Mother', inplace= True, keep= 'first')
		df.drop_duplicates('Father', inplace= True, keep= 'first')
                df.to_csv(output[3], sep= '\t', header= True, index= False)

rule format_sumstats:
	'Remove non-necessary rows from summary statistics.'
	input:
		'results/GWAS/sumstats_non_rel_lm/GWAS-miscarriage/allchr-fets.txt.gz'
	output:
		temp('results/effect_origin/aux/top_signals/miscarriage/fets-regions_to_extract.txt'),
	run:
		df_list= list()
		for infile in input:
			d= pd.read_csv(infile, sep = ' ', header= 0)
			d= d.loc[d.LOG10P> 3, :]
			d= d.loc[(d.A1FREQ> 0.005 ) & (d.A1FREQ<0.995), :]
			d.sort_values('LOG10P', inplace= True, ascending= False, ignore_index=True)
			d['GENPOS']= d.GENPOS.apply(int).apply(str)
			d['pos2']= d.GENPOS
			d['CHROM']= d.CHROM.apply(str)
			d['CHROM']= np.where(d.CHROM== '23', 'X', d.CHROM)
			df_list.append(d)
		d= pd.concat(df_list)
		d.sort_values(['CHROM', 'GENPOS'], inplace= True, ascending= True)
		d.to_csv(output[0], header= False, index= False, sep= '\t', columns= ['CHROM', 'GENPOS', 'pos2'])

rule get_DS_effect_origin:
        'Extract DS from VCF file for a subset of genetic variants.'
        input:
                'results/effect_origin/aux/top_signals/miscarriage/fets-regions_to_extract.txt',
                'results/effect_origin/aux/ids/fets_toextract.txt',
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/vcf/{CHR}.vcf.gz'
        output:
                temp('results/effect_origin/aux/DS/temp/miscarriage/fets_ds{CHR}')
        run:
                shell("bcftools query -S {input[1]} -R {input[0]} -f '%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n' {input[2]} -o {output[0]}")

rule add_header_DS_effect_origin:
        'Add header to genotype files.'
        input:
                'results/effect_origin/aux/ids/fets_toextract.txt',
                'results/effect_origin/aux/DS/temp/miscarriage/fets_ds{CHR}'
        output:
                'results/effect_origin/DS/miscarriage/fets_DS{CHR}'
        run:
                cols= ['chr','pos','ref','eff'] + [line.strip() for line in open(input[0], 'r')]
                d= pd.read_csv(input[1], header= None, names= cols, sep= '\t')
                d.drop_duplicates(['chr', 'pos'], keep=False, inplace= True)
                d.to_csv(output[0], sep= '\t', header= True, index= False)

rule poisson_regression:
	''
	input:
		'results/effect_origin/DS/miscarriage/fets_DS{CHR}',
		'results/pheno_non_rel/fets_phenotype_miscarriage.txt',
		'results/pheno_non_rel/fets_covars_miscarriage.txt'
	output:
		'results/effect_origin/results/temp/miscarriage/fets_poisson_{CHR}.txt'
	conda: '../envs/plots.yml'
	script:
		'../scripts/effect_origin_conditional.R'

rule concat_poisson:
	'Collect DS from all CHR.'
	input:
		expand('results/effect_origin/results/temp/miscarriage/fets_poisson_{CHR}.txt', CHR= CHROM)
	output:
		'results/effect_origin/poisson/miscarriage/fets_poisson.txt'
	shell:
                '''
                set +o pipefail;
                head -1 {input[0]} > {output[0]}
                cat {input} | grep -v 'term' >> {output[0]}
                '''

