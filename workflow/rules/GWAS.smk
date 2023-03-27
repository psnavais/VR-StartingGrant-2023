

rule filter_mac:
        'Filter mac for Step 1 of REGENIE'
        input:
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/subset/grm-high-quality-pruned/grm-high-quality-pruned.bed',
                'results/aux/ids/samples/{sample}_ids_{pheno}.txt'
        output:
                'results/GWAS/regenie/step1/snp_to_filter/{sample}_filtered_{pheno}.snplist'
        params:
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/subset/grm-high-quality-pruned/grm-high-quality-pruned',
                'results/GWAS/regenie/step1/snp_to_filter/{sample}_filtered_{pheno}'
        shell:
                '''
                /home/pol.sole.navais/soft/plink2 \
                --bfile {params[0]} \
                --mac 100 \
                --write-snplist \
                --keep {input[1]} \
                --out {params[1]}
                '''

rule list_bgen_samples:
        'List of bGEN samples.'
        input:
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/bgen/{CHR}.bgen'
        output:
                'results/aux/ids/samples/bgen/{CHR}_samples.txt'
        shell:
                '/home/pol.sole.navais/soft/qctool_v2.2.0/qctool -g {input[0]} -os {output[0]}'

rule chrX_to_diploid:
	'Use PLINK2 to convert BGEN haploid males to diplod (pgen file format).'
	input:
		'/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/vcf/X.vcf.gz'
	output:
		temp(multiext('results/GWAS/pgen/X', '.pgen', '.pvar', '.psam'))
	params:
		'results/GWAS/pgen/X'
	shell:
		"/home/pol.sole.navais/soft/plink2 --vcf {input[0]} dosage='DS' --double-id --make-pgen psam-cols=+fid --out {params[0]}"


rule REGENIE_step1:
        'Whole genome regression model is fit to the traits.'
        input:
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/subset/grm-high-quality-pruned/grm-high-quality-pruned.bed',
                'results/pheno/{sample}_phenotype_{pheno}.txt',
                'results/pheno/{sample}_covars_{pheno}.txt',
                'results/aux/ids/samples/{sample}_ids_{pheno}.txt',
                'results/GWAS/regenie/step1/snp_to_filter/{sample}_filtered_{pheno}.snplist'
        output:
                temp('results/GWAS/regenie/step1/results/{sample}/{pheno}_1.loco.gz'),
                temp('results/GWAS/regenie/step1/results/{sample}/{pheno}_pred.list')
        params:
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/subset/grm-high-quality-pruned/grm-high-quality-pruned',
                'results/GWAS/regenie/step1/results/{sample}/{pheno}',
                'results/GWAS/regenie/step1/results/{sample}/{pheno}_temp'
	threads: 30
	run:
		if wildcards.pheno == 'miscarriage_bin':
                	shell("""
		/home/pol.sole.navais/soft/regenie_v3.2.1.gz_x86_64_Linux \
                --step 1 \
                --threads {threads} \
                --gz \
                --bed {params[0]} \
                --covarFile {input[2]} \
                --phenoFile {input[1]} \
                --keep {input[3]} \
                --extract {input[4]} \
                --bsize 1000 \
                --bt --lowmem \
                --lowmem-prefix {params[2]} \
                --catCovarList cohort \
                --out {params[1]}
                """)
		else:
			shell("""
                /home/pol.sole.navais/soft/regenie_v3.2.1.gz_x86_64_Linux \
                --step 1 \
                --threads {threads} \
                --gz \
                --bed {params[0]} \
                --covarFile {input[2]} \
                --phenoFile {input[1]} \
                --keep {input[3]} \
                --extract {input[4]} \
                --bsize 1000 \
		--force-qt \
                --lowmem \
                --lowmem-prefix {params[2]} \
                --catCovarList cohort \
                --out {params[1]}
                """)

rule REGENIE_step2:
	'Whole genome regression model is fit to the traits.'
	input:
		'/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/bgen/{CHR}.bgen',
		'results/pheno/{sample}_phenotype_{pheno}.txt',
		'results/pheno/{sample}_covars_{pheno}.txt',
		'results/aux/ids/samples/{sample}_ids_{pheno}.txt',
		'results/GWAS/regenie/step1/results/{sample}/{pheno}_pred.list',
		'results/aux/ids/samples/bgen/{CHR}_samples.txt',
		'results/GWAS/regenie/step1/results/{sample}/{pheno}_1.loco.gz',
		multiext('results/GWAS/pgen/X', '.pgen', '.pvar', '.psam')
	output:
		temp('results/GWAS/regenie/step2/temp/{sample}/{CHR}_results_{pheno}.regenie')
	params:
		'results/GWAS/regenie/step2/temp/{sample}/{CHR}_results',
		'results/GWAS/pgen/X'
	threads: 5
	run:
		if wildcards.CHR != 'X':
			if wildcards.pheno == 'miscarriage_bin':
				shell('''
                /home/pol.sole.navais/soft/regenie_v3.2.1.gz_x86_64_Linux \
                --step 2 \
                --bgen {input[0]} \
                --covarFile {input[2]} \
                --keep {input[3]} \
                --phenoFile {input[1]} \
                --bsize 400 \
                --bt \
                --firth --approx \
                --minINFO 0.6 \
                --threads {threads} \
                --sample {input[5]} \
                --pred {input[4]} \
                --out {params[0]} \
		--af-cc \
                --catCovarList cohort \
                --verbose
		                ''')
			else:
				shell('''
                /home/pol.sole.navais/soft/regenie_v3.2.1.gz_x86_64_Linux \
                --step 2 \
                --bgen {input[0]} \
                --covarFile {input[2]} \
                --keep {input[3]} \
                --phenoFile {input[1]} \
                --bsize 400 \
                --firth --approx \
                --minINFO 0.6 \
                --threads {threads} \
                --sample {input[5]} \
                --pred {input[4]} \
		--force-qt \
                --out {params[0]} \
                --catCovarList cohort \
                --verbose
                                ''')
		else:
			if wildcards.pheno == 'miscarriage_bin':
                                shell('''
                /home/pol.sole.navais/soft/regenie_v3.2.1.gz_x86_64_Linux \
                --step 2 \
                --pgen {params[1]} \
                --covarFile {input[2]} \
                --keep {input[3]} \
                --phenoFile {input[1]} \
                --bsize 400 \
                --firth --approx \
		--bt \
		--af-cc \
                --minINFO 0.6 \
                --threads {threads} \
                --sample {input[5]} \
                --pred {input[4]} \
                --out {params[0]} \
                --catCovarList cohort \
                --verbose
                                ''')
			else:
				shell('''
                /home/pol.sole.navais/soft/regenie_v3.2.1.gz_x86_64_Linux \
                --step 2 \
                --pgen {params[1]} \
                --covarFile {input[2]} \
                --keep {input[3]} \
                --phenoFile {input[1]} \
                --bsize 400 \
                --firth --approx \
                --minINFO 0.6 \
                --threads {threads} \
		--force-qt \
                --pred {input[4]} \
                --out {params[0]} \
                --catCovarList cohort \
                --verbose
                		''')

rule concat_GWAS_results:
	'Concatenate results from regenie'
	input:
		expand('results/GWAS/regenie/step2/temp/{{sample}}/{CHR}_results_{{pheno}}.regenie', CHR= CHROM)
	output:
		temp('results/GWAS/sumstats/GWAS-{pheno}/temp/allchr/{sample}.txt')
	shell:
		'''
		head -1 {input[0]} > {output[0]}
		tail -n +2 -q {input} >> {output[0]}
		'''

rule gzip_results:
        'Gzip results.'
        input:
                'results/GWAS/sumstats/GWAS-{pheno}/temp/allchr/{sample}.txt'
        output:
                'results/GWAS/sumstats/GWAS-{pheno}/allchr-{sample}.txt.gz'
        shell:
                'gzip -c {input[0]} > {output[0]}'

rule check_results:
	''
	input:
		expand('results/GWAS/sumstats/GWAS-{pheno}/allchr-{sample}.txt.gz', pheno= pheno_file['phenotypes'], sample= fam_ids['fam_id'])
	output:
		'results/GWAS/checks/GWAS_performed.txt'
	shell:
		'touch {output[0]}'


### Non relatives

rule REGENIE_step1_non_rel:
        'Whole genome regression model is fit to the traits.'
        input:
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/subset/grm-high-quality-pruned/grm-high-quality-pruned.bed',
                'results/pheno_non_rel/{sample}_phenotype_{pheno}.txt',
                'results/pheno_non_rel/{sample}_covars_{pheno}.txt',
                'results/aux/ids_non_rel/samples/{sample}_ids_{pheno}.txt',
                'results/GWAS/regenie/step1/snp_to_filter/{sample}_filtered_{pheno}.snplist'
        output:
                temp('results/GWAS/regenie_non_rel/step1/results/{sample}/{pheno}_1.loco.gz'),
                temp('results/GWAS/regenie_non_rel/step1/results/{sample}/{pheno}_pred.list')
        params:
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/subset/grm-high-quality-pruned/grm-high-quality-pruned',
                'results/GWAS/regenie_non_rel/step1/results/{sample}/{pheno}',
                'results/GWAS/regenie_non_rel/step1/results/{sample}/{pheno}_temp'
	threads: 30
	run:
		shell("""
                /home/pol.sole.navais/soft/regenie_v3.2.1.gz_x86_64_Linux \
                --step 1 \
                --threads {threads} \
                --gz \
                --bed {params[0]} \
                --covarFile {input[2]} \
                --phenoFile {input[1]} \
                --keep {input[3]} \
                --extract {input[4]} \
                --bsize 1000 \
                --force-qt \
                --lowmem \
                --lowmem-prefix {params[2]} \
                --catCovarList cohort \
                --out {params[1]}
                """)



rule REGENIE_step2_non_rel:
        'Whole genome regression model is fit to the traits.'
        input:
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/bgen/{CHR}.bgen',
                'results/pheno_non_rel/{sample}_phenotype_{pheno}.txt',
                'results/pheno_non_rel/{sample}_covars_{pheno}.txt',
                'results/aux/ids_non_rel/samples/{sample}_ids_{pheno}.txt',
                'results/GWAS/regenie_non_rel/step1/results/{sample}/{pheno}_pred.list',
                'results/aux/ids/samples/bgen/{CHR}_samples.txt',
                'results/GWAS/regenie_non_rel/step1/results/{sample}/{pheno}_1.loco.gz'
        output:
                temp('results/GWAS/regenie_non_rel/step2/temp/{sample}/{CHR}_results_{pheno}.regenie')
        params:
                'results/GWAS/regenie_non_rel/step2/temp/{sample}/{CHR}_results',
                'results/GWAS/pgen/X'
	threads: 5
	run:
		shell('''
                /home/pol.sole.navais/soft/regenie_v3.2.1.gz_x86_64_Linux \
                --step 2 \
                --bgen {input[0]} \
                --covarFile {input[2]} \
                --keep {input[3]} \
                --phenoFile {input[1]} \
                --bsize 400 \
                --firth --approx \
                --minINFO 0.6 \
                --threads {threads} \
                --sample {input[5]} \
                --pred {input[4]} \
                --force-qt \
                --out {params[0]} \
                --catCovarList cohort \
                --verbose
                                ''')

rule concat_GWAS_results_non_rel:
        'Concatenate results from regenie'
        input:
                expand('results/GWAS/regenie_non_rel/step2/temp/{{sample}}/{CHR}_results_{{pheno}}.regenie', CHR= CHROM)
        output:
                temp('results/GWAS/sumstats_non_rel/GWAS-{pheno}/temp/allchr/{sample}.txt')
        shell:
                '''
                head -1 {input[0]} > {output[0]}
                tail -n +2 -q {input} >> {output[0]}
                '''

rule gzip_results_non_rel:
        'Gzip results.'
        input:
                'results/GWAS/sumstats_non_rel/GWAS-{pheno}/temp/allchr/{sample}.txt'
        output:
                'results/GWAS/sumstats_non_rel/GWAS-{pheno}/allchr-{sample}.txt.gz'
        shell:
                'gzip -c {input[0]} > {output[0]}'

rule check_results_non_rel:
        ''
        input:
                expand('results/GWAS/sumstats_non_rel/GWAS-{pheno}/allchr-{sample}.txt.gz', pheno= 'miscarriage', sample= fam_ids['fam_id'])
        output:
                'results/GWAS/checks_non_rel/GWAS_performed.txt'
        shell:
                'touch {output[0]}'

rule REGENIE_step2_non_rel_lm:
        'Whole genome regression model is fit to the traits.'
        input:
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/bgen/{CHR}.bgen',
                'results/pheno_non_rel/{sample}_phenotype_{pheno}.txt',
                'results/pheno_non_rel/{sample}_covars_{pheno}.txt',
                'results/aux/ids_non_rel/samples/{sample}_ids_{pheno}.txt',
                'results/aux/ids/samples/bgen/{CHR}_samples.txt',
        output:
                temp('results/GWAS/regenie_non_rel_lm/step2/temp/{sample}/{CHR}_results_{pheno}.regenie')
        params:
                'results/GWAS/regenie_non_rel_lm/step2/temp/{sample}/{CHR}_results',
                'results/GWAS/pgen/X'
	threads: 5
	run:
		if (wildcards.pheno != 'micsarriage_bin'):
			shell('''
                /home/pol.sole.navais/soft/regenie_v3.2.1.gz_x86_64_Linux \
                --step 2 \
                --bgen {input[0]} \
                --covarFile {input[2]} \
                --keep {input[3]} \
                --phenoFile {input[1]} \
                --bsize 400 \
                --firth --approx \
                --minINFO 0.6 \
                --threads {threads} \
                --sample {input[4]} \
                --force-qt \
		--ignore-pred \
                --out {params[0]} \
                --catCovarList cohort \
                --verbose
                               ''')
		else:
			shell('''
                /home/pol.sole.navais/soft/regenie_v3.2.1.gz_x86_64_Linux \
                --step 2 \
                --bgen {input[0]} \
                --covarFile {input[2]} \
                --keep {input[3]} \
                --phenoFile {input[1]} \
                --bsize 400 \
		--bt \
                --firth --approx \
                --minINFO 0.6 \
                --threads {threads} \
                --sample {input[4]} \
                --af-cc \
                --ignore-pred \
                --out {params[0]} \
                --catCovarList cohort \
                --verbose
                               ''')

rule concat_GWAS_results_non_rel_lm:
        'Concatenate results from regenie'
        input:
                expand('results/GWAS/regenie_non_rel_lm/step2/temp/{{sample}}/{CHR}_results_{{pheno}}.regenie', CHR= CHROM)
        output:
                temp('results/GWAS/sumstats_non_rel_lm/GWAS-{pheno}/temp/allchr/{sample}.txt')
        shell:
                '''
                head -1 {input[0]} > {output[0]}
                tail -n +2 -q {input} >> {output[0]}
                '''

rule gzip_results_non_rel_lm:
        'Gzip results.'
        input:
                'results/GWAS/sumstats_non_rel_lm/GWAS-{pheno}/temp/allchr/{sample}.txt'
        output:
                'results/GWAS/sumstats_non_rel_lm/GWAS-{pheno}/allchr-{sample}.txt.gz'
        shell:
                'gzip -c {input[0]} > {output[0]}'

rule check_results_non_rel_lm:
        ''
        input:
                expand('results/GWAS/sumstats_non_rel_lm/GWAS-{pheno}/allchr-{sample}.txt.gz', pheno= pheno_file['phenotypes'], sample= fam_ids['fam_id'])
        output:
                'results/GWAS/checks_non_rel_lm/GWAS_performed.txt'
        shell:
                'touch {output[0]}'

