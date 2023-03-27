
rule hail_GWAS:
	''
	input:
		'results/pheno_non_rel/{sample}_phenotype_miscarriage.txt',
		'results/pheno_non_rel/{sample}_covars_miscarriage.txt',
		'/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/imputed/all/bgen/{CHR}.bgen',
		'results/aux/ids/samples/bgen/{CHR}_samples.txt'
	output:
		'results/GWAS/poisson/miscarriage/{sample}-poisson-{CHR}.txt'
	params:
		'results/GWAS/poisson/temp/idx/{CHR}.idx2',
		'results/GWAS/poisson/temp/matrix-tables/{CHR}.mt'
	threads: 30
	script:
		'../scripts/poisson_lm.py'


rule concat_GWAS_results_poisson:
        'Concatenate results from regenie'
        input:
                expand('results/GWAS/poisson/miscarriage/{{sample}}-poisson-{CHR}.txt', CHR= CHROM)
        output:
                temp('results/GWAS/poisson/sumstats/GWAS-miscarriage/temp/allchr/{sample}.txt')
        shell:
                '''
                head -1 {input[0]} > {output[0]}
                tail -n +2 -q {input} >> {output[0]}
                '''

rule gzip_results_poisson:
        'Gzip results.'
        input:
                'results/GWAS/poisson/sumstats/GWAS-miscarriage/temp/allchr/{sample}.txt'
        output:
                'results/GWAS/poisson/sumstats/GWAS-miscarriage/allchr-{sample}.txt.gz'
        shell:
                'gzip -c {input[0]} > {output[0]}'

rule check_results_poisson:
        ''
        input:
                expand('results/GWAS/poisson/sumstats/GWAS-miscarriage/allchr-{sample}.txt.gz', sample= 'fets')
        output:
                'results/GWAS/poisson/checks/GWAS_performed.txt'
        shell:
                'touch {output[0]}'

