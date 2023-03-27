rule miscarriage_fets:
	''
	input:
		'results/pheno_non_rel/fets_phenotype_miscarriage.txt',
		'results/GWAS/sumstats_non_rel_lm/GWAS-miscarriage/allchr-fets.txt.gz',
		'results/HWE/miscarriage/fets_HWE.txt'
	output:
		'reports/miscarriage_fets.html'
	conda:
		'../envs/reports.yml'
	script:
		'../scripts/fetal_miscarriage.Rmd'
