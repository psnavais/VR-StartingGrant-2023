rule HWE:
	''
	input:
		'results/effect_origin/DS/miscarriage/fets_DS{CHR}',
		'results/pheno_non_rel/fets_phenotype_miscarriage.txt'
	output:
		'results/HWE/miscarriage/temp/fets_DS{CHR}'
	conda: '../envs/plots.yml'
	script:
		'../scripts/HWE.R'

rule concat_HWE:
        'Collect DS from all CHR.'
        input:
                expand('results/HWE/miscarriage/temp/fets_DS{CHR}', CHR= CHROM)
        output:
                'results/HWE/miscarriage/fets_HWE.txt'
        shell:
                '''
                set +o pipefail;
                head -1 {input[0]} > {output[0]}
                cat {input} | grep -v 'SNP' >> {output[0]}
                '''

