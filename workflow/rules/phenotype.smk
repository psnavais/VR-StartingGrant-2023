import pandas as pd
import numpy as np
import csv
import gzip

##### load config and sample sheets #####

configfile: "config/config.yaml"

fam_ids = pd.read_csv(config["fam_member"], sep= "\t").set_index("fam_id", drop=False)
pheno_file= pd.read_csv(config["pheno"], sep= '\t')
CHROM= [i.strip() for i in open(config["CHROM"], 'r')]

rule phenofile:
        'Pheno file for miscarriage.'
        input:
                '/mnt/work/p1724/v12/PDB1724_MFR_541_v12.csv',
                '/mnt/work/pol/MoBaGenetics-1.0/delivery/sample-ids.txt',
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/aux/flaglist-merged/mobagen-flaglist-n99259.txt',
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/aux/pca/ethnic-core-samples'
        output:
                temp('results/aux/pheno/temp/{pheno}_temp_moms.txt'),
                temp('results/aux/pheno/temp/{pheno}_temp_fets.txt'),
                temp('results/aux/pheno/temp/{pheno}_temp_dads.txt')
	conda:
		'../envs/plots.yml'
	script:
		'../scripts/pheno_file.R'

rule concat_phenos_PCA:
        'Concat pheno files, and add PCA.'
        input:
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/aux/pca/mobagen-total/mobagen-total-proj-pc',
                'results/aux/pheno/temp/{pheno}_temp_{sample}.txt',
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/aux/pedigree/mobagen-ethnic-core-samples.kin0'
        output:
                'results/pheno/{sample}_phenotype_{pheno}.txt',
                'results/pheno/{sample}_covars_{pheno}.txt',
                'results/aux/ids/samples/{sample}_ids_{pheno}.txt',
        run:
                d= pd.read_csv(input[1], header= 0, sep= '\t')
                pca= pd.read_csv(input[0], header= 0, sep= '\t')
                d= pd.merge(d, pca, how= 'inner', on= 'IID')
                d.fillna('NA', inplace= True)
                for elem in set(d['cohort']):
                        d[str(elem)]= d['cohort'] == elem
                d[['HARVEST', 'NORMENT', 'ROTTERDAM1', 'ROTTERDAM2']] *= 1
                d['FID']= d.IID
                d.to_csv(output[0], sep= '\t', header= True, index= False, columns= ['FID', 'IID', wildcards.pheno])
                d.to_csv(output[1], sep= '\t', header= True, index= False, columns= ['FID', 'IID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 'cohort', 'MOR_FAAR'])
                d.to_csv(output[2], sep= '\t', header= False, index= False, columns= ['FID', 'IID'])

rule phenofile_non_rel:
        'Pheno file for miscarriage.'
        input:
                '/mnt/work/p1724/v12/PDB1724_MFR_541_v12.csv',
                '/mnt/work/pol/MoBaGenetics-1.0/delivery/sample-ids.txt',
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/aux/flaglist-merged/mobagen-flaglist-n99259.txt',
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/aux/pca/ethnic-core-samples'
        output:
                temp('results/aux/pheno_non_rel/temp/{pheno}_temp_moms.txt'),
                temp('results/aux/pheno_non_rel/temp/{pheno}_temp_fets.txt'),
                temp('results/aux/pheno_non_rel/temp/{pheno}_temp_dads.txt')
        conda:
                '../envs/plots.yml'
        script:
                '../scripts/pheno_file.R'

rule concat_phenos_PCA_non_rel:
        'Concat pheno files, and add PCA.'
        input:
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/aux/pca/mobagen-total/mobagen-total-proj-pc',
                'results/aux/pheno_non_rel/temp/{pheno}_temp_{sample}.txt',
                '/mnt/archive/moba/geno/MOBAGENETICS_1.0/genotypes-base/aux/pedigree/mobagen-ethnic-core-samples.kin0'
        output:
                'results/pheno_non_rel/{sample}_phenotype_{pheno}.txt',
                'results/pheno_non_rel/{sample}_covars_{pheno}.txt',
                'results/aux/ids_non_rel/samples/{sample}_ids_{pheno}.txt'
	run:
                d= pd.read_csv(input[1], header= 0, sep= '\t')
                pca= pd.read_csv(input[0], header= 0, sep= '\t')
                d= pd.merge(d, pca, how= 'inner', on= 'IID')
                remove= selectUnrelated(input[2], d, d.IID)
		d= d.loc[~d.IID.isin(remove.IID), : ]
		d.fillna('NA', inplace= True)
                for elem in set(d['cohort']):
                        d[str(elem)]= d['cohort'] == elem
                d[['HARVEST', 'NORMENT', 'ROTTERDAM1', 'ROTTERDAM2']] *= 1
                d['FID']= d.IID
                d.to_csv(output[0], sep= '\t', header= True, index= False, columns= ['FID', 'IID', wildcards.pheno])
                d.to_csv(output[1], sep= '\t', header= True, index= False, columns= ['FID', 'IID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 'cohort', 'MOR_FAAR'])
                d.to_csv(output[2], sep= '\t', header= False, index= False, columns= ['FID', 'IID'])
