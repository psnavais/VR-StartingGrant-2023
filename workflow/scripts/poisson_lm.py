import shutil
import hail as hl
hl.init()

from hail.plot import show
from pprint import pprint

hl.index_bgen(path= snakemake.input[2], index_file_map= {snakemake.input[2]: snakemake.params[0]})

hl.import_bgen(snakemake.input[2] , entry_fields=['dosage'], sample_file= snakemake.input[3], index_file_map= {snakemake.input[2]: snakemake.params[0]}).write(snakemake.params[1], overwrite=True)

table = hl.import_table(snakemake.input[0], impute=True).key_by('IID')

table2 = hl.import_table(snakemake.input[1], impute=True).key_by('IID')

mt = hl.read_matrix_table(params[1])

mt = mt.annotate_cols(pheno = table[mt.s])
mt = mt.annotate_cols(covars = table[mt.s])

gwas = hl.poisson_regression_rows(y= mt.pheno.miscarriage,x= mt.dosage, covariates= [1.0, mt.covars.MOR_FAAR, mt.covars.PC1, mt.covars.PC2, mt.covars.PC3, mt.covars.PC4, mt.covars.PC5, mt.covars.PC6, mt.covars.PC7, mt.covars.PC8, mt.covars.PC9, mt.covars.PC10])

gwas.export(snakemake.output[0], delimiter= '\t', header = True)

shutil.rmtree(snakemake.params[0])
shutil.rmtree(snakemake.params[1])
