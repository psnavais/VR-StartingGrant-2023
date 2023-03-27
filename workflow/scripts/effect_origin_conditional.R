library(data.table)
library(dplyr)
library(tidyr)
library(broom)
library(MASS)

format_haps= function(hap){
variants= paste('X', hap$chr, hap$pos, hap$ref, hap$eff, sep ='_')
ids= names(hap)[5:ncol(hap)]
hap= as.data.frame(t(hap[, 5:ncol(hap)]))
names(hap)=  variants
hap$IID= ids

return(hap)
}

fets= fread(snakemake@input[[1]])

fets= format_haps(fets)

pheno= fread(snakemake@input[[2]])
covar= fread(snakemake@input[[3]])

pheno= inner_join(pheno, covar, by= 'IID') 

fets= inner_join(fets, pheno, by= 'IID')

print(nrow(fets))

results_list= lapply(names(fets)[grepl('X_', names(fets))], function(snp) {

print(snp)
fets_temp= data.frame(fets)[, c(snp, 'miscarriage', 'MOR_FAAR', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'cohort')]

names(fets_temp)[1]= 'fet'

m1= glm(miscarriage~ fet + cohort + PC1 + PC2 + PC3 + PC4 + PC5 + PC6, data= fets_temp, family= poisson)

temp_df= tidy(m1) %>% filter(row_number() == 2)
temp_df$term= snp
temp_df$model= 'poisson'
temp_df$loglik= logLik(m1)[1]

m2= glm.nb(miscarriage~ fet + cohort + PC1 + PC2 + PC3 + PC4 + PC5 + PC6, data= fets_temp)

temp_df2= tidy(m2) %>% filter(row_number() == 2)
temp_df2$term= snp
temp_df2$model= 'negative-binomial'
temp_df2$loglik= logLik(m2)[1]

temp_df= bind_rows(temp_df, temp_df2)

return(temp_df)

}

)

print('Analyses performed, saving data.')

x= do.call('rbind', results_list)

fwrite(x, snakemake@output[[1]], sep= '\t')
