library(data.table)
library(dplyr)
library(tidyr)

mfr= fread(snakemake@input[[1]])
ids= fread(snakemake@input[[2]])


flag= fread(snakemake@input[[3]])
flag= filter(flag, genotypesOK== T, phenoOK== T)

out= readLines(snakemake@input[[4]])

ids= pivot_longer(ids, c('Child', 'Mother', 'Father'), names_to= 'ROLE', values_to= 'IID')
ids= ids[!duplicated(ids[,c('PREG_ID_1724', 'IID')]), ]

ids= filter(ids, (IID %in% out), IID %in% flag$IID)
ids= group_by(ids, PREG_ID_1724, ROLE) %>% filter(row_number()== 1)

ids= spread(ids, key= ROLE, value= IID)

mfr= inner_join(mfr, ids, by= 'PREG_ID_1724')

mfr$miscarriage= with(mfr, ifelse(is.na(SPABORT_12_5), SPABORT_23_5, ifelse(is.na(SPABORT_23_5), SPABORT_12_5, SPABORT_12_5 + SPABORT_23_5)))

mfr$miscarriage_bin= with(mfr, ifelse(is.na(miscarriage), NA, ifelse(miscarriage> 1, 1, 0)))

mfr2= select(mfr, PREG_ID_1724, miscarriage, PARITET_5, Child, Mother, Father)
mfr2= filter(mfr2, !is.na(miscarriage), !is.na(PARITET_5))
mfr2$miscarriage_resid= lm(miscarriage ~ PARITET_5, mfr2)$resid
mfr2= select(mfr2, PREG_ID_1724, miscarriage_resid, Child, Mother, Father)
mfr= left_join(mfr, mfr2, by= c('PREG_ID_1724', 'Child', 'Mother', 'Father'))

mfr= arrange(mfr, desc(miscarriage))

if (grepl('bin', snakemake@output[[1]])){

moms= filter(mfr, !duplicated(PREG_ID_1724)) %>% select(Mother, miscarriage_bin, PREG_ID_1724, MOR_FAAR, BATCH_moms) %>% filter(!is.na(Mother), !is.na(miscarriage_bin))
fets= filter(mfr, !duplicated(PREG_ID_1724)) %>% select(Child, miscarriage_bin, PREG_ID_1724, MOR_FAAR, BATCH_fets) %>% filter(!is.na(Child), !is.na(miscarriage_bin))
dads= filter(mfr, !duplicated(PREG_ID_1724)) %>% select(Father, miscarriage_bin, PREG_ID_1724, MOR_FAAR, BATCH_dads) %>% filter(!is.na(Father), !is.na(miscarriage_bin))

} else if (grepl('_resid', snakemake@output[[1]])) {

moms= filter(mfr, !duplicated(PREG_ID_1724)) %>% select(Mother, miscarriage_resid, PREG_ID_1724, MOR_FAAR, BATCH_moms) %>% filter(!is.na(Mother), !is.na(miscarriage_resid))
fets= filter(mfr, !duplicated(PREG_ID_1724)) %>% select(Child, miscarriage_resid, PREG_ID_1724, MOR_FAAR, BATCH_fets) %>% filter(!is.na(Child), !is.na(miscarriage_resid))
dads= filter(mfr, !duplicated(PREG_ID_1724)) %>% select(Father, miscarriage_resid, PREG_ID_1724, MOR_FAAR, BATCH_dads) %>% filter(!is.na(Father), !is.na(miscarriage_resid))

} else {

moms= filter(mfr, !duplicated(PREG_ID_1724)) %>% select(Mother, miscarriage, PREG_ID_1724, MOR_FAAR, BATCH_moms) %>% filter(!is.na(Mother), !is.na(miscarriage))
fets= filter(mfr, !duplicated(PREG_ID_1724)) %>% select(Child, miscarriage, PREG_ID_1724, MOR_FAAR, BATCH_fets) %>% filter(!is.na(Child), !is.na(miscarriage))
dads= filter(mfr, !duplicated(PREG_ID_1724)) %>% select(Father, miscarriage, PREG_ID_1724, MOR_FAAR, BATCH_dads) %>% filter(!is.na(Father), !is.na(miscarriage))

}

moms= moms[!duplicated(moms$Mother, incomparables= NA), ]
fets= fets[!duplicated(fets$Child, incomparables= NA), ]
dads= dads[!duplicated(dads$Father, incomparables= NA), ]

names(moms)[ncol(moms)]= 'cohort'
names(fets)[ncol(fets)]= 'cohort'
names(dads)[ncol(dads)]= 'cohort'

moms$cohort= ifelse(grepl('NORM', moms$cohort), 'NORMENT', moms$cohort)
dads$cohort= ifelse(grepl('NORM', dads$cohort), 'NORMENT', dads$cohort)
fets$cohort= ifelse(grepl('NORM', fets$cohort), 'NORMENT', fets$cohort)

names(moms)[1]= 'IID'
names(fets)[1]= 'IID'
names(dads)[1]= 'IID'

fwrite(moms, snakemake@output[[1]], sep= '\t')
fwrite(fets, snakemake@output[[2]], sep= '\t')
fwrite(dads, snakemake@output[[3]], sep= '\t')
