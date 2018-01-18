
## filter + copy
## //////////////

system("cp ../ref/sub1KG-tiny/chr22.vcf.gz ../ref/small.vcf.sample.vcf.gz")
system("cp ../ref/sub1KG-tiny/chr22.vcf.gz.tbi ../ref/small.vcf.sample.vcf.gz.tbi")



##create 
library(readr)
library(tidyverse)
dat <- read_tsv("../gwas/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz")
ref <- read_tsv("../ref/small.vcf.sample.vcf.gz")

dat <- left_join(dat, ref %>% select(ID,POS,1) %>% rename(MarkerName = ID), by = "MarkerName") %>% rename("CHR" = !!names(.[10])) %>% filter(CHR == 22)


small.random.n.txt <- dat
small.random.txt <- dat %>% mutate(Z = b/sqrt(SE))%>% select(MarkerName,Allele1,Allele2,Z) ## remove N/CHR/POS
small.random.p.b.txt <- dat %>% select(MarkerName,Allele1,Allele2,b,p) ## add p/b

small.random.chr.pos.txt <- dat %>% mutate(MarkerName = NA) %>% select(-N) 

write_delim(small.random.n.txt, path = "../gwas/small.random.n.txt",delim = " ")
write_delim(small.random.txt, path = "../gwas/small.random.txt",delim = " ")
write_delim(small.random.p.b.txt, path = "../gwas/small.random.p.b.txt",delim = " ")
write_delim(small.random.chr.pos.txt, path = "../gwas/small.random.chr.pos.txt",delim = " ")

set.seed(3)
listofimputesnps <- ref %>% select(ID) %>% sample_n(10)
listoftagsnps <- dat %>% select(MarkerName) %>% sample_n(10)

write_delim(, path = "../gwas/listofimputesnps.txt",delim = " ",col_names=FALSE)
write_delim(, path = "../gwas/listoftagsnps.txt",delim = " ",col_names=FALSE)
write_delim(, path = "../gwas/filename.samples.txt",delim = " ",col_names=FALSE)
write_delim(, path = "../gwas/filename.samples.txt",delim = " ",col_names=FALSE)


`listofimputesnps.txt` contains SNP id's separated by new lines (no header).
`listoftagsnps.txt` contains SNP id's separated by new lines (no header).
`filename.samples.txt` contains sample id's separated by new lines (no header). >> should be able to leave this
Here, `filename.samples.txt` contains sample id's (`sample`) along with a second attribute (here `super_pop`) that has different values, among them is `EUR`, for which we separate. 
