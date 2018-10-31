## This file creates the small toy datasets used in example.md


## filter + copy
## //////////////

system("cp ../ref/sub1KG-tiny/chr22.vcf.gz ../ref/small.vcf.sample.vcf.gz")
system("cp ../ref/sub1KG-tiny/chr22.vcf.gz.tbi ../ref/small.vcf.sample.vcf.gz.tbi")



## create GWAS
## -------------
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



## chr x example
## ----------------
ref <- read_tsv("../ref/sub1KG-tiny/chrX.vcf.gz")

map.x <- ref %>% select(1,2,3,4,5) %>% rename("CHR" = !!names(.[1])) 

set.seed(4)
gwas.x <- map.x %>%  sample_n(nrow(ref)/3*2) %>% mutate(Z = runif(nrow(.), -4,4))
target.id <- setdiff( map.x$ID, gwas$ID)
target.pos <- setdiff( map.x$POS, gwas$POS)

write_delim(gwas.x, path = "../gwas/small.random.x.txt",delim = " ")


## listof...
## -------------
set.seed(3)
listofimputesnps <- ref %>% select(ID) %>% sample_n(10)
listoftagsnps <- dat %>% select(MarkerName) %>% sample_n(10)

write_delim(, path = "../gwas/listofimputesnps.txt",delim = " ",col_names=FALSE)
write_delim(, path = "../gwas/listoftagsnps.txt",delim = " ",col_names=FALSE)
write_delim(, path = "../gwas/filename.samples.txt",delim = " ",col_names=FALSE)
write_delim(, path = "../gwas/filename.samples.txt",delim = " ",col_names=FALSE)

## 
## -------------
sample.small <- read_table("../ref/filename.samples.small.txt",col_names=FALSE)
sample.pop <- read_delim("../ref/filename.samples.txt",delim="\t")

## make half european/gbr
sample.pop.small <- sample.pop %>% filter(sample %in% sample.small$X1) %>% select(sample, pop,super_pop, gender)
sample.pop.small[1:7,"pop"] <- "GBR"
sample.pop.small[1:7,"super_pop"] <- "EUR"

write_delim(sample.pop.small, path = "../gwas/filename.samples.pop.txt", delim = "\t")

