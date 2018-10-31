
## We impute a selection of tag SNPs. We call this "reimpute". This is only done for the first window. 

## When imputing a tag SNP, this SNP is excluded from the list of tag SNPs. 

## Checking the imputation of tag SNPs with the actual Z-statistics, can be useful as a sanity check. 
## For example, if your alleles in the gwas file as mislabelled, you will see this here. 

## lets run an example

library(tidyverse)

## run ssimp
system("../compiled/ssimp-osx-0.2 --gwas ../gwas/small.random.txt --ref ../ref/small.vcf.sample.vcf.gz --out ../output.txt")

## load data
dat.imp <- read_tsv("../output.txt")
dat <- read_delim("../gwas/small.random.txt",delim=" ")

## merge
dat.merge <- left_join(dat.imp %>% filter(source == "GWAS") %>% select(SNP, Allele1, Allele2, z_imp, r2.pred, Z_reimputed, r2_reimputed), dat %>% rename(SNP = MarkerName, Z.GWAS = Z, Allele1.GWAS = Allele1, Allele2.GWAS=Allele2) , by = "SNP")
dat.merge <- dat.merge %>% mutate(Z.GWAS = ifelse(Allele1.GWAS != Allele1 & Allele2.GWAS != Allele2, (-1) * Z.GWAS, Z.GWAS)) ## we only swap the sign of the Z stats
  
## lets flip the alleles and Z stats
qplot(Z.GWAS, Z_reimputed, data = dat.merge, color = r2_reimputed) + geom_abline(intercept = 0, slope = 1) ## plot the identity line

## note, that z.imp is the imputation of the tag SNP, while leaving the SNP in the set of tag SNPs. hence r2.pred = 1.
