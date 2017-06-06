rm(list = ls())


## set paths
## -------------------
path.out <- "/Users/admin/Documents/Work/Projects/tagging/ssimp_software/tests/test.wood.R/"
setwd(path.out)

## load function
## -------------------
source("../r-function.R")



## run ssimp function
## ---------------
ssimp(path.gwas = paste0("/Users/admin/Documents/Work/Projects/tagging/ssimp_software/gwas/", "GIANT_HEIGHT_Wood_somechr22snps.txt"),
      path.ref = paste0("/Users/admin/Documents/Work/Projects/tagging/ssimp_software/ref/","1000genomes_somechr22snps_EURsamples.vcf"), 
      lambda = "sqrt", 
      what.to.impute = "Z", 
      gwas.names = c("MarkerName", NA, NA, "Allele1","Allele2","b", "SE", NA, "N"),
      tag.snps = c("rs136382","rs136383", "rs136389"),
      target.snps = "rs16987627",
      impute.tags = FALSE
      )



