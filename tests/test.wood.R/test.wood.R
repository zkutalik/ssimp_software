rm(list = ls())


## set paths
## -------------------
path.outdir = Sys.getenv('COMMANDDIR')

## load function
## -------------------
source("tests/r-function.R")



## run ssimp function
## ---------------
ssimp(path.gwas = "gwas/GIANT_HEIGHT_Wood_somechr22snps.txt.appendZcolumn.txt",
      path.ref  = "ref/1000genomes_somechr22snps_EURsamples.vcf",
      path.outdir = paste0(path.outdir , '/'),
      lambda = "sqrt",
      what.to.impute = "Z",
      gwas.names = c("MarkerName", NA, NA, "Allele1","Allele2",NA,NA, 'Z', "N"),
      tag.snps = c("rs136382","rs136383", "rs136389"),
      #target.snps = "rs16987627", #it's biallelic, so we'll ignore it
      target.snps = c("rs5753220","rs5753231","rs5753236","rs5753259","rs5753260","rs5753263","rs5753268","rs5753271","rs5753272","rs5753281","rs5753284","rs5753285","rs5753290","rs5753298"),
      impute.tags = FALSE
      )



