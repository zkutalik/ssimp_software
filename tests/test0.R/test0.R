rm(list = ls())


## set paths
## -------------------

## load function
## -------------------
source("../r-function.R")



## run ssimp function
## ---------------
ssimp(path.gwas = paste0("../../gwas/", "small.UKBB.both.csv"),
      path.ref = paste0("../../ref/","small.vcf.sample.vcf"), 
      path.outdir = './',
      lambda = 0.01,
      what.to.impute = "Z", 
      gwas.names = c("rnpid", NA, NA, "a2","a1",NA, NA, "z.from.peff", "N"),
      tag.snps = c("rs138257042", "rs200923174", "rs62224618", "rs79847867"),
      target.snps = c("rs587697622", "rs587755077","rs587654921"),
      impute.tags = FALSE
)



## compare with aarons results
## -----------------------------
out.aaron <- read_tsv("../test0/imputations.txt", skip = 1, n_max = 3, col_names = FALSE)
out.aaron

out.sina.R <- read_tsv("imputations.txt")
out.sina.R

