rm(list = ls())


## set paths
## -------------------
path.out <- "/Users/admin/Documents/Work/Projects/tagging/ssimp_software/tests/test0.R/"
setwd(path.out)

## load function
## -------------------
source("../r-function.R")



## run ssimp function
## ---------------
ssimp(path.gwas = paste0("/Users/admin/Documents/Work/Projects/tagging/ssimp_software/gwas/", "small.UKBB.both.csv"),
      path.ref = paste0("/Users/admin/Documents/Work/Projects/tagging/ssimp_software/ref/","small.vcf.sample.vcf"), 
      lambda = "sqrt", 
      what.to.impute = "Z", 
      gwas.names = c("rnpid", NA, NA, "a2","a1",NA, NA, "z.from.peff", "N"),
      tag.snps = c("rs138257042", "rs200923174", "rs62224618", "rs79847867"),
      target.snps = c("rs587697622", "rs587755077","rs587654921"),
      impute.tags = FALSE
)



## compare with aarons results
## -----------------------------
out.aaron <- read_tsv("/Users/admin/Documents/Work/Projects/tagging/ssimp_software/tests/test0/imputations.txt", skip = 1, n_max = 3, col_names = FALSE)
out.aaron

out<- read_tsv("/Users/admin/Documents/Work/Projects/tagging/ssimp_software/tests/test0.R/imputations.txt")
out

