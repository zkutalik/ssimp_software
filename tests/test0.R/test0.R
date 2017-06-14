rm(list = ls())


## set paths
## -------------------

## load function
## -------------------
source("../r-function.R")



## run ssimp function
## ---------------
## a1 is reference allele and a2 alternative allele
ssimp(path.gwas = paste0("../../gwas/", "small.UKBB.both.csv"),
      path.ref = paste0("../../ref/","small.vcf.sample.vcf"), 
      path.outdir = './',
      lambda = "sqrt",
      what.to.impute = "Z", 
      gwas.names = c("rnpid", NA, NA, "a1","a2",NA, NA, "z.from.peff", "N"),
      tag.snps = c(
                   #, "rs200923174" # just used by SR
                     "rs62224618"
                   ,"rs138257042"
                   #, "rs79847867" # why wasn't I used this?
                   ),
      target.snps = c("rs587697622", "rs587755077","rs587654921"),
      impute.tags = TRUE
)



## compare with aarons results
## -----------------------------
out.aaron <- read_tsv("../test0/imputations.txt", skip = 1, n_max = 3, col_names = FALSE)

out.sina.R <- read_tsv("imputations.txt")
out.aaron
out.sina.R

## compare imputation of tags
## -----------------------------
out <- read_tsv("imputations_tag.txt")
out
