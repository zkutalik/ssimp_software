#####################################################################################
## Author: Sina Rueeger [sina *.* rueger *a*t* unil *.* ch]
## Project: 
##        
## Time-stamp: <[modify.R] by SR Wed 22/03/2017 13:09 (CET)>
##
## Description:
## 
##
## History:
## 
#####################################################################################


n.tsnps <- 500

library(readr)

## random gwas example
## ---------------------
tmp <- read_delim("/data/sgg/sina/data/project.eqtl/data/metal/inputfile_Fehrmann_H8v2.txt",delim =
                  " ")
tmp
tmp <- tmp[1:n.tsnps,]

write_delim(tmp, "gwas_example.txt")

## METAL
## ---------------------
tmp <- read_tsv("/data/sgg/sina/data/project.eqtl/data/metal/METAANALYSIS1.TBL")
tmp
tmp <- tmp[1:n.tsnps,]

write_tsv(tmp, "metal_example.TBL")

## snptest imp
## ---------------------
tmp <-
    read_delim("/data/sgg/sina/data/project.imputation/ukbb.comparison/gwas_rs138273386_2017-03-08.out",delim
               = " ", comment = "#")
tmp
tmp <- tmp[1:n.tsnps,]

write_delim(tmp, "snptest_genotype_imp.out")
