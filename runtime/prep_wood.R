#####################################################################################
## Author: Sina Rueeger [sina *.* rueger *a*t* unil *.* ch]
## Project: 
##        
## Time-stamp: <[prep_wood.R] by SR Wed 21/06/2017 19:03 (CEST)>
##
## Description:
## 
##
## History:
## 
#####################################################################################

## import wood
## -------
library(readr)
path.raw <- "/data/sgg/sina/public.data/meta.summaries/meta.giant/"
tmp <- read_tsv(paste0(path.raw,
                           "GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt"))

tmp$Z <- tmp$b / tmp$SE
## name wood with MarkerName, Z, a1 and a2,
## Allele1 is the effect allele
vec.giant.full <- c("MarkerName", "a2", "a1", "EAF", "b", "SE", "p", "N","Z")
names(tmp) <- vec.giant.full

## write out
write_tsv(tmp, path = paste0(path.raw, "wood.txt"))
