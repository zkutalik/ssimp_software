#####################################################################################
## Author: Sina Rueeger [sina *.* rueger *a*t* unil *.* ch]
## Project: 
##        
## Time-stamp: <[liftover_dbsnp.R] by SR Thu 29/06/2017 15:57 (CEST)>
##
## Description:
## 
##
## History:
## 
#####################################################################################

## checkout /data/sgg/sina/software/liftOver/howto_liftOver.txt
rm(list=ls())

data.dbsnp <- "/data/sgg/sina/public.data/dbsnp/"
liftover <- "/data/sgg/sina/software/liftOver/liftOver"

## file to lift over
## ------------------
library(readr)
#zcat b150_SNPChrPosOnRef_108.bcp.gz  | egrep -v NotOn | awk '{print "chr"$2"\t"$3"\t"($3+1)"\trs"$1"\t"$2}' > dbsnp_hg20.bed


## liftover hg20 > hg19
## ---------------------

system(paste0(liftover, " ",data.dbsnp,"dbsnp_hg20.bed /data/sgg/sina/software/liftOver/hg20ToHg19.over.chain ",data.dbsnp,"liftover_hg19.bed ",data.dbsnp,"unlifted_hg19.bed"))

## liftover hg20 > hg18
## ---------------------
system(paste0(liftover," ",data.dbsnp,"liftover_hg19.bed /data/sgg/sina/software/liftOver/hg19ToHg18.over.chain ",data.dbsnp,"liftover_hg18.bed ",data.dbsnp,"unlifted_hg18.bed"))

## merge files together
## --------------------
library(data.table)
hg18 <- fread(paste0(data.dbsnp, "liftover_hg18.bed"),select = 3:5)
hg19 <- fread(paste0(data.dbsnp, "liftover_hg19.bed"),select = 3:5)
hg20 <- fread(paste0(data.dbsnp, "dbsnp_hg20.bed"),select = 3:5)
names(hg18) <- c("Pos.hg18","rsid","Chr")
names(hg19) <- c("Pos.hg19","rsid","Chr")
names(hg20) <- c("Pos.hg20","rsid","Chr")

setkey(hg18,rsid,Chr)
setkey(hg19,rsid,Chr)
setkey(hg20,rsid,Chr)
comb = merge(hg18, merge(hg19, hg20, all = TRUE, sort = FALSE) , all = TRUE, sort = FALSE) 





















## OLD

tmp.18 <- (read_tsv(paste0(data.dbsnp,"output_hg18.bed"),col_names = FALSE))[,c(3:5)]
tmp.19 <- (read_tsv(paste0(data.dbsnp,"output_hg19.bed"),col_names = FALSE))[,c(3:5)]
tmp.20 <- (read_tsv(paste0(data.dbsnp,"toy_hg20.bed"),col_names = FALSE))[,c(3:5)]
names(tmp.18) <- c("Pos.hg18","rsid","Chr")
names(tmp.19) <- c("Pos.hg19","rsid","Chr")
names(tmp.20) <- c("Pos.hg20","rsid","Chr")

out <- merge(tmp.18, merge(tmp.19[,-3], tmp.20[,-3], by = "rsid"), by = "rsid")
out <- out[,c("rsid", "Chr","Pos.hg18", "Pos.hg19", "Pos.hg20")]


tmp <- read.table(paste0(data.dbsnp,"b150_SNPChrPosOnRef_108.bcp.gz"), nrow =10,
                  header = FALSE, sep="\t")

toyex <- data.frame(matrix(NA, ncol = 4, nrow = nrow(tmp)))
toyex[,1] <- paste0("chr",tmp$V2) ## chr#
toyex[,2] <- tmp$V3 ## pos-1
toyex[,3] <- tmp$V3+1L ## pos
toyex[,4] <- paste0("rs",tmp$V1) ## rs#
toyex[,5] <- tmp$V2  ## CHR

write_tsv(toyex, path = paste0(data.dbsnp, "toy_hg20.bed"), col_names = FALSE)
