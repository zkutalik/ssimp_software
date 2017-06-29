#####################################################################################
## Author: Sina Rueeger [sina *.* rueger *a*t* unil *.* ch]
## Project: 
##        
## Time-stamp: <[liftover_dbsnp.R] by SR Thu 29/06/2017 15:09 (CEST)>
##
## Description:
## 
##
## History:
## 
#####################################################################################

## checkout /data/sgg/sina/software/liftOver/howto_liftOver.txt

data.dbsnp <- "/data/sgg/sina/public.data/dbsnp/"
liftover <- "/data/sgg/sina/software/liftOver/liftOver"

## file to lift over
## ------------------
library(readr)
tmp <- read.table(paste0(data.dbsnp,"b150_SNPChrPosOnRef_108.bcp.gz"), nrow =10,
                  header = FALSE, sep="\t")

toyex <- data.frame(matrix(NA, ncol = 4, nrow = nrow(tmp)))
toyex[,1] <- paste0("chr",tmp$V2) ## chr#
toyex[,2] <- tmp$V3 ## pos-1
toyex[,3] <- tmp$V3+1L ## pos
toyex[,4] <- paste0("rs",tmp$V1) ## rs#
toyex[,5] <- tmp$V2  ## CHR

write_tsv(toyex, path = paste0(data.dbsnp, "toy_hg20.bed"), col_names = FALSE)



## liftover hg20 > hg19
## ---------------------

system(paste0(liftover, " ",data.dbsnp,"toy_hg20.bed /data/sgg/sina/software/liftOver/hg20ToHg19.over.chain ",data.dbsnp,"output_hg19.bed ",data.dbsnp,"unlifted_hg19.bed"))

## liftover hg20 > hg18
## ---------------------
system(paste0(liftover," ",data.dbsnp,"output_hg19.bed /data/sgg/sina/software/liftOver/hg19ToHg18.over.chain ",data.dbsnp,"output_hg18.bed ",data.dbsnp,"unlifted_hg18.bed"))


## merge files together
## --------------------
tmp.18 <- (read_tsv(paste0(data.dbsnp,"output_hg18.bed"),col_names = FALSE))[,c(3:5)]
tmp.19 <- (read_tsv(paste0(data.dbsnp,"output_hg19.bed"),col_names = FALSE))[,c(3:5)]
tmp.20 <- (read_tsv(paste0(data.dbsnp,"toy_hg20.bed"),col_names = FALSE))[,c(3:5)]
names(tmp.18) <- c("Pos.hg18","rsid","Chr")
names(tmp.19) <- c("Pos.hg19","rsid","Chr")
names(tmp.20) <- c("Pos.hg20","rsid","Chr")

out <- merge(tmp.18, merge(tmp.19[,-3], tmp.20[,-3], by = "rsid"), by = "rsid")
out <- out[,c("rsid", "Chr","Pos.hg18", "Pos.hg19", "Pos.hg20")]
