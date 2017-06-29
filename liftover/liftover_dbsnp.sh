#####################################################################################
## Author: Sina Rueeger [sina *.* rueger *a*t* unil *.* ch]
## Project: 
##        
## Time-stamp: <[liftover_dbsnp.R] by SR Thu 29/06/2017 17:02 (CEST)>
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
#zcat b150_SNPChrPosOnRef_108.bcp.gz  | egrep -v NotOn | egrep -v ChrUn | awk '{print "chr"$2"\t"$3"\t"($3+1)"\trs"$1"\t"$2}' > dbsnp_hg20.bed


## liftover hg20 > hg19
## ---------------------

system(paste0(liftover, " ",data.dbsnp,"dbsnp_hg20.bed /data/sgg/sina/software/liftOver/hg20ToHg19.over.chain ",data.dbsnp,"liftover_hg19.bed ",data.dbsnp,"unlifted_hg19.bed"))

## liftover hg20 > hg18
## ---------------------
system(paste0(liftover," ",data.dbsnp,"liftover_hg19.bed /data/sgg/sina/software/liftOver/hg19ToHg18.over.chain ",data.dbsnp,"liftover_hg18.bed ",data.dbsnp,"unlifted_hg18.bed"))

## merge files together
## --------------------
cd /data/sgg/sina/public.data/dbsnp
awk 'FNR==NR{a[$4]=$3;next}{if(a[$4]==""){a[$4]=0}; print $3,a[$4],$4,$5}' dbsnp_hg20.bed liftover_hg19.bed > tmp.txt

awk 'FNR==NR{a[$4]=$3;next}{if(a[$4]==""){a[$4]=0}; print a[$3],$1,$2,$3,$4}' liftover_hg18.bed tmp.txt > dbsnp_hg18_hg19_hg20.txt

rm tmp.txt


#FORMAT
#dbsnp_hg18_hg19_hg20.txt

#pos_hg18 pos_hg19 pos_hg20 rs  chr
#91617493 91779557 92150243 rs7 7
#...
