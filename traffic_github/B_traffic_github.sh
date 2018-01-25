#!/usr/bin/env Rscript




## ////////////////////
## (B) using the functions
## ////////////////////
sessionInfo()

args <- commandArgs(trailingOnly = TRUE)
ssimp.args <- (as.character(args[1])) ## your ssimp line
nbr.chunks <- as.numeric(as.character(args[2])) ## in how many chunks should the computation be split up?
out.name <- (as.character(args[3])) ## in how many chunks should the computation be split up?

cat(ssimp.args, "\n", nbr.chunks, "\n" , out.name,  "\n")

library(readr) ## yep - this one you need: to read and write data 
library(plyr) ## for ddply
library(dplyr) ## for filter

## download reference file that contains chr/pos
## ----------------------------------------------
url <-"https://drive.switch.ch/index.php/s/e3MuOYEudU1JbUo/download"
file <- "reference_panels/dbsnp_hg20_chr_pos_sorted.txt.gz"
if(!file.exists(file))
{
    download.file(url, file)
}
## file from my public.data folder, processed in the following way:
## 1) awk -F"\t" '{print $5"\t"$3}' dbsnp_hg20.bed > dbsnp_hg20_chr_pos.txt
## 2) remove all chr that are not 1:23 : cat dbsnp_hg20_chr_pos.txt | sed '/MT/d' > tmp_noMT.txt
## 2) sort -nk1,1 -nk 2,2 tmp_noMT.txt > dbsnp_hg20_chr_pos_sorted.txt  ## order pos and chr
## 3) gzip dbsnp_hg20_chr_pos_sorted.txt
## 4) copied to switch drive

cat("Start: loading large file\n")
tkg <- read_tsv(file, col_names = FALSE) ## this can take some time
cat("Finished: loading large file\n")
## columns
## X1 = chr
## X2 = pos - hg20

## remove MT and PAR and Y and X
tkg <- filter(tkg, X1 %in% 1:22)


## create file with ssimp chunks
## ----------------------------------
sessionInfo()

impute.range <- uk10k.chunks.from.to(ref.file=tkg, nbr.chunks = nbr.chunks)  ## returns "chr:pos.start-chr:pos.end"

print.chunks(ssimp.args = ssimp.args, nbr.chunks = nbr.chunks, ref.file = tkg, out.name=out.name)
cat("cat2\n")



## run this file with
## chmod a+x ssimp_chunks.sh
## ./ssimp_chunks.sh "bin/ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt" 100 run_ssimp_array.sh

## then integrate your shell script [out.name] into your scheduler.

## main R shell script skeleton from:
## https://github.com/gastonstat/tutorial-R-noninteractive/blob/master/04-shell-script.Rmd

