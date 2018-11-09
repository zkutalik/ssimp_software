#!/usr/bin/env Rscript

# These ranges can be computed automatically, by giving your main arguments [`args1`], the number of chunks [`args2`] and the name of the shell script returned [`args3`]. 
# 
# `./ssimp_chunks.sh "args1" args2 args3`
# 
# For example, run this in your terminal:
# `./ssimp_chunks.sh "bin/ssimp --gwas gwas/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz --ref ref/sub1KG-tiny/chr22.vcf.gz --out output.txt" 100 run_ssimp_array.sh`
# 
# Note, that you still need to integrate the shell script content into your scheduler.
# 

## (A) two functions
## ////////////////////

## impute.range
## ----------------------------------
## function that calculates the start
## and end positions for a number of
## chunks needed (nbr.chunks)
## ----------------------------------
uk10k.chunks.from.to <- function(ref.file = NULL, nbr.chunks = 1) ## returns "chr:pos.start-chr:pos.end"
{
  
  ## estimate the fraction each chromosome occupies
  #frac.chr <- ddply(ref.file, .(X1), function(x) nrow(x)/nrow(ref.file))

  ## translate the fraction into number of chunks per chromosome
  #frac.chr <- mutate(frac.chr, nbr.chunks = round(V1 * nbr.chunks))
  
  cat(nrow(ref.file), "\n")
  cat(nbr.chunks, "\n")
  
  ## build chunks
  ind <- round(seq(1, nrow(ref.file), length.out = nbr.chunks)) ## make sure its integers
  ind[length(ind)] <- nrow(ref.file)+1 ## replace by nrow, just to be sure that all pos are included from 1:nrow
  cat(paste(ind, collapse = " "), "\n")
  
  from <- ref.file[ind[-length(ind)],] %>% plyr::summarise(vector=paste(X1, X2, sep=":"))
  to <- ref.file[ind[-1]-1,] %>% plyr::summarise(vector=paste(X1, X2, sep=":")) ## chr:(pos-1)
  from.to <- paste0(c(from, recursive = TRUE), "-", c(to, recursive = TRUE))
  return(from.to) ## vector with "chr:pos.start-chr:pos.end" entries
  
}


## ssimp.chunks
## -----------------------------------------
## writes a shell script that has nbr.chunks
## lines, covering the whole genome
## -----------------------------------------

print.chunks <- function(ssimp.args = "fill me", nbr.chunks = 1, out.name = NULL, ref.file = NULL)
{
  ## ssimp.args: your ssimp line
  ## nbr.chunks: in how many chunks should the computation be split up?
  ## out.name:   path of the out file
  ## ref.file:   name of ref file (tkg)

  ## calculate chunks according to nbr.chunks
  impute.range <- uk10k.chunks.from.to(ref.file, nbr.chunks)  ## returns "chr:pos.start-chr:pos.end"
  
  ## pasting together the standard ssimp arguments and the from to positions.
  out <- paste0(ssimp.args, " --impute.range ", impute.range, "&")
  
  ## return
  write_tsv(data.frame(out), path = out.name, col_names = FALSE) ## writing out shell script	
  
  
}


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
url <-"https://drive.switch.ch/index.php/s/uOyjAtdvYjxxwZd/download"
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
