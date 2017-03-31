
setwd("/Users/admin/Documents/Work/Projects/tagging/ssimp_software/ref/")
library(readr)
vcf.raw <- read_tsv("1000genomes_90snps_100samples.vcf")

snps <- vcf.raw$ID[vcf]
n.snps <- length(snps)

n.tsnps <- 2/3*n.snps
set.seed(3)
tsnps <- sample(snps, n.tsnps)
ind <- which(snps %in% tsnps)
chr <- vcf.raw[ind,1]
pos <- vcf.raw[ind,2]
tsnps <- vcf.raw[ind,3]
