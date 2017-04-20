
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

## adapt the gwas's
## --------------------
path.raw <- "/Users/admin/Documents/Work/Projects/tagging/ssimp_software/gwas/gwas.original/"
path.gwas <- "/Users/admin/Documents/Work/Projects/tagging/ssimp_software/gwas/"

tmp <- read_tsv(paste0(path.raw, "gwas_example.txt"))
set.seed(3)
ind <- sample(1:nrow(tmp), n.tsnps)
tmp[ind,"CHR"] <- chr
tmp[ind,"SNP"] <- tsnps
tmp[ind,"POS"] <- pos

## store out
## -----------
write_tsv(tmp, paste0(path.gwas, "gwas_example.txt"))