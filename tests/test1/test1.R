
setwd("/Users/admin/Documents/Work/Projects/tagging/ssimp_software/ref/")
library(readr)

## ref
vcf.raw <- read_tsv("1000genomes_90snps_100samples.vcf")
snps <- vcf.raw$ID

## prep vcf
map <- vcf.raw[,1:5]
G.char <- t(as.data.frame(vcf.raw[,-(1:9)]))
G <- matrix(NA, ncol = ncol(G.char), nrow = nrow(G.char))
G[G.char == "0|0"] <- 0
G[G.char == "0|1"] <- 1
G[G.char == "1|0"] <- 1
G[G.char == "1|1"] <- 2
colnames(G) <- map$ID

## gwas
path.gwas <- "/Users/admin/Documents/Work/Projects/tagging/ssimp_software/gwas/"
gwas <- read_tsv(paste0(path.gwas, "gwas_example.txt"))

## snps2impute
snps2imp <- setdiff(snps, gwas$SNP)
tsnps <- intersect(snps, gwas$SNP)

## imputation


## prep C.mat
v.g <- apply(G, 2, var)
C.mat <- cor(G[,v.g > 0], use = "pairwise")
lambda <- 2/sqrt(nrow(G))

C.mat <- (C.mat * (1-lambda)) + diag(C.mat) * lambda

tsnps <- intersect(tsnps, colnames(C.mat))
snps2imp <- intersect(snps2imp, colnames(C.mat))

rho <- C.mat[tsnps, snps2imp]
C <- C.mat[tsnps, tsnps]

## gwas data
gwas.sub <- gwas[gwas$SNP %in% tsnps, ]
Z <- gwas.sub$BETA/gwas.sub$SE


## actual imputaton
res <- t(rho) %*% solve(C) %*% Z
out <- data.frame(SNP = rownames(res), Z.imp = res)

write_tsv(out, "test1.R.out")
