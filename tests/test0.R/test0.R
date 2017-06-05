## set paths
## -------------------
path.out <- "/Users/admin/Documents/Work/Projects/tagging/ssimp_software/tests/test0.R/"
path.ref <- "/Users/admin/Documents/Work/Projects/tagging/ssimp_software/ref/"
path.gwas <- "/Users/admin/Documents/Work/Projects/tagging/ssimp_software/gwas/"

setwd(path.out)

## load function
## -------------------
source("r-function.R")

## load gwas dataset
## -------------------
dat <- read.table(paste0(path.gwas, "small.UKBB.both.csv"), sep = ",", header = TRUE)
dat$SNP <- dat$rnpid
rownames(dat) <- dat$SNP
what.to.impute <- "z.from.peff"
## a1 effect allele

## region
## ---------------
chr <- 22
corew <- 1e6
sw <- 250e3
core.start <- 16500000
core.end <- core.start + corew
window.start <- core.start - sw
window.end <- core.end + sw


## load vcf
## ------
library(readr)
ref <- read_delim(paste0(path.ref, "small.vcf.sample.vcf"), delim = "\t", comment = "##")
G.raw <- ref[,10:ncol(ref)]
sm <- ref[,1:5]
sm$SNP <- sm$ID

## take only subset
## -------------------
#ind <- unique(c(1:100, which(sm$SNP %in% dat$SNP)))
#G.raw <- G.raw[ind,]
#sm <- sm[ind,]

## build G matrix
## ---------------
G <- matrix(NA, ncol = nrow(G.raw), nrow=ncol(G.raw))

f.translate <- function(x)
{
  ## turns
  ## 0|0 into 0 
  ## 1|0 into 1
  ## 0|1 into 1
  ## 1|1 into 2
  
	tmp <- unlist(strsplit(as.character(x), "|"))
	ret <- NA
	
	if(as.numeric(tmp[1])<= 1 & as.numeric(tmp[3]) <= 1)
	{
		ret <- as.numeric(tmp[1]) + as.numeric(tmp[3])	
	}
	
	return(ret)
	
}

## translate each snp into a dosage
library(parallel)
tmp <- mclapply(1:nrow(G.raw), function(i)
{
	 sapply(G.raw[i,], function(k) f.translate(k))
}, mc.cores = 4) ## each list entry is a SNP
G <- do.call(cbind.data.frame, tmp)
names(G) <- sm$SNP
stopifnot(ncol(G) == nrow(sm))
save(G, sm, file = "tmp_G_sm.RData")

## allele swaping
## --------------
dat.swap <- merge(dat, sm, by = "SNP")
rownames(dat.swap) <- dat.swap$SNP

ind.swap <- which(dat.swap$REF != dat.swap$a2 & dat.swap$ALT != dat.swap$a1)

if(length(ind.swap) > 0) dat.swap[ind,what.to.impute] <- dat.swap[ind,what.to.impute] * (-1)

## define tag snps and target snps
## -------------------------------
tag.snps <- intersect(as.character(dat.swap$SNP), sm$SNP)
target.snps <- c("rs587697622", "rs587755077","rs587654921")#setdiff(sm$SNP, tag.snps)


## estimate C and rho
## --------------------
C <- extract.C(G, tag.snps)
rho <- extract.rho(G, tag.snps, target.snps)

## remove NAs


## imputation
## -------------
svd.C <- f.svd(C)
estim <- tagging.speed(betah = dat.swap[tag.snps, ], rho = rho, C = svd.C, lambda = 2/sqrt(nrow(G)), tau = 1, to.imp = what.to.impute)

## same as 
## estim <- tagging(betah = dat.swap[tag.snps,], rho = rho, C = C, lambda = 2/sqrt(nrow(G)), tau = 1, to.imp = what.to.impute)

r2 <- f.impqual.r2.penal(C, rho, target.snps, tag.snps, lambda  = 1e-8, n = nrow(G))
#f.impqual.r2(C, rho, target.snps, tag.snps, lambda  = 1e-8)

out <- data.frame(SNP = target.snps, Z.imp = estim$upper, impqual = r2$r2.penal.b)

write_delim(out, path = paste0(path.out, "imputation.txt"), delim = "\t")


## compare with aarons results
## -----------------------------
out.aaron <- read_delim("/Users/admin/Documents/Work/Projects/tagging/ssimp_software/tests/test0/imputations.txt", skip = 1, n_max = 3, delim = "\t", col_names = FALSE)
out.aaron
out


