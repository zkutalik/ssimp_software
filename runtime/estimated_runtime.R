#####################################################################################
## Author: Sina Rueeger [sina *.* rueger *a*t* unil *.* ch]
## Project: 
##        
## Time-stamp: <[estimated_runtime.R] by SR Thu 22/06/2017 10:40 (CEST)>
##
## Description:
## 
##
## History:
## 
#####################################################################################


## take the runtime from chr 22, # of snps imputed and scale it up to the number of snps in each chr
## --------------------------
libary(readr)
nams <- grep(".txt", grep("runtime_",dir(), value = TRUE), value = TRUE)
rem <- grep("tmp", nams)
nams <- nams[-rem]

out.raw <- lapply(nams, function(x) read_tsv(x))
out <- do.call("rbind", out.raw)
out$runtime <- as.numeric(out$end.time -out$start.time)

mat <-data.frame(chr = c(1:22), bp=c( 1625210,1734854,1452031,1467148, 1303335,1312948,1190098, 1123847, 878596,1011303,1014519,971857,731711, 672316, 609085,656612, 576924, 576257, 474673, 454010, 293847, 287072))

out <- merge(out, mat, by = "chr",all=TRUE)

out$estim.chr1 <- NA
out$estim.chr22 <- NA
out$estim.chr1[out$chr != 1] <- out$runtime[out$chr == 1]/out$bp[out$chr == 1] * out$bp[out$chr != 1]
out$estim.chr1[out$chr == 1] <- out$runtime[out$chr == 1]
out$estim.chr22[out$chr != 22] <- out$runtime[out$chr == 22]/out$bp[out$chr == 22] * out$bp[out$chr != 22]
out$estim.chr22[out$chr == 22] <- out$runtime[out$chr == 22]
sum(out$estim.chr1)
sum(out$estim.chr22)
 
