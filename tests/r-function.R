## sina, 2017-06-06

## these are the functions that sina used originally for her work
## and are now used to test aarons functions
## =================================================

library(data.table) # for 'rbindlist'


## top level function to impute
## roughly the same arguments as aarons function
## -----------------------------------------------
ssimp <- function(path.gwas, 
                  path.ref, 
                  path.outdir,
                  lambda = "sqrt", 
                  what.to.impute = "Z", 
                  gwas.names = c("SNP", "Chr", "Pos", "ref.allele", "effect.allele", "b", "SE", "Z", "N"),
                  tag.snps = NULL,
                  target.snps = NULL,
                  impute.tags = FALSE
                  )
{
  ## !!! here only SNP id names, not Chr:pos (except if SNP is defined as Chr:Pos)
  ## ----------------------------------------------------
  ## lambda: "sqrt" or numeric
  ## what.to.impute: string, column name to impute
  ## path.gwas: string to a tab separated file
  ## path.vcf: string to a tab separated file
  ## ref allele:
  ## effect allele: 
  ## gwas.names: naming of columns that refer to "SNP", "Chr", "Pos", "ref.allele", "effect.allele", "b", "SE", "Z", "N"
  ## if not present, fill with NA
  ## tag.snps, target.snps: character vector with predefined snp names
  ## impute.tags: TRUE or FALSE, if TRUE, tagging.schur will run on 100 or less randomly picked tag snps
  ## ----------------------------------------------------
  
  
  ## load gwas dataset
  ## -------------------
  if(length(grep("csv", path.gwas))> 0)
  {
    dat <- read.table(path.gwas, sep = ",", header = TRUE)
  }else{
    library(readr)
    dat <- read_tsv(path.gwas)
  }

  # Two steps for applying 'gwas.names':
  # First, put the columns in the appropriate order:
  cbind.data.frame(
  lapply(gwas.names, function(col.name) {
      if(is.na(col.name)) NA else dat[,col.name]
  })) -> dat
  # ... Second step: override the column names
  names(dat) <- c("SNP", "Chr", "Pos", "ref.allele", "effect.allele", "b", "SE", "Z", "N")

  ## only for giant
  if(all(!(names(dat) %in% "Z")) & what.to.impute == "Z")
  {
    dat$Z <- dat$b/dat$SE
  }

  ## load vcf
  ## ------
  library(readr)
  ref <- read_tsv(path.ref,comment = "##")
  G.raw <- ref[,10:ncol(ref)]
  
  ## build summary "sm" object
  ## ---------------------------
  sm <- ref[,1:5] ## sm is the summary dataset with columns SNP, chr, pos, REF and ALT
  sm$SNP <- sm$ID
  
  ## if tag.snps and target snps already defined, only transform these
  ## -------------
  if(!is.null(tag.snps) & !is.null(target.snps))
  {
 
      ind <- which(sm$SNP %in% c(tag.snps, target.snps))
      G.raw <- G.raw[ind, ]  
      sm <- sm[ind, ]  
  
  }
  
  ## take only subset (for quick use)
  ## -------------------
  #ind <- unique(c(1:100, which(sm$SNP %in% dat$SNP)))
  #G.raw <- G.raw[ind,]
  #sm <- sm[ind,]
  
  ## build "G" matrix
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
#  save(G, sm, file = "tmp_G_sm.RData")
  
  
  ## remove monomorphic from refpanel
  ## remove any SNP that has NAs from refpanel
  ## -------------------------------
  v.g <- apply(G, 2, var)
  missing.g <- apply(is.na(G), 2, sum)
  rem <- which(v.g == 0 | missing.g > 0)
  if (length(rem) > 0)
  {
    G <- G[,-rem]
    sm <- sm[-rem, ]
  }

  ## swap alleles
  ## --------------
  dat.swap <- merge(dat, sm, by = "SNP")
  rownames(dat.swap) <- dat.swap$SNP
  
  ## swap alleles if REF != Allele1 and ALT != Allele2
  ## swap the ones that have the alleles swapped
  ind.swap <- which(dat.swap$REF != dat.swap$ref.allele & dat.swap$ALT != dat.swap$effect.allele)
  ## keep the ones that have the alleles exactly right  
  ind.keep <- which(dat.swap$REF == dat.swap$ref.allele & dat.swap$ALT == dat.swap$effect.allele)
  ## meaning - drop the ones that have only one allele right
  
  if(length(ind.swap) > 0) dat.swap[ind,what.to.impute] <- dat.swap[ind,what.to.impute] * (-1)
  
  dat.swap <- dat.swap[c(ind.swap, ind.keep),]
  
  ## define tag snps and target snps
  ## -------------------------------
  if(is.null(tag.snps)) 
  {
    tag.snps <- intersect(as.character(dat.swap$SNP), sm$SNP)
  }else{
    tag.snps <- intersect(tag.snps, sm$SNP)
  }
  if(is.null(target.snps))
  {
    target.snps <- setdiff(sm$SNP, tag.snps)
  }else{
    target.snps <- intersect(target.snps, sm$SNP)
  }
    
  ## estimate C and rho
  ## --------------------
  C <- extract.C(G, tag.snps)
  rho <- extract.rho(G, tag.snps, target.snps)
  
  ## remove NAs
  ## ---------------
  
  
  ## imputation
  ## -------------
  if (lambda == "sqrt")
  {
    lambda <-  2/sqrt(nrow(G))
  }

  
  ## impute
  ## ----------
  svd.C <- f.svd(C)
  estim <- tagging.speed(betah = dat.swap[tag.snps, ], rho = rho, C = svd.C, lambda = lambda, tau = 1, to.imp = what.to.impute)
  
  ## same as 
  ## estim <- tagging(betah = dat.swap[tag.snps,], rho = rho, C = C, lambda = 2/sqrt(nrow(G)), tau = 1, to.imp = what.to.impute)
  
  r2 <- f.impqual.r2.penal(C, rho, target.snps, tag.snps, lambda  = 1e-8, n = nrow(G))
  #f.impqual.r2(C, rho, target.snps, tag.snps, lambda  = 1e-8)
  
  out <- data.frame(SNP = target.snps, Z.imp = estim$imp, impqual = r2$r2.penal.b)
  write_tsv(out, path = paste0(path.outdir, "imputations.txt"))
  
  ## print log file
  ## ---------------
  sink(paste0(path.outdir,"log.txt"))

  
  cat("GWAS:\n")
  print(path.gwas)
  
  cat("\nreference panel used:\n")
  print(path.ref)

  cat("\ncolumn name in GWAS to impute:\n")
  print(what.to.impute)
  
  cat("\ntag.snps:\n")
  print(tag.snps)
  
  cat("\nsnps2impute:\n")
  print(target.snps)

  pp(lambda)
  
  cat("\nC matrix:\n")
  print(C)
  
  cat("\nc matrix:\n")
  print(rho)

  sink()
  
  ## impute tag snps
  ## ------------------
  if(impute.tags)
  {
    set.seed(100)
    tag.snps.fake <- sample(tag.snps, min(c(100, length(tag.snps))))
    estim.tag <- f.impute.schur(C, tag.snps.fake, betah = dat.swap[tag.snps, ], lambda = lambda, what.to.imp  = what.to.impute)
  
    # plot(dat.swap[tag.snps.fake, what.to.impute], estim.tag$imp)
    #abline(a = 0, b=1)
    out.tag <- data.frame(estim.tag, Z.real = dat.swap[tag.snps.fake, what.to.impute])
    write_tsv(out.tag, path = paste0(path.outdir, "imputations_tag.txt"))
  }
  
}
  
  




## --------------------
f.svd <- function(C)
{
  tmp <- svd(round(C, dig = 40)) ## to avoid errors in fortran lapack code
  list(D = tmp$d, U = tmp$u)
}


## --------------------
extract.rho <- function(db, rows, cols)
{
  out <- cor(db[,rows], db[,cols], use = "pairwise.complete.obs")
  out <- matrix(out, ncol = length(cols), nrow = length(rows))
  colnames(out) <- cols
  rownames(out) <- rows
  out
  
}

## --------------------
extract.C <- function(db, rows)
{
  if (length(rows) > 1)
  {
    cor(db[,rows], use = "pairwise.complete.obs")
  }else{
    diag(length(rows))
  }
}




## --------------------
##
tagging.speed <-function(betah = NULL, rho = NULL, C = NULL, lambda = NULL, tau = 1, to.imp = "eff")
{
  ## betah: gwas data with same order of tag snps than C matrix 
  ## C: matrix with tag snps as rows and columns
  ## rho: matrix with tag snps as rows and target snps as columns
  ## lambda: for C matrix
  ## to.imp: column name in betah for ss to impute
  ## tau: obsolete
  
  svd.C <- C
  
  n <- betah$N
  b <- betah[,to.imp] # b.st
  
  b.wide <- matrix(b, nrow = length(b), ncol = ncol(rho), byrow = FALSE)
  
  
  rho.lambda <- (1-lambda) * rho
  
  ## d
  d <- ((1-lambda) * svd.C$D + lambda)^(-1)
  U <- svd.C$U
  
  ## (rho * inv(C) * beta)/(rho * inv(C) * rho)
  ## upper
  ## lower
  uu <- t(U) %*% rho.lambda ## rho
  vu <- t(U) %*% b.wide ## beta
  upper <- t(d) %*% (uu * vu)
  lower <- t(d) %*% (uu * uu)
  
  
  N.mean <- mean(n)
  N.max <- max(n)
  N.min <- min(n)
  
  
  ## impqual
  ## --------
  lambda.r2 <- 1e-6
  rho.lambda <- (1-lambda.r2) * rho
  
  ## d
  d <- ((1-lambda.r2) * svd.C$D + lambda.r2)^(-1)
  U <- svd.C$U
  uu <- t(U) %*% rho.lambda ## rho
  r2 <- t(d) %*% (uu * uu)
  
  # t(sqrt(alpha.var))
  ret <- data.frame(imp = t(upper), lambda = lambda)    
  
  
  return(ret)
  
}


## --------------------
f.impqual.r2 <- function(C, rho, target.snps, tag.snps, lambda = 1e-8)
{
  
  ## C: matrix with tag snps as rows and columns
  ## rho: matrix with tag snps as rows and target snps as columns
  ## targets.snps, tag.snps: respective char vectors
  ## lambda: for C matrix
  rho <- rho * (1-lambda) 
  C <- C * (1-lambda) + diag(nrow(C)) * lambda

  C.inv <- solve(C)
  r2 <-  diag(t(rho) %*% C.inv %*% rho)
  
  ## out
  out <- data.frame(SNP = target.snps, r2.classic = r2)
  
  return(out)
}

## --------------------
f.impqual.r2.penal <- function(C, rho, target.snps, tag.snps,lambda = 1e-8,  n = NULL)
{
  ## n = nrow(G.ref)
  
  out.pasaniuc <- f.impqual.r2(C, rho, target.snps, tag.snps, lambda = 1e-1)
  out <- f.impqual.r2(C, rho, target.snps, tag.snps, lambda)
  ## penal version
  m <- length(tag.snps)
  m.eff <- f.eff.number.tests(C, TRUE)
  
  r2.penal.m <- out$r2.classic - m/n
  r2.penal.a <- out$r2.classic - m.eff/n
  r2.penal.b <- 1- (1 - out$r2.classic) * (n - 1)/(n - m.eff - 1)
  
  ## out
  out <- data.frame(out, r2.classic.pasaniuc = out.pasaniuc$r2.classic, r2.penal.m = r2.penal.m, r2.penal.a = r2.penal.a, r2.penal.b =
                      r2.penal.b, n.tag.eff = m.eff, n.tag = m)
  
  return(out)
}





## slow imputation function
## -----------------------------
tagging <- function(betah = NULL, rho = NULL, C = NULL, lambda = NULL, tau = 1,to.imp = "Z")
{
  ## betah: gwas data with same order of tag snps than C matrix 
  ## C: matrix with tag snps as rows and columns
  ## rho: matrix with tag snps as rows and target snps as columns
  ## lambda: for C matrix
  ## to.imp: column name in betah for ss to impute
  ## tau: obsolete
  
  b <- betah[,to.imp]
  
  ## lambda
  C.lambda <- (1-lambda) * C + diag(ncol(C)) * lambda
  rho.lambda <- (1-lambda) * rho
  
  ## matrix
  C.inv <- solve(C.lambda)
  upper <- t(rho.lambda) %*% C.inv %*% (b)

  ## impqual >> done in a seperate function
  ## lambda.r2 <- 1e-6  
  ## C.lambda <- (1-lambda.r2) * C + diag(ncol(C)) * lambda.r2
  ## rho.lambda <- (1-lambda.r2) * rho
  ## r2 <- diag(t(rho.lambda) %*% C.inv %*% rho.lambda)
  
  ret <- data.frame(imp = upper, lambda = lambda)    
  
  return(ret)
  
}

## imputation of tag snps using schur complements
## -----------------------------------------------

f.impute.schur <- function(C, tag.snps.imp, betah, lambda = 0.01, lambda.r2 = 1e-8, what.to.imp  = "Z")
{
  ## C: matrix with tag snps as rows and columns
  ## tag.snps.imp: tag snps to impute
  ## betah: gwas data with same order of tag snps than C matrix 
  ## lambda: for C matrix
  ## labmda.r2: lambda used for imputation quality
  ## what.to.imp: column name in betah for ss to impute
  
  ## penalize lambda
  C.lambda <- (1-lambda) * C + diag(ncol(C)) * lambda
  C.inv <- solve(C.lambda)
  
  ## for r2
  C.lambda.r2 <- (1-lambda.r2) * C + diag(ncol(C)) * lambda.r2
  C.inv.r2 <- solve(C.lambda.r2)
  
  ## for each tSNP
  x <- 1
  
  estim <- sapply(1:length(tag.snps.imp), function(x)
  {
    
    i <- which(colnames(C.inv) == tag.snps.imp[x])
    
    ## imputation
    ## ------------
    estim <- (-1)*C.inv[i,-i] %*% betah[-i, what.to.imp] /(C.inv[i,i])
    ## identical to
    # C.lambda[i,-i] %*% stolve(C.lambda[-i,-i]) %*% betah[-i, what.to.imp]
    ## identical to
    
    # svd.C <- f.svd(C[-i, -i])
    # tagging.speed(betah = betah[-i,], rho = as.matrix(C[-i, i], ncol = 1),
    # C =   svd.C, lambda = lambda, to.imp = what.to.imp)$imp
    return(estim)
    
  })
  
  impqual <- sapply(1:length(tag.snps.imp), function(x)
  {
    
    i <- which(colnames(C.inv) == tag.snps.imp[x])
    
    ## imputation quality
    ## ------------------
    impqual <- (-1)*C.inv.r2[i,-i] %*% C.lambda.r2[-i,i] /(C.inv.r2[i,i])
    #   f.impqual.r2(C = C[-i, -i], rho = C[-i,i], x, colnames(C.inv)[-i])
    
    return(impqual)
    
  })
  
  OUT <- data.frame(SNP = tag.snps.imp, imp = estim, lambda = lambda, r2.classic = impqual)
  
  return(OUT)
  
}


## effective number of tests
## >> needed to calc r2.penal
## ---------------------------
f.eff.number.tests <- function(mat, cor.true = FALSE)
{
  if(!cor.true)
  {
    cor.mat <- cor(mat, use = 'pairwise.complete.obs')
  }else{
    cor.mat <- mat
  }
  cor.mat[which(is.na(cor.mat))] <- 0
  svd.data <- svd(cor.mat)
  sum.tp         <- 0
  zhc.correction <- 0
  while (sum.tp/sum(svd.data$d) < 0.995)
  {
    zhc.correction     <- zhc.correction+1
    sum.tp <- sum.tp+svd.data$d[zhc.correction]
  }
  return(zhc.correction)
}


# ## not used right now
# ## --------------------
# select.tag.target.na <- function(C, tag.snps, target.snps = NULL)
# {
#   
#   
#   is.na(C) <- is.nan(C)
#   
#   if(!is.null(target.snps))
#   {
#     snps <- c(tag.snps, target.snps)
#   }else{
#     snps <- c(tag.snps)
#   }
#   
#   
#   ## exclude nas: tag.snps
#   ind <- which(is.na(C[snps, snps]), arr.ind = TRUE)
#   if(!is.null(nrow(ind)))
#   {
#     ind <-  matrix(ind[ind[,1] > ind[,2],], ncol = 2)
#     tab <- sort(table(ind), decreasing = TRUE)
#     while(nrow(ind) > 0)
#     {
#       excl <- as.numeric(names(tab[1]))
#       snps <-snps[-excl]
#       ind <- which(is.na(C[snps, snps]), arr.ind = TRUE)
#       ind <-  matrix(ind[ind[,1] > ind[,2],], ncol = 2)
#       tab <- sort(table(ind), decreasing = TRUE)
#       
#     }
#     
#   }
#   
#   tag.snps <- snps[snps %in% tag.snps]
#   
#   if(!is.null(target.snps))
#     target.snps <- snps[snps %in% target.snps]
#   
#   length(snps)
#   return(list(target.snps = target.snps, tag.snps = tag.snps))
# }
# 
# 
# select.tag.target.na.rho <- function(C, tag.snps, target.snps = NULL)
# {
#   
#   
#   is.na(C) <- is.nan(C)
#   
#   
#   ## exclude nas: tag.snps
#   ind <- which(is.na(C[tag.snps, target.snps]), arr.ind = TRUE)
#   if(!is.null(nrow(ind)))
#   {
#     tab <- sort(table(ind[,2]), decreasing = TRUE)
#     
#     while(nrow(ind) > 0 & is.matrix(ind))
#     {
#       excl <- as.numeric(names(tab[1]))
#       target.snps <- target.snps[-excl]
#       ind <- which(is.na(C[tag.snps, target.snps]), arr.ind = TRUE)
#       tab <- sort(table(ind[,2]), decreasing = TRUE)
#       
#     }
#     
#   }
#   
#   return(list(target.snps = target.snps, tag.snps = tag.snps))
# }
# 




pp <- function(...) {
	sl = substitute(list(...))
	N=length(sl)
	pf=parent.frame()

	old.options.width=options()$width # store the old width
	options(width=old.options.width-6) # subtract 6 to make space for the six spaces

	rbindlist( lapply(2:N, function(n) {
		paste0(deparse(sl[[n]],width.cutoff=500),collapse='') -> nm
		#my.deparse(sl[[n]]) -> nm
		utils::capture.output( eval (sl[[n]]
									, envir=pf
									) ) -> val
		l = length(val)
		if(l>1) {
			val = paste0('      ', val) # six spaces
		}
		list(nm=nm,str=paste0(collapse='\n',val),l=l)
	})) -> nms.and.vals

	options(width=old.options.width) # restore the old width

	nms.and.vals[  l==1 & substr(str,1,4) == "[1] "
				 , str := substr(str,5,100000)     ]
	nms.and.vals[l==1, if(length(nm)>0) {
		cat(sep=''
			,'{ '
			, paste0(collapse=',', nm)
			, ' }\t|=|  '
			)
		catn( paste0(collapse='\t|  ', str) )
		NULL
	}]
	nms.and.vals[, {
        if(l!=1) {
            catn(sep='', '{{ ', nm, ' }}:' )
            catn( str )
        }
    NULL} ,by=1:nrow(nms.and.vals)]
	invisible(NULL)
}
catn <- function(...) {
	cat(...)
	cat('\n')
}
