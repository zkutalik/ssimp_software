

f.svd <- function(C)
{
  tmp <- svd(round(C, dig = 40)) ## to avoid errors in fortran lapack code
  list(D = tmp$d, U = tmp$u)
}


extract.rho <- function(db, rows, cols)
{
  out <- cor(db[,rows], db[,cols], use = "pairwise.complete.obs")
  out <- matrix(out, ncol = length(cols), nrow = length(rows))
  colnames(out) <- cols
  rownames(out) <- rows
  out
  
}

extract.C <- function(db, rows)
{
  if (length(rows) > 1)
  {
    cor(db[,rows], use = "pairwise.complete.obs")
  }else{
    diag(length(rows))
  }
}




## including rho
tagging.speed <-function(betah = NULL, rho = NULL, C = NULL, lambda = NULL, tau = 1, to.imp = "eff")
{
  ## betah: vector with standardized effect sizes (N, eff needed)
  ## rho: 
  ## C: 
  ## eps: parameter for penalizing C
  ## tau: set to 1 if needed
  ## db.C for nrow sample size
  ## data.regular.offdiag: for when type regular off diag is set TRUE
  
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
  ret <- data.frame(upper = t(upper), lower = t(lower), lambda = lambda, N.mean = t(N.mean), N.max = t(N.max), N.min = t(N.min), r2 = t(r2))    
  
  
  return(ret)
  
}



f.impqual.r2 <- function(C, rho, target.snps, tag.snps, lambda = 1e-8)
{
  ## n = nrow(G.ref)
  
  C <- C * (1-lambda) + diag(nrow(C)) * lambda
  C.inv <- solve(C)
  r2 <-  diag(t(rho) %*% C.inv %*% rho)
  
  ## out
  out <- data.frame(SNP = target.snps, r2.classic = r2)
  
  return(out)
}


f.impqual.r2.penal <- function(C, rho, target.snps, tag.snps,lambda = 1e-8,  n = NULL)
{
  ## n = nrow(G.ref)
  
  out.pasaniuc <- f.impqual.r2(C, rho, target.snps, tag.snps, lambda = 1e-1)
  out <- f.impqual.r2(C, rho, target.snps, tag.snps, lambda)
  ## penal version
  m <- length(tag.snps)
  m.eff <- f.eff.number.tests(C)
  
  r2.penal.m <- out$r2.classic - m/n
  r2.penal.a <- out$r2.classic - m.eff/n
  r2.penal.b <- 1- (1 - out$r2.classic) * (n - 1)/(n - m.eff)
  
  ## out
  out <- data.frame(out, r2.classic.pasaniuc = out.pasaniuc$r2.classic, r2.penal.m = r2.penal.m, r2.penal.a = r2.penal.a, r2.penal.b =
                      r2.penal.b, n.tag.eff = m.eff, n.tag = m)
  
  
  #, r2.varsel =
  #out2$r2.varsel, r2.varsel.nvars = out2$r2.varsel.nvars)
  return(out)
}












tagging <- function(betah = NULL, rho = NULL, C = NULL, lambda = NULL, tau = 1,to.imp = "Z")
{
  ## betah: vector with standardized effect sizes
  ## rho: 
  ## C: 
  ## eps: parameter for penalizing C
  ## tau: set to 1 if needed
  ## db.C for nrow sample size
  ## data.regular.offdiag: for when type regular off diag is set TRUE
  
  n <- betah$N
  b <- betah[,to.imp]
  
  ## lambda
  C.lambda <- (1-lambda) * C + diag(ncol(C)) * lambda
  rho.lambda <- (1-lambda) * rho
  
  ## matrix
  C.inv <- solve(C.lambda)
  upper <- t(rho.lambda) %*% C.inv %*% (b)
  lower <- diag(t(rho.lambda) %*% C.inv %*% rho.lambda)
  
  ## impqual
  lambda.r2 <- 1e-6
  C.lambda <- (1-lambda.r2) * C + diag(ncol(C)) * lambda.r2
  rho.lambda <- (1-lambda.r2) * rho
  
  ## d
  r2 <- diag(t(rho.lambda) %*% C.inv %*% rho.lambda)
  
  
  N.mean <- NA
  N.max <- NA
  ret <- data.frame(upper = upper, lower = lower, lambda = lambda, N.mean = N.mean, N.max = N.max, r2 = r2)    
  
  return(ret)
  
}


f.impute.schur <- function(C, tag.snps.imp, betah, lambda = 0.01, what.to.imp  = "Z")
{
  ## columns and rows of C same than length and order of betah
  ## beta vector
  
  C.lambda <- (1-lambda) * C + diag(ncol(C)) * lambda
  C.inv <- solve(C.lambda)
  
  ## for r2
  lambda.r2 <- 1e-8
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
    # C =   svd.C, lambda = lambda, to.imp = what.to.imp)$upper
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
  
  OUT <- data.frame(SNP = tag.snps.imp, estim = estim, impqual = impqual)
  
  return(OUT)
  
}
## test for f.impute.schur
set.seed(3)
n <- 1000
p <- 4
G <- matrix(sample(0:2, n * p, replace = TRUE), ncol = p)
P <- G[,3] * 0.02 + rnorm(n)
Z <- sapply(1:ncol(G), function(x) summary(lm(formula("pheno ~ SNP"), data = data.frame(SNP = G[,x], pheno = P)))$coef[2,3])
C <- cor(G)
rsid <- paste0("rs",1:p)
colnames(C) <- rsid
rownames(C) <- rsid
betah <- data.frame(SNP = rsid, Z = Z)

f.impute.schur(C, rsid, betah)

f.gridsearch.20160613 <- function(ax = 1e-8, cx = 1, betah, C.tilde, n = NULL, miss = FALSE, snps.fake =
                                    NULL, mc.cores = 1, what.to.imp = "eff")
{
  
  ## if ax given, but cx null, then xseq is ax
  
  if(!is.null(cx))
  {
    xseq1 <- 10^(seq(log10(cx), log10(ax), length.out = 40))
    xseq2 <- seq(ax, cx, length.out = 40)
    xseq <- unique(c(0.1, xseq1, xseq2, 2/sqrt(n)))
  }else{
    xseq <- ax
  }
  
  ## for each lambda
  MSE <- NULL
  BETA.IMP <- NULL
  BETA.TRUE <- NULL
  
  l <- 1
  for (l in 1:length(xseq))
  {
    lambda2 <- xseq[l]
    C.tilde2 <- (1-lambda2) * C.tilde + diag(ncol(C.tilde)) * lambda2
    C.inv <- solve(C.tilde2)
    
    ## for each tSNP
    x <- 1
    out.mse <- mclapply(1:length(snps.fake), function(x)
    {
      
      i <- which(colnames(C.inv) == snps.fake[x])
      
      # if(rho)
      # {
      estim <- (-1)*C.inv[i,-i] %*% betah[-i, "Z"] /(C.inv[i,i])
      
      # }else{
      #estim <- (-1)* C.tilde[i,-i] %*% t(ginv(C.tilde2[i,-i])) %*% C.inv[i,-i] %*% t(t(betah[-i, what.to.imp]))/(C.inv[i,i])
      
      #  }
      
      
      #  svd.C <- f.svd(C.tilde[-i, -i])
      #  tmp <- tagging.speed(betah = betah[-i,], rho = as.matrix(C.tilde[-i,
      #  i], ncol = 1), C =   svd.C, lambda = lambda2, to.imp = what.to.imp)
      
      
      # estim <- f.v1(upper)
      
      beta.true <- betah[i, what.to.imp]
      mse <- f.mse(estim, beta.true) 
      
      
      return(c(mse, estim, beta.true))
    }, mc.cores = mc.cores)
    
    ## should be the same as
    # rho2 <- rho.tilde[-1,1]
    # C2 <- C.tilde[-1,-1]
    # svd.C <- f.svd(C2)
    # tagging.speed(betah = betah2[-1,], rho = as.matrix(rho2, ncol = 1), C = svd.C, lambda = lambda, miss = FALSE)
    
    MSE <- rbind(MSE, unlist(sapply(out.mse, function(x) x[1])))
    BETA.IMP <- rbind(BETA.IMP,  unlist(sapply(out.mse, function(x) x[2])))
    BETA.TRUE <- rbind(BETA.TRUE,  unlist(sapply(out.mse, function(x) x[3])))
  }
  
  is.na(MSE) <- which(is.infinite(MSE), arr.ind = TRUE)
  
  ## calc mean of each row (for each lambda)
  # mse.mean <- apply(MSE, 1, mean, na.rm = TRUE)
  mse.mean <- apply(MSE, 1, mean, na.rm = TRUE)
  ind <- which(mse.mean == min(mse.mean, na.rm = TRUE))
  lambda.out <- xseq[ind]
  
  
  # return(list(xseq, mse.mean))
  return(list(lambda.opt = lambda.out, lambdas = xseq, mse = mse.mean, MSE.mat = MSE, BETA.IMP =
                BETA.IMP, BETA.TRUE = BETA.TRUE, snps.fake
              = snps.fake))
  
}


select.tag.target.na <- function(C, tag.snps, target.snps = NULL)
{
  
  
  is.na(C) <- is.nan(C)
  
  if(!is.null(target.snps))
  {
    snps <- c(tag.snps, target.snps)
  }else{
    snps <- c(tag.snps)
  }
  
  
  ## exclude nas: tag.snps
  ind <- which(is.na(C[snps, snps]), arr.ind = TRUE)
  if(!is.null(nrow(ind)))
  {
    ind <-  matrix(ind[ind[,1] > ind[,2],], ncol = 2)
    tab <- sort(table(ind), decreasing = TRUE)
    while(nrow(ind) > 0)
    {
      excl <- as.numeric(names(tab[1]))
      snps <-snps[-excl]
      ind <- which(is.na(C[snps, snps]), arr.ind = TRUE)
      ind <-  matrix(ind[ind[,1] > ind[,2],], ncol = 2)
      tab <- sort(table(ind), decreasing = TRUE)
      
    }
    
  }
  
  tag.snps <- snps[snps %in% tag.snps]
  
  if(!is.null(target.snps))
    target.snps <- snps[snps %in% target.snps]
  
  length(snps)
  return(list(target.snps = target.snps, tag.snps = tag.snps))
}


select.tag.target.na.rho <- function(C, tag.snps, target.snps = NULL)
{
  
  
  is.na(C) <- is.nan(C)
  
  
  ## exclude nas: tag.snps
  ind <- which(is.na(C[tag.snps, target.snps]), arr.ind = TRUE)
  if(!is.null(nrow(ind)))
  {
    tab <- sort(table(ind[,2]), decreasing = TRUE)
    
    while(nrow(ind) > 0 & is.matrix(ind))
    {
      excl <- as.numeric(names(tab[1]))
      target.snps <- target.snps[-excl]
      ind <- which(is.na(C[tag.snps, target.snps]), arr.ind = TRUE)
      tab <- sort(table(ind[,2]), decreasing = TRUE)
      
    }
    
  }
  
  return(list(target.snps = target.snps, tag.snps = tag.snps))
}


#mat <- snps.na.num
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

