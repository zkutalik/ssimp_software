testname = Sys.getenv('testname')
readDT(paste0('tests/test',testname, '.R/imputations.txt'),'\t',T)  -> fromR
readDT(paste0('tests/test',testname,   '/imputations.txt'),'\t',T) [,.(SNP=SNP,z.CXX=z_imp,iq.CXX=r2.pred)]    -> fromCXX
stopifnot(nrow(fromR) == nrow(fromCXX))
join = fromR[fromCXX,on='SNP']
z.cor  = join[,abs(Z.imp-z.CXX)]
iq.cor = join[,abs(impqual- iq.CXX)]
stopifnot(min( z.cor ) < 1e-5)
stopifnot(min(iq.cor ) < 1e-5)
catn('OK, both correlations are close to 1.0')
