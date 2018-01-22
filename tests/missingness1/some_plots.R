cts = readr::cols(r2.pred = 'd', r2_reimputed='d', N.imp = 'd')
none    = readDT('none/imputations.txt', '\t', T, col_types=cts)[source=='GWAS']
dep     = readDT('dep/imputations.txt' , '\t', T, col_types=cts)[source=='GWAS']
ind     = readDT('ind/imputations.txt' , '\t', T, col_types=cts)[source=='GWAS']

stopifnot(none[,.(chr,pos,Allele1,Allele2,z_imp)] == dep[,.(chr,pos,Allele1,Allele2,z_imp)])
stopifnot(none[,.(chr,pos,Allele1,Allele2,z_imp)] == ind[,.(chr,pos,Allele1,Allele2,z_imp)])

none2 = none[,.(chr,pos,z_obs = z_imp, zN_re = Z_reimputed, zN_re_r2 = r2_reimputed)]
 dep2 =  dep[,.(                       zD_re = Z_reimputed, zD_re_r2 = r2_reimputed)]
 ind2 =  ind[,.(                       zI_re = Z_reimputed, zI_re_r2 = r2_reimputed)]

pp(none2, dep2, ind2)
j = cbind(none2, dep2, ind2)
pp(j)

plot.with.cor = function(x,y) {
    plot(x,y, main=cor(x,y)
         ,xlab = deparse(substitute(x))
         ,ylab = deparse(substitute(y))
         )
}

par(mfrow=c(3,3))
j[,{
    plot(zN_re_r2, zI_re_r2); abline(0:1)
    plot(zN_re_r2, zD_re_r2); abline(0:1)
    plot(zD_re_r2, zI_re_r2); abline(0:1)

    plot.with.cor(z_obs, zN_re)
    plot.with.cor(z_obs, zD_re)
    plot.with.cor(z_obs, zI_re)

    plot(zN_re, zI_re); abline(0:1)
    plot(zN_re, zD_re); abline(0:1)
    plot(zD_re, zI_re); abline(0:1)
}]

