cts = readr::cols(r2.pred = 'd', r2_reimputed='d', N.imp = 'd')
none    = readDT('none/imputations.txt', '\t', T, col_types=cts)[source=='SSIMP']
dep     = readDT('dep/imputations.txt' , '\t', T, col_types=cts)[source=='SSIMP']
ind     = readDT('ind/imputations.txt' , '\t', T, col_types=cts)[source=='SSIMP']

stopifnot(none[,.(chr,pos,Allele1,Allele2)] == dep[,.(chr,pos,Allele1,Allele2)])
stopifnot(none[,.(chr,pos,Allele1,Allele2)] == ind[,.(chr,pos,Allele1,Allele2)])

pp(none)

none2 = none[,.(chr,pos,zN_imp = z_imp)]
 dep2 =  dep[,.(chr,pos,zD_imp = z_imp)]
 ind2 =  ind[,.(chr,pos,zI_imp = z_imp)]



pp(none2, dep2, ind2)
j = cbind(none2, dep2, ind2)
pp(j)
plot.with.cor = function(x,y) {
    plot(x,y, main=cor(x,y)
         ,xlab = deparse(substitute(x))
         ,ylab = deparse(substitute(y))
         )
}

par(mfrow=c(1,3))
j[,{
    plot.with.cor(zN_imp, zD_imp); abline(0:1)
    plot.with.cor(zN_imp, zI_imp); abline(0:1)
    plot.with.cor(zI_imp, zD_imp); abline(0:1)
#
    #plot.with.cor(z_obs, zN_re)
    #plot.with.cor(z_obs, zD_re)
    #plot.with.cor(z_obs, zI_re)
#
    #plot(zN_re, zI_re); abline(0:1)
    #plot(zN_re, zD_re); abline(0:1)
    #plot(zD_re, zI_re); abline(0:1)
}]

