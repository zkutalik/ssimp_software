wood %just.once% { readDT('../../../gwas/too.big/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz','\t',T)[,.(SNP=MarkerName, N)] %>%setkey(SNP) }

read.imps <- function(fn) { readDT  (fn, '\t', T
                                     ,  col_types = "iinccccnnnnn"
                                     ,  na='-inf'
                                     )[source=='GWAS'] }
none    %just.once% read.imps('chrALL.MISSnone.TAGimputations.txt')
dep     %just.once% read.imps('chrALL.MISSdep.TAGimputations.txt')
ind     %just.once% read.imps('chrALL.MISSind.TAGimputations.txt')
pp(ind, dep, none)

stopifnot(none[,.(chr,pos,Allele1,Allele2,z_imp)] == dep[,.(chr,pos,Allele1,Allele2,z_imp)])
stopifnot(none[,.(chr,pos,Allele1,Allele2,z_imp)] == ind[,.(chr,pos,Allele1,Allele2,z_imp)])

none2 = none[,.(chr,pos,SNP,z_obs = z_imp,
                                       zN_re = Z_reimputed, zN_re_r2 = r2_reimputed)]
 dep2 =  dep[,.(                       zD_re = Z_reimputed, zD_re_r2 = r2_reimputed)]
 ind2 =  ind[,.(                       zI_re = Z_reimputed, zI_re_r2 = r2_reimputed)]

pp(none2, dep2, ind2)
j = cbind(none2, dep2, ind2)

j = j[complete.cases(j)] # 2'111'509 out of 2'121'899
j = na.fail(wood[j,on='SNP'])



pp(j)

pp(j[,.( minmax( zN_re_r2 )
       , minmax( zD_re_r2 )
       , minmax( zI_re_r2 )
     )])

pp(j[order(abs(zN_re_r2))])
pp(j[order(abs(zD_re_r2))])
pp(j[order(abs(zI_re_r2))])

if(1){
    par(mfrow=c(1,3))
    j[,{
        plot(N, zN_re_r2);
        plot(N, zD_re_r2);
        plot(N, zI_re_r2);
    }]
}

plot.with.cor = function(x,y) {
    plot(x,y, main=cor(x,y)
         ,xlab = deparse(substitute(x))
         ,ylab = deparse(substitute(y))
         )
}

if(0){
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
}
