#!/bin/bash

compiled/ssimp-osx-0.2 --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output1.txt


compiled/ssimp-osx-0.2 --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output2.txt --impute.range 22:18000000-22:18100075

compiled/ssimp-osx-0.2 --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output3.txt --impute.range 22


compiled/ssimp-osx-0.2 --gwas gwas/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz ref/sub1KG-tiny/chr{CHRM}.vcf.gz  --out output4.txt --impute.range 20-22

compiled/ssimp-osx-0.2 --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output5.txt --impute.snp gwas/listofimputesnps.txt

compiled/ssimp-osx-0.2 --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output6.txt --impute.snp <(echo rs2587101 rs2277831 | tr ' ' '\n')


compiled/ssimp-osx-0.2 gwas/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz ref/sub1KG-tiny/chr{CHRM}.vcf.gz output7.txt --sample.names ref/link.to.1kg.data/integrated_call_samples_v3.20130502.ALL.panel/sample/super_pop=EUR --impute.snps <(echo rs148911000 rs111659000 rs183059100 rs76979500 rs150095300 rs115012100 rs187649300 rs560286600 rs78808100 | tr ' ' '\n')

compiled/ssimp-osx-0.2 --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output8.txt --tag.snp gwas/listoftagsnps.txt

compiled/ssimp-osx-0.2 --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output9.txt --tag.snp <(echo rs2305001 rs10854521 | tr ' ' '\n')


compiled/ssimp-osx-0.2 --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output10.txt --sample.names ref/filename.samples.small.txt

compiled/ssimp-osx-0.2 --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output11.txt —sample.names ref/filename.samples.txt/sample/super_pop=EUR

compiled/ssimp-osx-0.2 --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output12.txt --lambda 0.01

compiled/ssimp-osx-0.2 --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output13.txt --impute.maf 0.05

compiled/ssimp-osx-0.2 --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output14.txt --tag.maf 0.05

compiled/ssimp-osx-0.2 --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output15.txt --window.width 500000 --flanking.width 100000

compiled/ssimp-osx-0.2 --gwas gwas/small.random.n.txt --ref ref/small.vcf.sample.vcf.gz --out output16.txt --missingness dep

compiled/ssimp-osx-0.2 --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output17.txt --log my_ssimp_logfile
