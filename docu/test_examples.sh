#!/bin/bash

compiled/ssimp-osx-0.3 --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt

compiled/ssimp-osx-0.3 --gwas gwas/small.random.txt.gz --ref ref/small.vcf.sample.vcf.gz --out output.txt


compiled/ssimp-osx-0.3 --gwas gwas/small.random.p.b.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt

compiled/ssimp-osx-0.3 --gwas gwas/small.random.chr.pos.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt


compiled/ssimp-osx-0.3 --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt --impute.range 22:18000000-22:18100075

compiled/ssimp-osx-0.3 --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt --impute.range 22


compiled/ssimp-osx-0.3 --gwas gwas/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz --ref ref/sub1KG-tiny/chr{CHRM}.vcf.gz  --out output.txt --impute.range 20-22

compiled/ssimp-osx-0.3 --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt --impute.snp gwas/listofimputesnps.txt

compiled/ssimp-osx-0.3 --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt --impute.snp <(echo rs2587101 rs2277831 | tr ' ' '\n')


compiled/ssimp-osx-0.3 gwas/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz ref/sub1KG-tiny/chr{CHRM}.vcf.gz output.txt --sample.names ref/link.to.1kg.data/integrated_call_samples_v3.20130502.ALL.panel/sample/super_pop=EUR --impute.snps <(echo rs148911000 rs111659000 rs183059100 rs76979500 rs150095300 rs115012100 rs187649300 rs560286600 rs78808100 | tr ' ' '\n')

compiled/ssimp-osx-0.3 --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt --tag.snp gwas/listoftagsnps.txt

compiled/ssimp-osx-0.3 --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt --tag.snp <(echo rs2305001 rs10854521 | tr ' ' '\n')


compiled/ssimp-osx-0.3 --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt --sample.names ref/filename.samples.small.txt

compiled/ssimp-osx-0.3 --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt --sample.names ref/filename.samples.pop.txt/sample/super_pop=EAS

compiled/ssimp-osx-0.3 --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt --lambda 0.01

compiled/ssimp-osx-0.3 --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt --impute.maf 0.05

compiled/ssimp-osx-0.3 --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt --tag.maf 0.05

compiled/ssimp-osx-0.3 --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt --window.width 500000 --flanking.width 100000

compiled/ssimp-osx-0.3 --gwas gwas/small.random.n.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt --missingness dep

compiled/ssimp-osx-0.3 --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt --log my_ssimp_logfile






compiled/ssimp-osx-0.3 --gwas gwas/small.random.x.txt --ref ref/sub1KG-tiny/chrX.vcf.gz --out output.txt

compiled/ssimp-osx-0.3 --gwas gwas/small.random.x.txt --ref ref/sub1KG-tiny/chrX.vcf.gz --out output.txt --impute.snp <(echo rs183055800 rs146115300 rs150092800 | tr ' ' '\n')

compiled/ssimp-osx-0.3 --gwas gwas/small.random.x.txt --ref ref/sub1KG-tiny/chrX.vcf.gz --out output.txt --impute.snp <(echo X:18008535 X:18010084 X:18037863 | tr ' ' '\n')

compiled/ssimp-osx-0.3 --gwas gwas/small.random.x.txt --ref ref/sub1KG-tiny/chrX.vcf.gz --out output.txt --impute.snp <(echo 23:18008535 23:18010084 23:18037863 | tr ' ' '\n')

compiled/ssimp-osx-0.3 --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt --impute.snp <(echo X:18008535 X:18010084 X:18037863 | tr ' ' '\n')


compiled/ssimp-osx-0.3 --gwas gwas/small.random.txt --ref 1KG/EUR --out output.txt

compiled/ssimp-osx-0.3 --gwas gwas/small.random.txt ~/reference_panels/1KG/ALL.chr{CHRM}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz output.txt --sample.names ~/reference_panels/1KG/integrated_call_samples_v3.20130502.ALL.panel/sample/super_pop=EUR
