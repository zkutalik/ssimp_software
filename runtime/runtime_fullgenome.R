list(rm = ls())

## chromosome to impute
chr <- 22

## record start time
## ---------------
start.time <- Sys.time()

## print log file
sink(paste0("log_start_chr",chr,".txt"))
cat("starting time:\n")
print(start.time)
sink()


## run ssimp (small example)
## -----------
## system("../bin/ssimp --gwas ../gwas/small.random.csv --impute.snp impute.snp --tags.snp tags.snp --wind 250000 --flan 250000 --lambda 0.01 --ref ../ref/small.vcf.sample.vcf.gz --out runtime_test.txt")

## bigger example, but still small reference panel >> not working 
## system("../bin/ssimp --gwas wood.txt --ref ../ref/sub1KG-tiny/chr{CHRM}.vcf.gz --impute.snp impute.snp2 --wind 250000 --flan 250000")

## another bigger examples, with default reference panel
## system("../bin/ssimp --gwas wood.txt --impute.snp impute.snp2 --wind 250000 --flan 250000")

## "wrong" headers:: check if original headers (without Z, swapped a1 and a2) work
##system("../bin/ssimp --gwas /data/sgg/sina/public.data/meta.summaries/meta.giant/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt --impute.snp impute.snp2 --wind 250000 --flan 250000 --out test.ssimp.txt")

## impute wood with 1KG, but only chr 22
system(paste0("../bin/ssimp --gwas /data/sgg/sina/public.data/meta.summaries/meta.giant/wood.txt --impute.range ",chr))

## record end time
## ---------------
end.time <- Sys.time()

## print log file
sink(paste0("log_chr",chr,".txt"))
cat("starting:\n")
print(start.time)
cat("ending:\n")
print(end.time)
sink()
