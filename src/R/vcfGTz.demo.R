source('src/R/vcfGTz.reader.R')

filename='ref/compressed/50mb.vcfGTz'

num.SNPs = vcfGTz_num_SNPs(filename)
PP(num.SNPs)

