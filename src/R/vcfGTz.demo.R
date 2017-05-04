source('src/R/vcfGTz.reader.R')

filename='ref/compressed/50mb.vcfGTz'

num.SNPs = vcfGTz_num_SNPs(filename)
PP(num.SNPs)

file_offsets = vcfGTz_get_internal_offsets(filename, 0:3 %c% -1:3)
PP(file_offsets)

vcfGTz_get_CHROM_from_internal_offsets(filename, file_offsets)

