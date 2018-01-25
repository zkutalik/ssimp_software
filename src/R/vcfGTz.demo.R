source('src/R/vcfGTz.reader.R')

args = commandArgs(trailing=T)
filename=args[1]

num.SNPs = vcfGTz_num_SNPs(filename)
PP(num.SNPs)

file_offsets = vcfGTz_get_internal_offsets(filename, 0:3 %c% -1:3 %c% num.SNPs %c% (num.SNPs+1))
PP(length(file_offsets))

df = data.frame(lapply(1:7, function(which.field) {
    vcfGTz_get_1field_from_internal_offsets(filename, file_offsets, which.field)
}))
colnames(df) <- ''
pp(data.table(df))

m = vcfGTz_get_012calls_from_internal_offsets(filename, file_offsets)
pp(m[                               ,1:15])
