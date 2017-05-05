source('src/R/vcfGTz.reader.R')

args = commandArgs(trailing=T)
filename=args[1]

num.SNPs = vcfGTz_num_SNPs(filename)
PP(num.SNPs)

file_offsets = vcfGTz_get_internal_offsets(filename, 0:3 %c% -1:3 %c% num.SNPs %c% (num.SNPs+1))
PP(file_offsets)

for(which.field in 1:7) {
    x = vcfGTz_get_1field_from_internal_offsets(filename, file_offsets, which.field)
    PP(x)
}

m = vcfGTz_get_012calls_from_internal_offsets(filename, file_offsets)
pp(m[                               ,1:15])
pp(m['rs62224609'                   ,1:15])
pp(m[(!is.na(rownames(m))) & 'rs62224609' ==rownames(m)   ,1:15])
pp(m[(!is.na(rownames(m))) & '.'          ==rownames(m)   ,1:15])
