
rm(list= ls())

## checkout the data quickly
## -----------
library(readr)
vcf.raw <- read_table("ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf",  n_max=10 ,skip = 200)

head(vcf.raw)

## un-zip
## -----------
system("gzip -d ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")

## remove header
## -------------
system("grep -v '##' ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf > temp")

## remove snps (rows)
## ---------------------
system("wc -l temp")

##1103548
##region starting points: 16877135
## -------------
system("tail -n +400000 temp | head -n 100000 - > temp1")
#tail -n +400000 temp | head -n 30 - > temp2
#tail -n +700000 temp | head -n 30 - > temp3
system("tail -n +1 temp | head -n 1 - > temp.header")

## paste it together
## -------------------
temp1 <- read_tsv("temp1", col_names = FALSE)
#temp2 <- read_tsv("temp2", col_names = FALSE)
#temp3 <- read_tsv("temp3", col_names = FALSE)
temp.header <- read_tsv("temp.header", col_names = FALSE)

temp.all <- temp1
names(temp.all) <- unlist(temp.header)
dim(temp.all) ## should be 1e5

## extract only some individuals (columns)
ped <- read_tsv("integrated_call_samples_v3.20130502.ALL.panel")
set.seed(3)
id <- ped$sample[ped$super_pop == "EUR"]


temp.subset <- temp.all[,c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", 
                           "FORMAT", id)]
dim(temp.subset)


## take only a subsample according to giatn wood chr 22
## ---------------

ind <- which(temp.subset$POS >= 30132167 & temp.subset$POS <= 31632167)

out <- temp.subset[ind, ]


## store out
## -------------
#write_tsv(temp.subset, "1000genomes_90snps_100samples.vcf")
write_tsv(out, "1000genomes_somechr22snps_EURsamples.vcf")


## library(fastVCF)

## tmp <- vcf_range_of_positions(root = "/data/sgg/aaron/shared/fastVCF/1KG", chr = 1, first.pos = 30, last.pos = 1e5)       


