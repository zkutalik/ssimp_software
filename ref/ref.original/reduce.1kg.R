
library(readr)
vcf.raw <- read_table("ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",  n_max=10 ,skip = 200)

head(vcf.raw)

## ungz
gzip -d ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

## remove header
grep -v "##" ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf > temp

## remove snps (rows)
wc -l temp
1103548
region starting points: 1e5, 4e5, 7e5

tail -n +100000 temp | head -n 30 - > temp1
tail -n +400000 temp | head -n 30 - > temp2
tail -n +700000 temp | head -n 30 - > temp3
tail -n +1 temp | head -n 1 - > temp.header

## paste it together
temp1 <- read_tsv("temp1", col_names = FALSE)
temp2 <- read_tsv("temp2", col_names = FALSE)
temp3 <- read_tsv("temp3", col_names = FALSE)
temp.header <- read_tsv("temp.header", col_names = FALSE)

temp.all <- rbind(temp1, temp2, temp3)
names(temp.all) <- unlist(temp.header)
dim(temp.all) ## should be 90

## extract only some individuals (columns)
ped <- read_tsv("integrated_call_samples_v3.20130502.ALL.panel")
set.seed(3)
id <- sample(ped$sample[ped$super_pop == "EUR"], 100)


temp.subset <- temp.all[,c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", 
"FORMAT", id)]
dim(temp.subset)

## store out
write_tsv(temp.subset, "1000genomes_90snps_100samples.vcf")
