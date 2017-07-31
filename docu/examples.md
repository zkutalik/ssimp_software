[//]: ==================================
# Examples
[//]: ==================================

All exmples work with test data from the `ref` and `gwas` folder, therefore use significantly smaller GWAS dataset and reference panels to limit computation time. 

### Minimal example
[//]: -------------------------------
The minimal requirements are:
1. `--gwas`: path to the *GWAS summary statistics* text file, containing at least the following columns SNP-id, Z-statistic, reference allele and risk allele and at least one row, 
2. `--ref`: path to the *reference panel*,
3. `--out`: path for output file.

`ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt`

For more info regarding automatical column names recognition of the GWAS file, see section `GWAS dataset` in [detailed manual](https://github.com/sinarueeger/ssimp_software/blob/master/docu/manual.md).


### Impute a bp range on a specific chromosome
`ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--impute.range 22:16000000-22:16050075`**


### Impute chromosome 22
`ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--impute.range 22`**


### Impute chromosome 20 to 22
`ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--impute.range 20-22`**


### Impute specific set of SNPs
`ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--impute.snp listofimputesnps.txt`**

`listofimputesnps.txt` contains SNP id's separated by new lines (no header).

If it is only a handful of SNPs it might be easier to use:

`ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--impute.snp <(echo rs587755077 rs587697622 | tr ' ' '\n')`**


### Use a set of tag SNPs
`ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--tag.snp listoftagsnps.txt`**

`listoftagsnps.txt` contains SNP id's separated by new lines (no header).

If it is only a handful of SNPs it might be easier to use:

`ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--tag.snp <(echo rs62224618 rs333 | tr ' ' '\n')`**


### Select individuals from reference panel (short version)

`ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--sample.names ref/filename.samples.small.txt`**

filename.samples.txt contains sample id's separated by new lines (no header). 

`ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--sample.names ref/filename.samples.txt/sample/super_pop=EUR`**

in this case `filename.samples.txt` contains sample id's (`sample`) along with a second attribute (here `super_pop`) that has different values, among them is `EUR`, for which we separate. 


### Shrinking by 2/sqrt(n)
`ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--lambda "2/sqrt(n)"`**


### Impute SNPs with MAF > 0.05
`ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--impute.maf 0.05`**


### Use tag SNPs with MAF > 0.05
`ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--tag.maf 0.05`**


### Use a smaller core window and flanking region
`ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--window.width 500000 --flanking.width 100000`**


### Use the dependent missingness approach
`ssimp my_gwas.txt ~/ref_panels/my_reference_panel.vcf output.txt` **`--missingness dep`**


### Store a copy of all console output to a file
`ssimp my_gwas.txt ~/ref_panels/my_reference_panel.vcf output.txt` **`--log my_ssimp_logfile`**


### Paralellize computation: use 10 cores 
```diff 
- (this option is not yet integrated)
```
`ssimp my_gwas.txt` **`--n.cores 10`**


### Use your own reference panel
Follow the instructions on `ssimp --help`.
