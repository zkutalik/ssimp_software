[//]: ==================================
# Examples
[//]: ==================================

All examples work with test data from the `ref` and `gwas` folder, therefore use significantly smaller GWAS dataset and reference panels to limit computation time. 

### Minimal example (default: imputation of Z-statistics)
[//]: -------------------------------
The minimal requirements are:
1. `--gwas`: path to the *GWAS summary statistics* text file, containing at least the following columns SNP-id, Z-statistic, reference allele and risk allele and at least one row, 
2. `--ref`: path to the *reference panel*,
3. `--out`: path for output file.

`ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt`

For more info regarding automatically column names recognition of the GWAS file, see section `GWAS dataset` in [detailed manual](https://github.com/sinarueeger/ssimp_software/blob/master/docu/manual.md).


### Download the reference panel
[//]: -------------------------------

For detailed instructions and explanations see `docu/usage.txt`. 

A **quick solution** is to run `ssimp` without reference panel, but with a shortcut indicating 1KG and a preferred population. This will create a folder called `refpanel` and download 1KG.

`ssimp --gwas gwas/small.random.csv --ref 1KG/EUR --out output.txt`


### You don not have the Z-statistics available, but you have...
[//]: -------------------------------

#### ... P-values and betas
No special argument needed, but GWAS input file needs to contain effect size `b` along with the P-value `p`. 

So the header of the GWAS file would be: `SNP  a1  a2  b   p`

`ssimp --gwas gwas/small.random.p.b.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt`


#### ... odds ratios
No special argument needed, but GWAS input file needs to contain Z-statistics of the odds ratios. 

Alternatively, if no Z-statistics is available, provide `b=log(OR)` and the P-value.

So the header of the GWAS file would either be: `SNP  a1  a2  Z`

Or: `SNP  a1  a2  b   p`


### Chr and Pos instead of SNP as SNP-identifier
[//]: -------------------------------
Provide chromosome and position instead.

`ssimp --gwas gwas/small.random.chr.pos.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt`


### Impute a bp range on a specific chromosome
[//]: -------------------------------

`ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--impute.range 22:16000000-22:16050075`**


### Impute chromosome 22
[//]: -------------------------------

`ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--impute.range 22`**


### Impute chromosome 20 to 22
[//]: -------------------------------

`ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--impute.range 20-22`**


### Impute specific set of SNPs
[//]: -------------------------------

`ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--impute.snp listofimputesnps.txt`**

`listofimputesnps.txt` contains SNP id's separated by new lines (no header).

If it is only a handful of SNPs it might be easier to use:

`ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--impute.snp <(echo rs587755077 rs587697622 | tr ' ' '\n')`**


### Use a set of tag SNPs
[//]: -------------------------------

`ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--tag.snp listoftagsnps.txt`**

`listoftagsnps.txt` contains SNP id's separated by new lines (no header).

If it is only a handful of SNPs it might be easier to use:

`ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--tag.snp <(echo rs62224618 rs333 | tr ' ' '\n')`**


### Select individuals from reference panel (short version)
[//]: -------------------------------

`ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--sample.names ref/filename.samples.small.txt`**

filename.samples.txt contains sample id's separated by new lines (no header). 

`ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--sample.names ref/filename.samples.txt/sample/super_pop=EUR`**

in this case `filename.samples.txt` contains sample id's (`sample`) along with a second attribute (here `super_pop`) that has different values, among them is `EUR`, for which we separate. 


### Shrinking by 2/sqrt(n)
[//]: -------------------------------

`ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--lambda "2/sqrt(n)"`**


### Impute SNPs with MAF > 0.05
[//]: -------------------------------

`ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--impute.maf 0.05`**


### Use tag SNPs with MAF > 0.05
[//]: -------------------------------

`ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--tag.maf 0.05`**


### Use a smaller core window and flanking region
[//]: -------------------------------

`ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--window.width 500000 --flanking.width 100000`**


### Use the dependent missingness approach
[//]: -------------------------------

`ssimp --gwas gwas/small.random.n.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--missingness dep`**


### Store a copy of all console output to a file
[//]: -------------------------------

`ssimp --gwas gwas/small.random.csv --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--log my_ssimp_logfile`**


### Paralellize computation (TBD)
[//]: -------------------------------


### Use your own reference panel
[//]: -------------------------------
Follow the instructions in docu/usage.txt
