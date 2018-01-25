[//]: ==================================
# Examples
[//]: ==================================

All examples work with test data from the `ref` and `gwas` folder, therefore use much smaller GWAS datasets and reference panels than in reality to limit computation time. 

### Minimal example to impute Z-statistics
[//]: -------------------------------
The minimal requirements are:
1. `--gwas`: path to the *GWAS summary statistics* text file, containing at least the following columns SNP-id, Z-statistic, reference allele and risk allele and at least one row, 
2. `--ref`: path to the *reference panel*,
3. `--out`: path for output file.

`ssimp --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt`

For more info regarding automatically column names recognition of the GWAS file, see section `GWAS dataset` in the [manual](https://github.com/sinarueeger/ssimp_software/blob/master/docu/manual.md).


### --gwas dataset is gzipped
[//]: -------------------------------

`ssimp` **`--gwas gwas/small.random.txt.gz`** ` --ref ref/small.vcf.sample.vcf.gz --out output.txt`

`.gz` works, `.zip` does not work.

### Downloading the reference panel
[//]: -------------------------------

For detailed instructions and explanations see detailed instructions in [usage-text](https://github.com/sinarueeger/ssimp_software/blob/master/docu/usage.txt).

A **quick solution** is to run `ssimp` without reference panel, but with a shortcut indicating 1KG and a preferred population. This will create a folder called `refpanel` and download 1KG (all populations, not only the selected one).

`ssimp --gwas gwas/small.random.txt` **`--ref 1KG/EUR`** ` --out output.txt`


### You do not have the Z-statistics available, but you have...
[//]: -------------------------------

#### ... P-values and betas
No special argument needed, but GWAS input file needs to contain effect size `b` along with the P-value `p`. 

The header of the GWAS file should be: `SNP  a1  a2  b   p`

`ssimp` **`--gwas gwas/small.random.p.b.txt`** ` --ref ref/small.vcf.sample.vcf.gz --out output.txt`


#### ... odds ratios
No special argument needed, but GWAS input file needs to contain Z-statistics of the odds ratios. 

Alternatively, if no Z-statistic is available, provide `b=log(OR)` and the P-value.

The header of the GWAS file should either be `SNP  a1  a2  Z` or `SNP  a1  a2  b   p`.


### Impute b, and not Z
[//]: -------------------------------
First, impute Z-statistics as shown above, then transform the Z-statistic into `se(b)` and `b` using this [R-function](https://github.com/sinarueeger/ssimp_software/blob/master/transform_z_to_b.R).


### Chr and Pos instead of SNP-identifier
[//]: -------------------------------
Provide chromosome and position instead, but include a SNP column that is empty too.

`ssimp` **`--gwas gwas/small.random.chr.pos.txt`** ` --ref ref/small.vcf.sample.vcf.gz --out output.txt`


### Impute a bp range on a specific chromosome
[//]: -------------------------------

`ssimp --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--impute.range 22:18000000-22:18100075`**


### Impute chromosome 22
[//]: -------------------------------

`ssimp --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--impute.range 22`**


### Impute chromosome 20 to 22
[//]: -------------------------------

`ssimp --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--impute.range 20-22`**


### Impute specific set of SNPs
[//]: -------------------------------

`ssimp --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--impute.snp gwas/listofimputesnps.txt`**

`listofimputesnps.txt` contains SNP id's separated by new lines (no header).

If it is only a few of SNPs it might be easier to use:

`ssimp --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--impute.snp <(echo rs2587101 rs2277831 | tr ' ' '\n')`**

A more complex example with SNPs from two different chromosomes:

`ssimp gwas/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz ref/sub1KG-tiny/chr{CHRM}.vcf.gz output.txt --sample.names ref/link.to.1kg.data/integrated_call_samples_v3.20130502.ALL.panel/sample/super_pop=EUR --impute.snps <(echo rs148911000 rs111659000 rs183059100 rs76979500 rs150095300 rs115012100 rs187649300 rs560286600
rs78808100 | tr ' ' '\n')`



### Imputing SNPs on the X-chromosome
[//]: -------------------------------

This was fixed in version 0.3. 

`ssimp --gwas gwas/small.random.x.txt --ref ref/sub1KG-tiny/chrX.vcf.gz --out output.txt`

Impute single SNPs (rsid)
`ssimp --gwas gwas/small.random.x.txt --ref ref/sub1KG-tiny/chrX.vcf.gz --out output.txt --impute.snp <(echo rs183055800 rs146115300 rs150092800 | tr ' ' '\n')`

Impute range of SNPs (X:pos)
`ssimp --gwas gwas/small.random.x.txt --ref ref/sub1KG-tiny/chrX.vcf.gz --out output.txt --impute.range 23:18000000-23:18100075`


### Use a set of tag SNPs
[//]: -------------------------------

`ssimp --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--tag.snp gwas/listoftagsnps.txt`**

`listoftagsnps.txt` contains SNP id's separated by new lines (no header).

If it is only a few of SNPs it might be easier to use:

`ssimp --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--tag.snp <(echo  rs2305001 rs10854521 | tr ' ' '\n')`**


### Select individuals from reference panel (short version)
[//]: -------------------------------

#### 1)
`ssimp --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--sample.names ref/filename.samples.small.txt`**

`filename.samples.small.txt` contains sample id's separated by new lines (no header). 

#### 2)
`ssimp --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--sample.names ref/filename.samples.pop.txt/sample/super_pop=EAS`**
Here, `filename.samples.txt` contains sample id's (`sample`) along with a second attribute (here `super_pop`) that has different values, among them is `EUR`, for which we separate. 


### Shrinking LD matrix by 0.01
[//]: -------------------------------

`ssimp --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--lambda 0.01`**


### Impute SNPs with MAF > 0.05
[//]: -------------------------------

`ssimp --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--impute.maf 0.05`**


### Use tag SNPs with MAF > 0.05
[//]: -------------------------------

`ssimp --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--tag.maf 0.05`**


### Use a smaller core window and flanking region
[//]: -------------------------------

`ssimp --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--window.width 500000 --flanking.width 100000`**


### Use the dependent missingness approach
[//]: -------------------------------

`ssimp` **`--gwas gwas/small.random.n.txt`** `--ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--missingness dep`**

Your GWAS summary statistics needs to have a sample size column. 

### Store a copy of all console output to a file
[//]: -------------------------------

`ssimp --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt` **`--log my_ssimp_logfile`**


### Sanitycheck with reimputed tag SNPs
[//]: -------------------------------
Check R-file `docu/sanitycheck_reimputed.R`


### Parallelize computation (TBD)
[//]: -------------------------------


### Use your own reference panel
[//]: -------------------------------
Follow the instructions in [docu/usage](https://github.com/sinarueeger/ssimp_software/blob/master/docu/usage.txt).
