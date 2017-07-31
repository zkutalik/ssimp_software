[//]: ==================================
# Examples
[//]: ==================================


### Minimal example
[//]: -------------------------------
The minimal requirements are *GWAS summary statistics* stored in a text file with at least the following columns SNP-id, Z-statistic, reference allele and risk allele and at least one row.

`ssimp my_gwas.txt` will generate a file `my_gwas.txt.ssimp` containing the imputation results.

For more info regarding automatical column names recognition in `my_gwas.txt`, see section `GWAS dataset` in [detailed manual](https://github.com/sinarueeger/ssimp_software/blob/master/docu/manual.md).

### Paralellize computation: use 10 cores 
```diff 
- (this option is not yet integrated)
```
`ssimp my_gwas.txt` **`--n.cores 10`**


### Point to your own reference panel

Follow the instructions on `ssimp --help`.

### Select individuals from reference panel (short version)

`ssimp my_gwas.txt ~/ref_panels/my_reference_panel.vcf output.txt` **`--sample.names ref/filename.samples.txt`**

filename.samples.txt contains sample id's separated by new lines (no header). 

`ssimp my_gwas.txt ~/ref_panels/my_reference_panel.vcf output.txt` **`--sample.names ref/filename.samples.txt/sample/super_pop=EUR`**

in this case `filename.samples.txt` contains sample id's (`sample`) along with a second attribute (here `super_pop`) that has different values, among them is `EUR`, for which we separate. 


### Impute range of SNPs on a specific chromosome
`ssimp my_gwas.txt ~/ref_panels/my_reference_panel.vcf output.txt` **`--impute.range 11:1000-11:300000`**


### Impute chromosome 6
`ssimp my_gwas.txt ~/ref_panels/my_reference_panel.vcf output.txt` **`--impute.range 6`**


### Impute chromosome 4 to 12
`ssimp my_gwas.txt ~/ref_panels/my_reference_panel.vcf output.txt` **`--impute.range 4-12`**


### Impute specific set of SNPs
`ssimp my_gwas.txt ~/ref_panels/my_reference_panel.vcf output.txt` **`--impute.snp listofimputesnps.txt`**

`listofimputesnps.txt` contains SNP id's separated by new lines (no header).

If it is only a handful of SNPs it might be easier to use:

`ssimp my_gwas.txt ~/ref_panels/my_reference_panel.vcf output.txt` **`--impute.snp <(echo rs1 rs2 | tr ' ' '\n')`**


### Use a set of tag SNPs
`ssimp my_gwas.txt ~/ref_panels/my_reference_panel.vcf output.txt` **`--tag.snp listoftagsnps.txt`**

`listoftagsnps.txt` contains SNP id's separated by new lines (no header).

If it is only a handful of SNPs it might be easier to use:

`ssimp my_gwas.txt ~/ref_panels/my_reference_panel.vcf output.txt` **`--tag.snp <(echo rs3 rs4 rs5 | tr ' ' '\n')`**


### Shrinking by 0.01
`ssimp my_gwas.txt ~/ref_panels/my_reference_panel.vcf output.txt` **`--lambda 0.01`**


### Impute SNPs with MAF > 0.05
`ssimp my_gwas.txt ~/ref_panels/my_reference_panel.vcf output.txt` **`--impute.maf 0.05`**


### Use tag SNPs with MAF > 0.05
`ssimp my_gwas.txt ~/ref_panels/my_reference_panel.vcf output.txt` **`--tag.maf 0.05`**


### Use a smaller core window and flanking region
`ssimp my_gwas.txt ~/ref_panels/my_reference_panel.vcf output.txt` **`--window.width 500000 --flanking.width 100000`**


### Use the dependent missingness approach
`ssimp my_gwas.txt ~/ref_panels/my_reference_panel.vcf output.txt` **`--missingness dep`**




### Store a copy of all console output to a file

`ssimp my_gwas.txt ~/ref_panels/my_reference_panel.vcf output.txt` **`--log my_ssimp_logfile`**


























	
