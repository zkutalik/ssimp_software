[//]: ==================================
# Examples
[//]: ==================================


### Minimal example
[//]: -------------------------------
The minimal requirements are: (1) GWAS summary statistics stored in a text file with at least the following columns SNP-id (`MarkerName`), Z-statistic (`Z`), reference allele (`a1`) and risk allele (`a2`) and at least one row, and (2) the path to the reference panel. 

`ssimp my_gwas.txt ~/ref_panels/my_reference_panel.vcf output.txt` will generate a file `output.txt` containing the imputation results.

<sup>If P-values are provided instead Z-statistics, there needs to be an extra column containing the effect sizes (the P-value will be turned into a Z-statistic, and therefore needs a negative sign if the effect size is negative). </sup>


### Store a log file

`ssimp my_gwas.txt ~/ref_panels/my_reference_panel.vcf output.txt` **`--log my_ssimp_logfile`**


### Use 10 cores

`ssimp my_gwas.txt ~/ref_panels/my_reference_panel.vcf output.txt` **`--n.cores 10`**


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





























	