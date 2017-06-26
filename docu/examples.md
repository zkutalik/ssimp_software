[//]: ==================================
# Examples
[//]: ==================================


	


## Minimal example
[//]: -------------------------------
The minimal requirements are: (1) GWAS summary statistics stored in a text file with at least the following columns SNP-id (`MarkerName`), Z-statistic (`Z`), reference allele (`a1`) and risk allele (`a2`) and at least one row, and (2) the path to the reference panel. 

`bin/ssimp --gwas data/my_gwas.txt --ref ref/my_reference_panel.vcf` will generate a file `data/my_gwas.txt.ssimp.txt`, containing the imputation results.

<sup>If P-values are provided instead Z-statistics, there needs to be an extra column containing the effect sizes (the P-value will be turned into a Z-statistic, and therefore needs a negative sign if the effect size is negative). </sup>

## Change filename of imputed results
[//]: -------------------------------

`bin/ssimp --gwas data/my_gwas.txt --ref ref/my_reference_panel.vcf --out pathtoresults/filename`

## Store a log file
[//]: -------------------------------

`bin/ssimp --gwas data/my_gwas.txt --ref ref/my_reference_panel.vcf --log pathtolog/filename`


## use 16 cores
[//]: -------------------------------
`bin/ssimp --gwas data/my_gwas.txt --ref ref/my_reference_panel.vcf --n.cores 16`


## Point to your own reference panel
[//]: -------------------------------
`bin/ssimp --gwas data/my_gwas.txt --ref ref/my_reference_panel.vcf`

`ref` should contain the vcf and tbi files. 


## Select individuals from reference panel
[//]: -------------------------------
`bin/ssimp --gwas data/my_gwas.txt --ref ref/my_reference_panel.vcf --sample.names ref/filename.samples.txt`

filename.samples.txt contains sample id's separated by new lines (no header). 

`bin/ssimp --gwas data/my_gwas.txt --ref ref/my_reference_panel.vcf --sample.names ref/filename.samples.txt/sample/super_pop=EUR`

in this case `filename.samples.txt` contains sample id's (`sample`) along with a second attribute (here `super_pop`) that has different values, among them is `EUR`, for which we separate. 


## impute range of SNPs on a specific chromosome
[//]: -------------------------------
`bin/ssimp --gwas data/my_gwas.txt --ref ref/my_reference_panel.vcf --impute.range 1:1000-1:300000`

## impute chromosome 1
[//]: -------------------------------
`bin/ssimp --gwas data/my_gwas.txt --ref ref/my_reference_panel.vcf --impute.range 1`

## impute chromosome 1 to 22
[//]: -------------------------------
`bin/ssimp --gwas data/my_gwas.txt --ref ref/my_reference_panel.vcf --impute.range 1-22`

## impute specific set of SNPs
[//]: -------------------------------
`bin/ssimp --gwas data/my_gwas.txt --ref ref/my_reference_panel.vcf --impute.snp listofimputesnps.txt`

`listofimputesnps.txt` contains SNP id's separated by new lines (no header).

If it is only a handful of SNPs it might be easier to use:

`bin/ssimp --gwas data/my_gwas.txt --ref ref/my_reference_panel.vcf --impute.snp <(echo rs1 rs2 | tr ' ' '\n')`

## use a set of tag SNPs
[//]: -------------------------------
`bin/ssimp --gwas data/my_gwas.txt --ref ref/my_reference_panel.vcf --tag.snp listoftagsnps.txt`

`listoftagsnps.txt` contains SNP id's separated by new lines (no header).

If it is only a handful of SNPs it might be easier to use:

`bin/ssimp --gwas data/my_gwas.txt --ref ref/my_reference_panel.vcf --tag.snp <(echo rs3 rs4 rs5 | tr ' ' '\n')`


## shrinking by 0.01
[//]: -------------------------------
`bin/ssimp --gwas data/my_gwas.txt --ref ref/my_reference_panel.vcf --lambda 0.01`

## impute SNPs with MAF > 0.05
[//]: -------------------------------
`bin/ssimp --gwas data/my_gwas.txt --ref ref/my_reference_panel.vcf --impute.maf 0.05`

## use tag SNPs with MAF > 0.05
[//]: -------------------------------
`bin/ssimp --gwas data/my_gwas.txt --ref ref/my_reference_panel.vcf --tag.maf 0.05`


## use a smaller core window and flanking region
[//]: -------------------------------
`bin/ssimp --gwas data/my_gwas.txt --ref ref/my_reference_panel.vcf --window.width 500000 --flanking.width 100000`


## use the dependent missingness approach
[//]: -------------------------------
`bin/ssimp --gwas data/my_gwas.txt --ref ref/my_reference_panel.vcf --missingness dep`





























	