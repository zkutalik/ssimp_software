[//]: ========================================
# Software for Summary Statistics Imputation
[//]: ========================================

This command-line software enables summary statistics imputation (SSimp) for GWAS summary statistics.

## Installation
[//]: -------------------------------

### Source code
1. Clone the github folder. 
2. Run  `stu bin/ssimp`

### Compiled version
Download
* [ssimp 0.1 - Mac OS X]()
* [ssimp 0.1 - Ubuntu 12.04]()

## Documentation
[//]: -------------------------------
Run ssimp with no arguments to see the [usage message](https://github.com/sinarueeger/ssimp_software/blob/master/docu/usage.txt). 

Check-out [examples](https://github.com/sinarueeger/ssimp_software/blob/master/docu/examples.md)

We also provide a [detailed manual](https://github.com/sinarueeger/ssimp_software/blob/master/docu/usage.txt) that contains information not present in the usage text.

## Minimal example
[//]: -------------------------------
The minimal requirements are: 

(1) GWAS summary statistics stored in a text file with at least the following columns SNP-id (`MarkerName`), Z-statistic (`Z`), reference allele (`a1`) and risk allele (`a2`) and at least one row, and 

(2) the path to the reference panel with a `vcf` and `tbi`. 

`bin/ssimp --gwas data/my_gwas.txt` will impute the Z-statistics and generate a file `data/my_gwas.txt.ssimp.txt`.
