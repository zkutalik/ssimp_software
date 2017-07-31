[//]: ========================================
# Software for Summary Statistics Imputation
[//]: ========================================

This command-line software enables summary statistics imputation (SSimp) for GWAS summary statistics. 

The only input needed from the user are the **GWAS summary statistics**. The reference panel (needed for LD computation) is downloaded automatically during the installation of the software.

## Installation
[//]: -------------------------------

### Source code
1. Clone the github folder. 
2. Run  `stu bin/ssimp`

### Compiled version
Download
* [ssimp 0.1 - Mac OS X]()
* [ssimp 0.1 - Ubuntu 12.04]()

## Minimal example
[//]: -------------------------------

`bin/ssimp --gwas data/my_gwas.txt` will impute the Z-statistics and generate a file `data/my_gwas.txt.ssimp.txt`. `data/my_gwas.txt` contains at least the following columns: SNP-id, Z-statistic, reference allele and risk allele, and at least one row. 

my_gwas.txt contains the  with 



## Documentation
[//]: -------------------------------
Run ssimp with no arguments to see the [usage message](https://github.com/sinarueeger/ssimp_software/blob/master/docu/usage.txt). 

Check-out [examples](https://github.com/sinarueeger/ssimp_software/blob/master/docu/examples.md).

We also provide a [detailed manual](https://github.com/sinarueeger/ssimp_software/blob/master/docu/usage.txt) that contains information not present in the usage message.

