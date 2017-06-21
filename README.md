[//]: ========================================
# Software for Summary Statistics Imputation
[//]: ========================================

This repositorys' purpose is to facilitate the software developement of the SSimp and its documentation. 

The software should be included in quicktest and therefore written in C++.

## Documentation
[//]: -------------------------------
The folder `docu` contains the user manual.

## Minimal example
[//]: -------------------------------
The minimal requirements are: 

(1) GWAS summary statistics stored in a text file with at least the following columns SNP-id (`MarkerName`), Z-statistic (`Z`), reference allele (`a1`) and risk allele (`a2`) and at least one row, and 

(2) the path to the reference panel with a `vcf` and `tbi`. 

`bin/ssimp --gwas data/my_gwas.txt` will impute the Z-statistics and generate a file `data/my_gwas.txt.ssimp.txt`.

## Test data
[//]: -------------------------------
The folder `gwas` contains examples of summary statistics datasets that are used for test purposes.

## Reference panel
[//]: -------------------------------
The folder `ref` contains reference panels for test purposes. 
