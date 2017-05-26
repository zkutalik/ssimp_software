# Manual of SSimp

This command line software tool enables summary statistics imputation (SSimp) for GWAS data. 
	
## Minimal example
The minimal requirements are: (1) GWAS summary statistics file with at least a SNP id, Z, reference alleles and risk alleles, and (2) the path to the reference panel. 

`../bin/ssimp --gwas ../data/gwas1.txt --ref ../ref/small.vcf.sample.vcf` will generate a `gwas1.imp.log` file and a `gwas1.imp.out` file.
	

## Arguments
[//]: <> (This is also a comment.)
Arguments can be shortend, e.g. `--wind` instead of `--window.width`

`--gwas [no default]`, path to GWAS dataset, extension does not matter, nor does text separation.  Columns should be named. Missings have to be marked as `NA`. quicktest, snptest, metal and plink output files will be automatically recognised. Minimal columns are SNP-id, Z, reference allele and risk allele. More info see below.

`--ref [no default]` path to vcf file (same folder should contain the tbi file)

`--out [gwasfilename]` string, if not define will be the gwasfilename (without extension ) + .imp

`--log [gwasfilename]` string, if not define will be the gwasfilename (without extension ) + .log, if set to FALSE, then no log file is stored

`--impute.range [no default]` should have the form of `chr:pos.start-chr:pos.end`, if just chr:NA, then the whole chr is imputed

`--impute.snp [no default]` strings

`--impute.snps [NULL]` txt file (e.g. snps2impute.txt) with a set of SNPs to impute (no header): rsid or chr:pos.hg19 (no quotes)

`--lambda [2/sqrt(n)]` numeric value or string (`2/sqrt(n)`, `optimize`), n are the number of individuals in the reference panel. `optimize` not yet implemented.

`--impute.maf [0]` numeric value, lower limit for variants to be imputed: everything above and equal this threshold will be imputed

`--tag.maf [0]` numeric value, lower limit for tag SNPs: everything above and equal this threshold will be used as tag SNPs. 

`--window.width [1e6]` numeric value, core window length

`--flanking.width [250e3]` numeric value, flanking space left and right side of the core window
		
`--missingness [TRUE]` enable variable sample size approach. This is automatically set to `FALSE` if `N` is not provided or `N` is set to `NA`.

### Note	
- If no `impute.range`, not `impute.snps` is chosen, then all variants in the refpanel are imputed that are not provided as tag SNPs.
- if your positions are not on hg19, you can use LiftOver as a command line tool: http://genome.ucsc.edu/cgi-bin/hgLiftOver
- The option `missingness` is automatically set to FALSE if `N` is not provided or `N` is set to `NA`.

## Format of GWAS dataset
Column names are automatically recognized. The log file will indicate which columns are used and recognized. Positions should be coded as the same as the reference panel. It is not needed, but N should be reported. 

## Imputation
For more details on the imputation, please see us (2017). 
`insert here a brief recap`

- if SNP ids are present in the GWAS, then pos copied from reference panel. 
- if SNP ids are not present, then the combination of Chr:Pos:A1:A2 are taken as identifier
- Because either SNP ids or Chr:Pos:A1:A2 are used as identifier, it is also possible to impute indels.
- Z statistics is imputed, along with the `N.imp` (an estimate for the sample size) and `r2.pred` (imputation quality)

 
## Output
### log file
The `.imp.log` file provides a summary of the imputation done and the results of sanity checks. 

### out file
The `.imp.out` file has
the following columns:

- `SNP` SNP id (if present)
- `Chr` Chromosome (only 1 to 22 right now)
- `Pos` Position (HG19)
- `Z.imp` Z
- `N.imp` estimation of N (only if missingness was set to `TRUE`)
- `r2.pred` imputation quality (adjusted as in ...)
- `A1` reference allele
- `A2` effect allele
- `window.nbr` in which window imputed
- `lambda` lambda used to impute

`.out.not` lists all the SNPs that were not found in the reference panel 
`.out.tnot` lists all the tSNPs that were not used as tSNPs (but needed): not found in the reference panel or had an NA in Z

if `Z.imp NA` and `r2.pred 0` means that there was not tag SNP.
