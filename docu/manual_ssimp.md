[//]: ==================================
# Manual of SSimp
[//]: ==================================

This command line software enables summary statistics imputation (SSimp) for GWAS summary statistics. 

<sup>Details to the method can be found in *us (2017)*. (+ *Pasaniuc et al. (2014)*, or more recently in a review on the use of summary statistics by *Pasaniuc & Price (2016)*). </sup>

	


## Minimal example
[//]: -------------------------------
The minimal requirements are: (1) GWAS summary statistics stored in a text file with at least the following columns SNP-id, Z-statistic, reference allele and risk allele and at least one row, and (2) the path to the reference panel. 

<sup>All other arguments have defaults defined (see below).</sup>

`bin/ssimp --gwas data/my_gwas.txt --ref ref/my_reference_panel.vcf` will generate a file `imputation.txt`, containing the imputation, and a log file called `output.actual`.

<sup>If P-values are provided instead Z-statistics, there needs to be an extra column containing the effect sizes (the P-value will be turned into a Z-statistics, and therefore needs a negative sign if the effect size is negative too). </sup>
	

## Arguments
[//]: -------------------------------
Here are all arguments listed, with a defined default in brackets, and some formal... . Note that arguments can be shortend, e.g. `--wind` instead of `--window.width`

`--gwas [no default]`, path to GWAS dataset, in quotes, extension (e.g. `.txt`) does not matter, nor does text separator (e.g. `\t`).  Columns need to be named after common conventions (see file `../header_translation.md`). Missings have to be marked as `NA`. QUICKTEST, SNPTEST, METAL AND PLINK output files will be automatically recognised. The minimal set of columns that should be provided, are: SNP-id, Z-statistics, reference allele and risk allele. For more info on possible sets of columns, see `Note` below.

`--ref [no default]` path to vcf file, in quotes (same folder should contain the `tbi` file).

`--out [gwasfilename]` string, in quotes. If not define will be the gwas filename (without extension ) + `.imp`. Path and name of the imputation results.

`--log [gwasfilename]` string, in quotes. If not define will be the gwas filename (without extension ) + `.log`. If set to `FALSE`, then no log file is produced. Path and name of the log file.

`--impute.range [no default]` text. Should have the form of `chrCHR:pos.start-chrCHR:pos.end`. If `chrCHR`, then the whole chromosome `CHR` is imputed. For `chr1-chr5` chromosome 1 to chromosome 5 are imputed. 

`--impute.snp [no default]` text. E.g. `(rs1, rs2)`. List of tag SNPs. 

`--impute.snps [NULL]` text. To define the SNPs to impute there are two options: (1) provide a text file (e.g. `snps2impute.txt`) with a set of SNP names (no header), give a list of SNP names, e.g. `(rs3, rs4)`. SNP names can be either rsid's or `chr:pos.hg19` (no quotes). 

`--lambda [2/sqrt(n)]` numeric value or string (`2/sqrt(n)`, `optimize`), n are the number of individuals in the reference panel. `optimize` not yet implemented. Lambda controls the shrinking of the correlation matrix (lambda = 0 applies no shrinking, lambda = 1 turns the correlation matrix into the identity matrix).

`--impute.maf [0]` numeric value. Lower MAF limit for variants to be imputed: everything above and equal this threshold will be imputed.

`--tag.maf [0]` numeric value. Lower MAF limit for tag SNPs: everything above and equal this threshold will be used as tag SNPs. 

`--window.width [1e6]` numeric value. Core window length.

`--flanking.width [250e3]` numeric value. Flanking space left and right side of the core window.
		
`--missingness [TRUE]` logical. Enables variable sample size approach. This is automatically set to `FALSE` if `N` is not provided or `N` is set to `NA`.

### Note	
[//]: -------
- If `impute.range` and `impute.snps` are not defined, then all variants in the refpanel are imputed that are not provided as tag SNPs.
- If the GWAS positions do not match the reference panel positionn, use use LiftOver as a command line tool: http://genome.ucsc.edu/cgi-bin/hgLiftOver
- The option `missingness` is automatically set to FALSE if `N` is not provided or `N` is set to `NA`.
- The minimal columns required are `SNP`, `A1`, `A2`, `Z`. If `Z` is not present, but `P` and `b` are, `Z` is calculated through `P` and `b`. Alternatively, if `b` and `SE` are present, then it is also possible to calculate `Z` via `b` and `SE`. 
- Odds ratios need to be provided as Z-statistics or log-transformed into effect sizes.

## Format of GWAS dataset
[//]: -------------------------------
Column names are automatically recognized using commonly used names. The `log` file will indicate which columns are used and recognized. Positions should match the positions in the reference panel (e.g. both hg19). It is recommended to provide the sample size (N), as incorporating missingness leads to a more accurate estimate. 

## Reference panel
[//]: -------------------------------


## Imputation
[//]: -------------------------------
For more details on the summary statistics imputation method, please see our-paper (2017). 

`insert here a brief recap of the method`

- If SNP-ID are present in the GWAS, then pos copied from reference panel. 
- If SNP-ID are not present, then the combination of Chr:Pos:A1:A2 are taken as identifier
- Because either SNP-ID or Chr:Pos:A1:A2 are used as identifier, it is also possible to impute indels.
- Z statistics are imputed, along with `N.imp` (an estimate for the sample size) and `r2.pred` (adjusted imputation quality).
- To speed up computation, we use a sliding window approach (`--window.width` and `--flanking.width`). SNPs to be imputed are assigned to one window. The window number is reported in the `.imp.out` file. 

 
## Output
[//]: -------------------------------

### log file
[//]: -------
The `.imp.log` file provides a summary of the imputation done and the sanity check results. 

### out file
[//]: -------
The `.imp.out` file has the following columns:

- `SNP` SNP-ID
- `Chr` Chromosome (only 1 to 22 right now)
- `Pos` Position (hg19, or same as in reference panel)
- `Z.imp` Imputed Z statistics
- `N.imp` Estimation of N (only if missingness was set to `TRUE`)
- `r2.pred` Imputation quality (adjusted as in ...)
- `A1` Reference allele
- `A2` Effect allele
- `window.nbr` in which window imputed
- `lambda` lambda used to impute

if `Z.imp NA` and `r2.pred 0` means that there was not tag SNP.

## References
[//]: -------
** Pasaniuc, B. and Price, A. L. (2016).** *Dissecting the genetics of complex traits using summary association statistics.* bioRxiv, 18(2):117–127.

** Pasaniuc, B., Zaitlen, N., Shi, H., Bhatia, G., Gusev, A., Pickrell, J., Hirschhorn, J., Strachan, D. P., Patterson, N., and Price, A. L. (2014). ** *Fast and accurate imputation of summary statistics enhances evidence of functional enrichment.* Bioinformatics, 30(20).

