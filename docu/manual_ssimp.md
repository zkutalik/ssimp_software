[//]: ==================================
# Manual of SSimp
[//]: ==================================

This command line software enables summary statistics imputation (SSimp) for GWAS summary statistics. 

<sup>Details to the method can be found in *our-paper (2017)* (as well as *Pasaniuc et al. (2014)*, or more recently in a review on the use of summary statistics by *Pasaniuc & Price (2016)*). </sup>

	


## Minimal example
[//]: -------------------------------
The minimal requirements are: (1) GWAS summary statistics stored in a text file with at least the following columns SNP-id, Z-statistic, reference allele and risk allele and at least one row, and (2) the path to the reference panel. 

<sup>All other arguments have defaults defined (see below).</sup>

`bin/ssimp --gwas data/my_gwas.txt --ref ref/my_reference_panel.vcf` will generate a file `my_gwas.txt.ssimp.txt`, containing the imputation results, and a log file called `my_gwas.log`.

<sup>If P-values are provided instead Z-statistics, there needs to be an extra column containing the effect sizes (the P-value will be turned into a Z-statistics, and therefore needs a negative sign if the effect size is negative). </sup>
	

## Arguments
[//]: -------------------------------
Here are all parameters listed. Each argument has: a default value defined (in `[brackets]`), valid options listed and a definition of the argument given. Note that arguments can be shortend, e.g. `--wind` instead of `--window.width`

`--gwas [no default]`, filename of the GWAS dataset, extension (e.g. `.txt`) does not matter, nor does text separator (e.g. `\t`).  Columns need to be named after common conventions (see file `../header_translation.md` **@Aaron: is this done already?**). Missings have to be marked as `NA` or left empty. QUICKTEST, SNPTEST, METAL AND PLINK output files will be automatically recognised. The minimal set of columns that should be provided, are: SNP-id, Z-statistics, reference allele and risk allele. For more info on possible sets of columns, see section `GWAS dataset` below.

`--ref [no default]` path to vcf file (same folder should contain the `tbi` file). **Important** >> currently, there is a default to use the uk10k on HPC1.

`--sample.names [no default]`  The argument can be used in two ways: (1) providing a text file (no header) with sample id's separated by new lines `--sample.names filename.samples.txt` or (2) providing a file with at least a sample id column and a column to constrain on: `--sample.names  x/f/e=v` does a lookup in file 'x' (which has a header), but filters on field 'e' being equal to value 'v', and then
use the sample names in column 'f'. An example of the latter is: `/data/sgg/aaron/shared/ref_panels/1kg/integrated_call_samples_v3.20130502.ALL.panel/sample/super_pop=EUR`. **Important** >> currently this argument is defaulted to `super_pop=EUR`. 

`--out [gwasfilename.ssimp.txt]` string. Filename in which to store the imputation results. If not defined, it will be the gwas filename + `.ssimp.txt`. 

`--log [gwasfilename.log]` string. Filename in which to store the log file. If not defined it will be the gwas filename + `.log`. If set to `FALSE`, then no log file is produced. (**TBD** >> this option is not supported yet)

`--impute.range [no default]` Should have the form of `CHR:pos.start-CHR:pos.end`, with `CHR` being the chromosome number, `pos.start` the start position and `pos.end` the end position, e.g. `1:10000-1:30000`. If `CHR`, then the single chromosome `CHR` is imputed. For `CHR-CHR`, a range of chromosomes are imputed, e.g. `1-5` chromosome 1 to chromosome 5 are imputed. (**TBD** >> For chromosome `X`, `Y` and `MT`, text or numbers (23, 24, 25) can be used.)

`--tag.snp [no default]` filename with list of tags (no header). For magic in bash see `Note` below.

`--impute.snp [NULL]` filename to define SNPs to impute (no header). For magic in bash see `Note` below.

`--lambda [2/sqrt(n)]` numeric value or string (`2/sqrt(n)`, `optimise`), n are the number of individuals in the reference panel. `optimise` not yet implemented. Lambda controls the shrinking of the correlation matrix (lambda = 0 applies no shrinking, lambda = 1 turns the correlation matrix into the identity matrix). (**TBD** >> string not implemented yet, so if you want to have 2/sqrt(2) you have to compute it yourself)

`--impute.maf [0]` numeric value. Lower MAF limit for SNPs to be imputed: everything above and equal this threshold will be imputed.

`--tag.maf [0]` numeric value. Lower MAF limit for tag SNPs: everything above and equal this threshold will be used as tag SNPs. 

`--window.width [1e6]` numeric value. Core window length.

`--flanking.width [250e3]` numeric value. Flanking space left and right side of the core window.
		
`--missingness [TRUE]` logical. Enables variable sample size approach. This is automatically set to `FALSE` if `N` is not provided or `N` is set to `NA`. (**TBD** >> missingness is not implemented it)

### Note	
[//]: -------
- If `impute.range` and `impute.snps` are not defined, then all variants in the reference panel are imputed that are not provided as tag SNPs.
- The option `missingness` is automatically set to `FALSE` if `N` is not provided or `N` is set to `NA`.
- Odds ratios need to be provided as Z-statistics or, alternatively, be log-transformed into effect sizes.
- Magic tipp in bash to produce a file within the command line: `--impute.snp <(echo rs5753220 rs5753231 rs5753236 rs5753259 rs5753260 rs5753263 rs5753268 rs5753271 rs5753272 rs5753281 rs5753284 rs5753285 rs5753290 rs5753298 | tr ' ' '\n')`

## GWAS dataset
[//]: -------------------------------
- Column names are automatically recognised using commonly used names (see `../header_translation.md`). Missing values should be marked as `NA` or left empty. 
- The minimal columns required are `SNP`, `A1`, `A2`, `Z`. If `Z` is not present, but `P` and `b` are, `Z` is calculated through `P` and `b`. Alternatively, if `b` and `SE` are present, then it is also possible to calculate `Z` via `b` and `SE`. 
- Positions should match the positions in the reference panel (e.g. both hg19). 
- It is recommended to provide the sample size (N), as incorporating missingness leads to a more accurate estimate. 
- SNP names should be named so they match the SNP-id in the reference panel. 
- If case positions in the GWAS file do not match the reference panel positions, use use LiftOver as a command line tool: http://genome.ucsc.edu/cgi-bin/hgLiftOver

## Reference panel
[//]: -------------------------------
Filename specified as `ref/chr{CHRM}.vcf.gz`, with `CHRM` as the placeholder if the vcf.gz files are split up for each chromosome. The same folder should contain also the `.tbi` file(s). **@Aaron** is this correct?

## Run-time
[//]: -------------------------------
**TBD**

## Technical aspects of summary statistics imputation
[//]: -------------------------------
For more details on the summary statistics imputation method, please see our paper (2017). 

`insert here a brief recap of the method`

- If SNP-ID are present in the GWAS, then pos copied from reference panel. 
- If SNP-ID are not present, then the combination of Chr:Pos:A1:A2 are taken as identifier
- Because either SNP-ID or Chr:Pos:A1:A2 are used as identifier, it is also possible to impute indels.
- Z-statistics are imputed, along with `N_imp` (an estimate for the sample size) and `r2_pred` (adjusted imputation quality).
- To speed up computation, we use a sliding window approach (`--window.width` and `--flanking.width`). SNPs to be imputed are assigned to one window.

## Output
[//]: -------------------------------

### log file
[//]: -------
The `.log` file is a copy of what is printed to the console. 

### out file
[//]: -------
The `.imp.out` file has the following columns:

- `SNP` SNP-ID
- `Chr` Chromosome (only 1 to 22 right now)
- `Pos` Position (same build as reference panel)
- `Z_imp` Imputed Z statistics (see below)
- `N_imp` Estimation of N (only if missingness was set to `TRUE`)
- `r2_pred` Imputation quality (adjusted as in ...)
- Reference allele (same column name as in the GWAS file)
- Effect allele (same column name as in the GWAS file)
- `lambda` lambda used to impute
- `origin` GWAS or SSimp, depending if the SNP was a tag SNP or an imputed SNP. 
- `Z_reimputed` for the first window, all tag SNPs will be re-imputed (sanity check).

`Z_imp` reports the imputed Z-statistics for SNPs that were imputed (`origin = SSimp`), as well as the GWAS Z-statistics for tag SNPs (`origin = GWAS`). 

## References
[//]: -------
**Pasaniuc, B. and Price, A. L. (2016).** *Dissecting the genetics of complex traits using summary association statistics.* Nature Reviews Genetics.

**Pasaniuc, B., Zaitlen, N., Shi, H., Bhatia, G., Gusev, A., Pickrell, J., Hirschhorn, J., Strachan, D. P., Patterson, N., and Price, A. L. (2014).** *Fast and accurate imputation of summary statistics enhances evidence of functional enrichment.* Bioinformatics.

