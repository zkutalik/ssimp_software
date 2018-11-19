[//]: ==================================
# Manual of SSimp
[//]: ==================================


## Minimal example
[//]: -------------------------------
The **minimal requirements** for `ssimp` to run are:
1. GWAS summary statistics stored in a text file with at least the following columns SNP-id (e.g. `MarkerName`), Z-statistic (e.g. `Z`), reference allele (e.g. `a1`) and risk allele (e.g. `a2`) and at least one row. 
2. The path to the reference panel.
3. The path to the output file.

`ssimp my_gwas.txt ~/reference_panels/my_reference_panel.vcf output.txt` 

will generate a file `output.txt`, containing the imputation results. This is identical to 

`ssimp --gwas my_gwas.txt --ref ~/reference_panels/my_reference_panel.vcf --out output.txt`
	
An **executable example** is:

`ssimp --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt`

Where `ssimp` is the newest version in the `compiled/` or `bin/`folder. 

## Options
[//]: -------------------------------
The options `--gwas`, `--ref` and `--out` are required arguments. Some options (e.g. `--lambda`) have `[defaults]` defined, and other options (e.g. `--log`) are entirely optional.

Note that arguments can be shortenend, e.g. `--wind` instead of `--window.width`.

`--download.build.db` downloads the build database, and exits.

`--download.1KG` downloads the 1000 genomes reference panel (18GB) and the build database, and then exits.

`--gwas [no default]`, path to the GWAS dataset. The file's extension (e.g. `.txt`) does not matter. The delimiter (comma, space or tab) is detected automatically. Common column names are recognised automatically (for details see section `GWAS dataset`). The minimal set of columns that should be provided, are: SNP-id, Z-statistics, reference allele and risk allele. Missings have to be marked as `NA` or left empty.

`--ref [no default]` path to `.vcf` and [`.tbi` file(s)](https://genome.ucsc.edu/goldenpath/help/vcf.html). Filename specified as `chr{CHRM}.vcf.gz`, with `CHRM` as the placeholder if the vcf.gz files are split up for each chromosome, e.g. `--ref ~/reference_panel/1000genomes/ALL.chr{CHRM}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz` (for this example, you need to have the `tbi` files `~/reference_panel/1000genomes/ALL.chr{CHRM}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi` available too).

`--out [no default]` string. Filename to store the imputation results. 

`--sample.names [no default]`  The argument can be used in two ways: (1) providing a text file (no header) with sample id's separated by new lines `--sample.names filename.samples.txt` or (2) providing a file with at least a sample id column and a column to constrain on: `--sample.names  x/f/e=v` does a lookup in file 'x' (which has a header), but filters on field 'e' being equal to value 'v', and then
use the sample names in column 'f'. An example of the latter is: `integrated_call_samples_v3.20130502.ALL.panel/sample/super_pop=EUR`. 

`--log [no default]` string. Filename in which to store the log file - this is simply a copy of whatever is printed to the console.

`--impute.range [no default]` Should have the form of `CHR:pos.start-CHR:pos.end`, with `CHR` being the chromosome number, `pos.start` the start position and `pos.end` the end position, e.g. `1:10000-1:30000`. If `CHR`, then the single chromosome `CHR` is imputed. For `CHR-CHR`, a range of chromosomes are imputed, e.g. `1-5` chromosome 1 to chromosome 5 are imputed.

`--impute.maf [0]` numeric value. Lower MAF limit for SNPs to be imputed: everything above and equal this threshold will be imputed.

`--impute.snp [NULL]` filename to define SNPs to impute (each SNP has a new line, no header). For magic in bash see `Note` below.

`--tag.range [no default]` same as `impute.range`, but applied to tag SNPs.

`--tag.maf [0]` same as `impute.maf`, but applied to tag SNPs.

`--tag.snp [no default]`same as `impute.snp`, but applied to tag SNPs.

`--lambda [2/sqrt(n)]` numeric value or string (`2/sqrt(n)`), n are the number of individuals in the reference panel. Lambda (????) controls the shrinking of the correlation matrix (lambda = 0 applies no shrinking, lambda = 1 turns the correlation matrix into the identity matrix).

`--window.width [1000000]` numeric value. Core window length. See illustration below.

`--flanking.width [250000]` numeric value. Flanking space left/right of the core window. See illustration below.
		
`--missingness [none]` string: `ind` (recommended), `dep`. Enables variable sample size approach. `ind` stands for independent, and `dep` for dependent. 

`--reimpute.tags` reimpute every tag in every window. Otherwise, only 100 tags are reimputed.


### Multiprocessing
[//]: -------
Multiprocessing mode is not implemented, hence to speed up computation we recommend splitting up the job to smaller chromosomal chunks using the option `--impute.range` for smaller chunks or even chromosomes:

`ssimp --gwas gwas/small.random.txt --ref ref/small.vcf.sample.vcf.gz --out output.txt --impute.range 22:18000000-22:18100075`

There is an untested bash script called [extra/ssimp_chunks.sh`](../extra/ssimp.sh) that computes chunks with an R script. 


### Tips	
[//]: -------
- If `impute.range` and `impute.snps` are not defined, then all variants in the reference panel are imputed (including the tag SNPs of the first window for sanity checks, see section `output` below).
- Magic tip in bash to produce a file within the command line: `--impute.snp <(echo rs5753220 rs5753231 | tr ' ' '\n')`. Have a look at the [examples](examples.md).


## GWAS dataset
[//]: -------------------------------
- Column names are automatically recognised using commonly used names. See subsection below.
- Odds ratios need to be provided as Z-statistics or, alternatively, be log-transformed into effect sizes (`b`).
- The minimal columns required are `SNP`, `A1`, `A2`, `Z`. If `Z` is not present, but `P` and `b` are, `Z` is calculated through `P` and `b`. 
- Missing values should be marked as `NA` or left empty. Sample size cannot contain `NA`. 
- Positions should match the positions in the reference panel (e.g. both hg19). 
- It is recommended to provide the sample size (N) and set `--missingness dep`, as incorporating missingness leads to a more accurate estimate. 
- SNP names should be named so they match the SNP-id in the reference panel. E.g. if the reference panel uses `chr:pos`, the GWAS should have the same SNP identifier.
- In case no SNP identifiers are present, use `chromosome` and `position`, but provide a `SNP` column with empty entries.
- If case positions in the GWAS file do not match the reference panel positions, use use LiftOver as a command line tool: http://genome.ucsc.edu/cgi-bin/hgLiftOver.

### Automatic header recognition

Column names listed above are not case sensitive. E.g. P-value column can be named `p` or `P`.

The column identifying the
- **SNP-id** should be named: `ID`, `rnpid`, `snpid`, `rsid`, `MarkerName`, `snp`, `id1` or `marker`.
- **reference allele** should be named: `REF`, `a1`, `Allele1`, `AlleleA` or `other_allele`.
- **effect allele** should be named: `ALT`, `a2`, `Allele2`, `AlleleB` or `effect_allele`.
- **Z-statistic** should be named: `z.from.peff`, `z`, `stat`, `zscore` or `z.score`. 
- **P-value** should be named: `p`, `P-value`, `P.value`, `PVALUE`, `frequentist_add_pvalue` or `normal.score.p`.
- **Effect size** should be named: `b`, `beta`, `ALT_EFFSIZE`, `frequentist_add_beta_1` or `normal.score.beta`.
- **Standard error** should be named: `se`, `frequentist_add_se_1` or `normal.score.se`.
- **Sample size** should be named: `n`, `NMISS` or `all_total`
- **Chromosome** should be named: `chr`, `chromosome`, `chrm` or `#CHROM`.
- **Position** should be named: `POS`,`position` or `BP`.


## --window.width and --flanking.width
[//]: -------------------------------
To speed up computation, we use a sliding window approach (`--window.width` and `--flanking.width`). SNPs to be imputed are assigned to one window.
![Caption for the picture.](visuals/visualisation_width.jpeg)

On each chromosome you can only impute variants that are in the range between **min(position)+flanking.width** to **max(position)-flanking.width**, the position being the position of all tag SNPs on a specific chromosome.

## Reference panel
[//]: -------------------------------

More info on handling reference panel data can be found in [usage message](doc/usage.txt).


## Technical aspects
[//]: -------------------------------

- If SNP-ID are present in the GWAS, then positions (bp) are copied from reference panel. 
- If SNP-ID are not present (e.g. filled with `.`), then the combination of Chr:Pos:A1:A2 are taken as identifier
- Because either SNP-ID or Chr:Pos:A1:A2 are used as identifier, it is also possible to impute indels.

## Output
[//]: -------------------------------

### log file
[//]: -------
The `.log` file is a copy of what is printed to the console. (not supported yet)

### out file
[//]: -------
The output file specified in `--out` file has the following columns:
- `chr` Chromosome (1-22, X)
- `pos` Position (same build as reference panel)
- `z_imp` Imputed Z statistics (see below)
- `source` GWAS or SSimp, depending if the SNP was a tag SNP or an imputed SNP. 
- `SNP` SNP-ID
- Reference allele (same column name as in the GWAS file)
- Effect allele (same column name as in the GWAS file)
- `maf` minor allele frequency in reference panel
- `r2.pred` Imputation quality (as defined in the Method outline)
- `lambda` lambda used to penalise.
- `Z_reimputed` imputed Z-statistics for tag SNPs for the first window (sanity check).
- `r2_reimputed` imputation quality for the imputed tag SNPs of for the first window (sanity check).
- `N.imp` Effective sample size after imputation (maximal sample size times the imputation quality (`r2.pred`))
- `P.imp` Imputed P-value [2 x CDF(-|`z_imp`|)] (CDF = cumulative distribution function of the normal distribution)
- `bst.imp` Imputed standardised effect size [`z_imp`/sqrt(`N.imp`)]

Note that `z_imp` reports the imputed Z-statistics for SNPs that were imputed (`origin = SSimp`), as well as the GWAS Z-statistics for tag SNPs (`origin = GWAS`). 

The imputation quality (`r2.pred`) can become zero, which results in `N.imp=0` and `bst.imp=inf`.

Check the R-file [`extra/sanitycheck_reimputed.R`](../extra/sanitycheck_reimputed.R) to find out more about the columns `Z_reimputed`, `r2_reimputed` and `source`.	

## Method outline
[//]: -------------------------------
The principle of *summary statistics imputation* is to combine the available summary statistics from a set of markers and the fine-scale LD structure in the same region in order to infer association summary statistics of unobserved variants. We can formally write this using the conditional expectation of a multivariate normal distribution. 

![Summary statistics equation](visuals/eq_main.jpg)

Here we aim to impute the **Z-statistic** of an untyped SNP *u*, given the Z-statistics of a set of tag SNPs called *M* (LHS of the equation). The RHS of the equation contains **c** (representing the correlations between SNP *u* and all the tag SNPs *M*), **C** (the pairwise correlations among the tag SNPs), and the Z-statistics of a set of tag SNPs *M*. Both, **c** and **C** are regularised using the option `--lambda`. *u* can be extended to a vector. SNPs to impute are in a core window (of size `--window.width`) and *M* is the set of tag SNPs within the core window +/- a flanking region (of size `--window.width`).

### Imputation quality
The estimated imputation quality is a measure that varies between 0 and 1, with 0 reflecting poor and 1 perfect imputation. We use an adjusted R<sup>2</sup> estimation that is additionally corrected by the effective number of tag SNVs `p_eff` (`n`=number of individuals in the reference panel).

![Imputation quality](visuals/eq_impqual.jpg)

### Variable missingness
To account for variable sample size in summary statistics of tag SNVs, we use an approach to down-weight entries in the **C** and **c** matrices for which summary statistics was estimated from a GWAS sample size lower than the maximum sample size in that dataset. Using such down-weighting, turns matrices **C** and **c** into **D** and **d**. Subsequently, the matrices are used for imputation of summary statistics as well as the calculation of the imputation quality. 

### Transform Z-statistics to `b` and `se(b)`
Imputed Z-statistics can be transformed into effect sizes (`b`) and standard errors of effect sizes (`se(b)`).

The calculation depends on the type of model used for your GWAS summary statistics. 

If your Z-statistics originated from a linear regression model, the estimated `se(b)` can be approximated by:

![se(b)](visuals/eq_seb.jpg)

with `q_u` being the allele frequency and `N_u` the sample size of SNP `u`. We can approximate the sample size of an untyped SNP by multiplying the maximum sample size with its imputation quality `N_u = N_max * r2_pred`. 

We provide an R-function called [`extra/transform_z_to_b.R`](../extra/transform_z_to_b.R) that automatically transforms imputed Z-statistics into imputed effect sizes.

### More background on method

For more details on *summary statistics imputation*, see [Rüeger et al. (2018A)](https://wjournals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007371). [Rüeger et al. (2018B)](https://doi.org/10.1101/203927) describes a method to manipulate the LD matrix in order to account for ancestry admixture. 

Most of our extended method builds on [Pasaniuc et al. (2014)](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btu416). 

We also recommend reading the review on the use of summary statistics by [Pasaniuc & Price (2016)](https://www.nature.com/articles/nrg.2016.142).

## How to cite SSimp
[//]: ----------

**Rüeger, S., McDaid, A., Kutalik, Z. (20178).** *Evaluation and application of summary statistic imputation to discover new height-associated loci* PLOS Genetics. https://doi.org/10.1371/journal.pgen.1007371

## References
[//]: -------

**Pasaniuc, B., Zaitlen, N., Shi, H., Bhatia, G., Gusev, A., Pickrell, J., Hirschhorn, J., Strachan, D. P., Patterson, N., and Price, A. L. (2014).** *Fast and accurate imputation of summary statistics enhances evidence of functional enrichment* Bioinformatics. [https://doi.org/10.1093/bioinformatics/btu416](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btu416)

**Pasaniuc, B. and Price, A. L. (2016).** *Dissecting the genetics of complex traits using summary association statistics* Nature Reviews Genetics. [doi:10.1038/nrg.2016.142](https://www.nature.com/articles/nrg.2016.142)

**Rüeger, S., McDaid, A., Kutalik, Z. (2018).** *Evaluation and application of summary statistic imputation to discover new height-associated loci* PLOS Genetics. https://doi.org/10.1371/journal.pgen.1007371.

**Rüeger, S., McDaid, A., Kutalik, Z. (2018).** *Improved imputation of summary statistics for realistic settings* bioRxiv. https://doi.org/10.1101/203927.
