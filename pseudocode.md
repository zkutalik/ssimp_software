
	
## Minimal example
`ssimp --data gwas.txt --refpanel EUR,EAS` will generate a `gwas.imp.log` file and a `gwas.imp.out` file.
	
	
	
	
	
	
	
	

## Input 
`--data [filename.txt]`, named columns "SNP", "Z", "N", "Chr", "Pos", "EAF" (Z is numeric, N is numeric, SNP is a character, either rsid or chr:pos.hg19, chr is numeric, pos.hg19 is numeric) space separated. If column names differ, define them in `--names`. Missings have to be marked as `NA`. quicktest, snptest,metal and plink output files will be automatically recognised.

`--pop.1kg [pop]` needs to be defined: abbreviations of 1000genomes populations EUR, ASN, see http://www.internationalgenome.org/faq/which-populations-are-part-your-study/

`--out [filename.imp]` no extension needed, if not defined it will be `filename.imp`

`--names [SNP, Z, N, ref.allele, effect.allele]` vector of names of columns SNP, Z, N, ref.allele, effect.allele

`--refpanel [1kg]` (default) or path to own LD structure (see below)

`--toimp [Z]` what kind of summary stats to impute: p-value (`P`), Z-statistic (`Z`), standardised beta (`b`)

`--lambda [1/sqrt(n)]` if LD is prestored and n is known, or a number or "optimize" (loo)

`--impute.maf [1/n]` lower limit for variants to be imputed: everthing above 1/n will be imputed

`--impute.range [chr:pos.start-chr:pos.end]` NULL (default), if just chr:NA, then the whole chr is imputed

`--impute.snps [snps2impute.txt]` NULL (default) or txt file with a set of SNPs to impute (no header): rsid or chr:pos.hg19 (no quotes)

`--window.core.length [1e6]` main window

`--window.flanking.length [250e3]` flanking region each side
		
`--missingness [TRUE]` enable variable sample size approach. This is automatically set to `FALSE` if `N` is not provided or `N` is set to `NA`.

### Note	
- If no `impute.range`, not `impute.snps` is chosen, then all variants in the refpanel are imputed that are not provided as tag SNPs.
- if your positions are not on hg19, you can use LiftOver as a command line tool: http://genome.ucsc.edu/cgi-bin/hgLiftOver
- The option `missingness` is automatically set to FALSE if `N` is not provided or `N` is set to `NA`.








## Output
The `.log` file provides the output and possible warning messages. The `.out` file has
the following columns:

- `SNP`
- `Chr` Chromosome (only 1 to 22 right now)
- `Pos` Position (HG19)
- `SNP2` Chr:Pos
- `bst.imp` standardized effect sizes
- `P.imp` P-value
- `Z.imp` z
- `N.est` estimation of N (if `missingness` option was set to `TRUE`)
- `r2.pred` impqual
- `A1` reference allele
- `A2` effect allele
- `window.nbr` in which window imputed

`.out.not` lists all the SNPs that were not found in the reference panel 
`.out.tnot` lists all the tSNPs that were not used as tSNPs (but needed): not found in the reference panel or had an NA in Z

if `P.imp NA` and `r2.pred 0` means that there was not tag SNP.

## Pseudocode

#### 1. Initial checks
- check if SNP is rsid or chr:pos
- Input checks: 
	* if `data` is separated by space and a header, missings as NA (no dash, etc.)
	* chr numeric, we cannot impute chr 23
	* if `!is.null(names)` check if column names exist, if `is.null(names)` check if SNP, Pos, Chr, Z, N exists.	
	* `pop.1kg`: check if str is among the ones in http://www.internationalgenome.org/faq/which-populations-are-part-your-study/
	* if `!is.null(refpanel)`: if valid path, and valid LD structure
	* `pop.1kg`: if P, Z or b
	* if `!is.null(lambda)`: numeric, between 0 and 1.
	* if `!is.null(impute.maf)`: if between 0 and 1 
	* if `is.null(impute.range) & is.null(impute.snps)`: impute all
	* if `!is.null(impute.range) & is.null(impute.snps)`: check if range valid
	* if `is.null(impute.range) & !is.null(impute.snps)`: check if valid txt file, no duplicates
	* if `!is.null(impute.range) & !is.null(impute.snps)`: set impute.range to NULL
	* `window.core.length`: numeric and large enough
	* `window.flanking.length`: numeric
	* write small protocol to `.log` file
- rename columns to standard or
- make data input conform (e.g. it its data.snptest, turn it into a common data set)
- check for duplicates in `data`
- if Z, P or b missing, compute
- exclude all lines with NA either in SNP, Z or N, ref.allele, effect.allele
- everything to upper case

#### 2. What is a tag SNP, which variants to impute, calculate windows
1. Vector of valid tag SNPs, we define it as `tag.snps`
2. Vector of SNPs to impute (consider arguments --impute.snps` and `--impute.range`), we define it as variable `imp.snps`
3. Calculate the set of `--window.core.length` of 1e6 core length

#### 3. Calculate LD structure
1. What reference panel to use (`--refpanel`)
2. Build correlation structure for selected SNPs (`tag.snps` and `imp.snps` and `--impute.maf`)`
3. Document in `.log` which tag and imp SNPs could not be used

#### 4. Swap sign's
1. Based on the Alleles of the reference panel, swap the sign of the summary statistics of interest (if `b` or `Z`).

#### 5. For each window `k`
1. do imputation


#### 6. wrap up all imputations

#### 7. End checks
- `0 >= r2.pred.adj <= 1` make sure imputation quality respects the range between 0 and 1.
- impute some known tag SNPs and check whether that worked

#### 8. document all steps in log file


## How to store your own LD structure?
I think its best if we provide the LD structure and MAF summaries for each population in http://www.internationalgenome.org/faq/which-populations-are-part-your-study/, then we can calculate the LD matrix for, say, EUR. 

Alternatively, users can provide a path to their own LD structure of a reference panel or their choice, and the summary of the variants (for allele swapping needed)











~~Strike through~~
[] task lists

```javascript
function
```

	code1
	code2
	function
