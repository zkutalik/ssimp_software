
	
## Minimal example
`ssimp --data gwas.txt --refpanel EUR,EAS` will generate a `gwas.imp.log` file and a `gwas.imp.out` file.
	
	
	
	
	
	
	
	

## Input 
`--data [filename.txt]`, named columns "SNP", "Z", "N" (Z is numeric, N is numeric, SNP is a character, either rsid or chr:pos.hg19) space separated. If column names differ, define them in `--names`

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
	
**If no impute.range, not impute.snps is chosen, then all variants in the refpanel are imputed that are not provided as tag SNPs.**













## Output
There is an `.log` file providing the output and possible warning messages. The `.out` file has
the following columns:

- `SNP`
- `Chr` Chromosome (only 1 to 22 right now)
- `Pos` Position (HG19)
- `SNP2` Chr:Pos
- `bst.imp` standardized effect sizes
- `P.imp` P-value
- `Z.imp` z
- `N.est` estimation of N
- `r2.pred` impqual
- `A1` reference allele
- `A2` effect allele
- `window.nbr` in which window imputed

`.out.not` lists all the SNPs that were not found in the reference panel 
`.out.tnot` lists all the tSNPs that were not used as tSNPs (but needed): not found in the reference panel or had an NA in Z

if `P.imp NA` and `r2.pred 0` means that there was not tag SNP.

## Pseudocode

#### Initial checks
- check if SNP is rsid or chr:pos
- exclude all lines with NA either in SNP, Z or N, ref.allele, effect.allele
- everthing to upper case

#### What is a tag SNP, what variants to impute, calculate windows
1. Vector of valid tag SNPs, we define it as `tag.snps`
2. Vector of SNPs to impute (consider arguments --impute.snps` and `--impute.range`), we define it as variable `imp.snps`
3. Calculate `--window.core.length` of 1e6 core length

#### LD structure
1. What reference panel to use (`--refpanel`)
2. Build correlation structure for selected SNPs (`tag.snps` and `imp.snps` and `--impute.maf`)`

#### Swap sign's
1. Based on the Alleles of the reference panel, swap the sign of the summary statistics of interest (if `b` or `Z`).

#### For each window `k`
1. do imputation


#### wrap up all imputations

#### checks





## checks
[] names columns input

[] tSNP not in reference panel

[]  SNP not in reference panel

[]  ref alleles not swaped as in reference panel

[] A1 and A2 are not matching with reference panel

[] NA in Z of tSNP

[]  two tSNPs or SNPs same Pos and chr

[]  !is.na(SNP) but is.na(SNP2)

[]  types of missings: NA, "-"



## How to store your own LD structure?
I think its best if we provide the LD structure and MAF summaries for each population in http://www.internationalgenome.org/faq/which-populations-are-part-your-study/, then we can calculate the LD matrix for, say, EUR. 

Alternatively, users can provide a path to their own LD structure (estimated from...), and the summary of the variants (for allele swapping needed)











~~Strike through~~
[] task lists

```javascript
function
```

	code1
	code2
	function
