
## Things we need to think about (or already did)

### SNP or CHR:POS

### Duplicated CHR:POS (but SNP1 and SNP2 are different)
what if two SNP1 and SNP2 were merged into one single SNP, therefore have identical CHR and POS, but the Z-scores don't agree with each other?

### reference panels and LD
For now (30 March 2017) we decided to use VCF files as reference panels (vcf = as commonly used, although bgen would be better). Also, we will recalculate the LD structure always from scratch. Otherwise we would need to define windows. 

### Windows: fixed or flexible
- predefined, fixed windows: no
- overlapping windows: ?



## Pseudocode (as Sina imagine's it)

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

