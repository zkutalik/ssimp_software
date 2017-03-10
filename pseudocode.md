~~Strike through~~
[] task lists

```javascript
function
```

	code1
	code2
	function
	
## Example
	ssimp --data gwas.txt --out res.txt --refpanel EUR,EAS

## Input 
	--data filename, named columns "SNP", "Z", "N" (Z is numeric, N is numeric, SNP is a character, either rsid or chr:pos)
	--out filename
	--names names of columns SNP, Z, N, ref.allele, effect.allele
	--toimp what kind of summary stats to impute: p-value, Z-statistic, beta
	--refpanel 1kg (default) or path to own LD structure (see below)
	--pop.1kg (nodefautl) EUR, ASN, ...
	--lambda 1/sqrt(n) (default) or number or "optimize" (loo)
	--impute.maf 1/n (default): everthing above 1/n will be imputed
	--impute.range NULL (default) chr:pos.start-chr:pos.end (if just chr:NA, then the whole chr is imputed)
	--impute.snps NULL (default) txt file with a set of SNPs to impute (no header): rsid or chr:pos (no quotes)
	--window.core.length 1e6 (default) main window
	--window.flanking.length 250e3 (default) flanking region each side

	pos on hg19
	
## checks
[] names columns input

[] tSNP not in reference panel
\item SNP not in reference panel
\item ref alleles not swaped as in reference panel
\item A1 and A2 are not matching with reference panel
\item NA in Z of tSNP
\item two tSNPs or SNPs same Pos and chr
\item !is.na(SNP) but is.na(SNP2)
\item types of missings: NA, "-"

## do
- check if SNP is rsid or chr:pos
- exclude all lines with NA either in SNP, Z or N, ref.allele, effect.allele
- everthing to upper case



## how to store own LD structure?

## Output

There is an \texttt{.log} file providing the output and possible warning messages. The .out file has
the following columns:

- `SNP`
- Chr Chromosome (only 1 to 22 right now)
- Pos Position (HG19)
- SNP2 Chr:Pos
- bst.imp standardized effect sizes
- P.imp P-value
- Z.imp z
- N.est estimation of N
- r2.pred impqual
- A1 reference allele
- A2 effect allele
- window.nbr in which window imputed


\texttt{.out.not} lists all the SNPs that were not found in the reference panel 
\texttt{.out.tnot} lists all the tSNPs that were not used as tSNPs (but needed): not found in the reference panel or had an NA in Z

if P.imp NA and r2.pred 0 means that there was not tag SNP.

