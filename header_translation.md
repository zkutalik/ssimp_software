#  Header translation

### There are at least four generic input types that we should consider:
- `METAL` 
- `plink`
- `snptest`
- `quicktest`

### possible **lowercase** entries for the GWAS column names are:
- `snp` >> `rsid`, `markername`, `id1`, `marker`
- `chr` >> `CHR`, `chromosome`
- `pos` >> `BP`, `position`
- `a1` >> `allele.ref`, `allele1`, `allelea`
- `a2` >> `allele.eff`, `allele2`, `alleleb`
- `b` >> `frequentist_add_beta_1`, `normal.score.beta`, `beta`
- `se` >> `frequentist_add_se_1`, `normal.score.se`
- `Z` >> `STAT`, `zscore`
- `p` >> `frequentist_add_pvalue`, `normal.score.p`, `p-value`, `p.value`
- `N` >> `N, all_total, NMISS, (meanAA + meanAB + meanBB)
[//]: `EAF` >> `all_maf`

### Things to consider
- Plink only operates with A1
- Plink only has NMISS
[//]: - Plink has no EAF by default
[//]: - quicktest no EAF

[//]: Maybe we should operate without EAF, as it is hardly ever reported...
