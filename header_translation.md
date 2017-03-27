#  Header translation

### There are at least four generic input types that we should consider:
- `METAL` : tab-separated
- `plink`: tab-separated
- `snptest`: space-separated
- `quicktest`: space-separated

### possible entries for the crucial variables are:
- `SNP` rsid, MarkerName, snp, SNP, id1, Marker, MARKER
- `Chr` CHR, chromosome
- `Pos` BP, position
- `Allele.ref` Allele1, AlleleA, A1
- `Allele.eff` Allele2, AllleleB, A2 (?plink??)
- `b`  frequentist_add_beta_1, normal.score.beta, BETA, b, B
- `se` frequentist_add_se_1, normal.score.se, se, SE
- `Z` STAT, Zscore, Z
- `p` frequentist_add_pvalue, normal.score.p, p, P, P-value, p.value, P.value, p-value
- `N`  N, all_total, NMISS
- `EAF` all_maf

### Things to consider
- Plink only operates with A1
- Plink only has NMISS
- Plink has no EAF by default
- quicktest no N, no EAF

Maybe we should operate without EAF, as it is hardly ever reported...
