#####################################################################################
## Author: Sina Rueeger [sina *.* rueger *a*t* unil *.* ch]
## Project: 
##        
## Time-stamp: <[checkout.R] by SR Wed 22/03/2017 13:10 (CET)>
##
## Description:
## 
##
## History:
## 
#####################################################################################
dir()

ss.random <- read.table("gwas.example.txt", nrow = 10, header = TRUE)

ss.qt <- read.table("quicktest1d.out", nrow = 10, header = TRUE)
ss.st1 <- read.table("snptest_genotype_gen.out", nrow = 10, header = TRUE)
ss.st2 <- read.table("snptest_genotype_bed.out", nrow = 10, header = TRUE)
ss.st3 <- read.table("snptest_genotype_imp.out", nrow = 10, header = TRUE)
ss.metal <- read.table("metal_example.TBL", nrow = 10, header = TRUE)

RANDOM

            MARKER Allele1 Allele2 ProbeName        Z   N            P
   5:96252589:2142       T       C      2142 16.42881 214 1.189705e-60

QUICKTEST
 id1  id2 pos alleleA alleleB      meanAA meanAB meanBB  rSqHat
1   1 snp1   1       A    T 1.42109e-14 12.761 62.239 0.86665
  normal.mean.beta normal.mean.se normal.mean.p normal.score.beta
1        -0.510611       0.338152      0.135359         -0.513378
  normal.score.se normal.score.p normal.score.info
1        0.344322       0.135965          0.788329

SNPTEST
 alternate_ids          rsid chromosome position alleleA alleleB index
1             6 Affx-28523749         NA 33734763       A       C     1
  average_maximum_posterior_call info cohort_1_AA cohort_1_AB cohort_1_BB
1                              1    1        4974       39654       75111
  cohort_1_NULL all_AA all_AB all_BB all_NULL all_total  all_maf
1           347   4974  39654  75111      347    120086 0.207125
  missing_data_proportion frequentist_add_pvalue frequentist_add_info
1               0.0014448              0.0010839                   NA
  frequentist_add_beta_1 frequentist_add_se_1 comment
1              0.0165465           0.00506338      NA

METAL
        MarkerName Allele1 Allele2 Weight  Zscore    P.value         Direction
1    8:40631655:539       a       g   7757 -13.832  1.627e-43 -??--???-?----?+-



    SNP or SNP2 or (CHR + POS)
    Zscore or (P + Beta) or (beta + se)
    N
    Allele1, allele2

