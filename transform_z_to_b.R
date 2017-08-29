## Transform the imputed Z-statistics to imputed effect sizes (b)
## --------------------------------------------------------------
## Assumptions: - if your gwas Z-stats originates from a linear regression with a standardised outcome (variance(outcome) = 1)
##              - if the effect sizes are small 

## This function calculates b_imp and se_imp
f.z2b <- function(z, af, n)
{
    # z = imputed z statistics
    # af = allele frequency
    # n = sample size (effective)
  
    se.b <- 1/sqrt(2* af * (1-af) * n)
  
    b <- z * se.b
    
    return(c(b, se.b))
  
  
}

## example
## ---------
output <- read_tsv("output.txt") ## reading ssimp output
N.max <- 1000 ## maximum sample size in the study
output$N_imp <- output$r2.pred * N.max ## calculating the effective sample size


## applying z-to-b function from above
output.b.se <- lapply(1:nrow(output), function(x) 
  f.z2b(output$z_imp[x], output$maf[x], output$N_imp[x])
)

## assigning b_imp and se_imp to dataframe
output$b_imp <- sapply(output.b.se, function(x) x[1])
output$se_imp <- sapply(output.b.se, function(x) x[2])

## store out extended output file
write_tsv(output, "output_b.txt")
