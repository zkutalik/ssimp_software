## Transform the imputed Z-statistics to imputed effect sizes
## --------------------------------------------------------------
## Assumptions: - if your gwas Z-stats originates from a linear regression with a standardised outcome (variance = 1)
##              - if the effect sizes are small 

## this function calculates b_imp and se_imp
f.z2b <- function(z, af, n)
{
  
    se.b <- 1/sqrt(2* af * (1-af)* n)
  
    b <- z * se.b
    
    return(c(b, se))
  
  
}

## example
## ---------
output <- read.table()
N.max <- 1000
output$N_imp <- output$r2.pred * N.max

lapply(1:nrow(output), function(x) 
  f.z2b(output$z_imp[x], output$maf[x], output$N_imp[x])
)
