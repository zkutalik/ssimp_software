

head(output)


## adds b_imp and se_imp
## linear regression!
f.z2b <- function(z, af, n)
{
  
  se <- 1/sqrt(2* af * (1-af)* n)
  
    b <- z * se
    
    return(c(b, se))
  
  
}


N.max <- 1000
output$N_imp <- output$r2.pred * N.max

lapply(1:nrow(output), function(x) 
  f.z2b(output$z_imp[x], output$maf[x], output$N_imp[x])
)