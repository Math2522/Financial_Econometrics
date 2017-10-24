rm(list = ls())
library(dplyr)
library(magrittr)

## Exercise 1 ## 

S <- 10000 # simulations 


AR_generate <- function(Burn, rho, T, sigma, y_start = 1) {
  
  
  y <- matrix(y_start, T + Burn, 1)
  
  for(i in 2:nrow(y)) {
    
    y[i,1] <- y[(i-1),1]*rho + rnorm(n = 1, mean = 0, sd = sqrt(sigma))
    
  }
  
  return(y[(Burn+1):nrow(y),])
  
  
}

rho <- matrix(0,S,1)


for(i in 1:S) {
  
  y <- AR_generate(Burn = 1000, rho = 0.4, T = 1000, sigma = 0.5)
  
   rho[i,1] <- (lm(y~lag(y))$coefficients %>% as.matrix)[2,1]
  
}

# Plotting the kernel density of rho

density(rho) %>% plot

##### Exercise 2 ########### 

Theta <- sqrt(10000)*(rho-0.4)
density(Theta) %>% plot

#### Exercise 3 ######

var(Theta)

# The empirical variance is much higher


##### Exercise 4 #######

rho_larger_sigma <- matrix(0,S,1)


for(i in 1:S) {
  
  y <- AR_generate(Burn = 1000, rho = 0.4, T = 1000, sigma = 2)
  
  rho_larger_sigma[i,1] <- (lm(y~lag(y))$coefficients %>% as.matrix)[2,1]
  
}

density(rho_larger_sigma) %>% plot()


# The distribution is still centered around the true value. There is no apparent change in the variance. 

######### Exercise 5 ########

rho_unit_root <- matrix(0,S,1)


for(i in 1:S) {
  
  y <- AR_generate(Burn = 1000, rho = 1, T = 1000, sigma = 0.5)
  
  rho_unit_root[i,1] <- (lm(y~lag(y))$coefficients %>% as.matrix)[2,1]
  
}

density(rho_unit_root) %>% plot()


