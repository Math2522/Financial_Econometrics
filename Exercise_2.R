rm(list = ls())

library(dplyr)
library(magrittr)
library(readxl)
library(statmod)


TurnBAC <- read_excel("~/Desktop/Financial_Econometrics/Exercises/data/TurnBAC.xlsx")


loglik = function(para, data) {
  
  
  omega <- para[1]
  alpha <- para[2]
  beta <- para[3]
  nu <- para[4]
  
  T <- data %>% as.matrix %>% nrow()
  mu <- omega/(1-alpha-beta)*identity(T)
  
  llk <- matrix(0, T, 1)
  
 for(t in 2:T) {
   
   mu[t] = omega + alpha*data[(t-1)] + beta*mu[(t-1)]
   epsilon = data[t]/mu[t]
   llk[t] = -log(data[t]) + nu*log(epsilon*nu) - epsilon*nu - log(gamma(nu))
   
 }  
  
  return(colSums(llk))
  
}

data_generate = function(para, data) {
  
  omega <- para[1]
  alpha <- para[2]
  beta <- para[3]
  nu <- para[4]
  T <- data %>% as.matrix %>% nrow()
  
  mu <- omega/(1-alpha-beta)*identity(T)
  
  for(t in 2:T) {
    
    mu[t] = omega + alpha*data[(t-1)] + beta*mu[(t-1)]
    epsilon = data[t]/mu[t]
  }
  
  return(mu)
}


para = matrix(c(1, 0.2, 0.2, 1))


TurnBAC$Date <- TurnBAC %>%
  select(Date) %>%
  transform(Date = as.Date(as.character(Date), "%Y%m%d"))

ui = matrix(c(1,0,0,0,
            0,1,0,0,
            0,0,1,0,
            0,0,0,1,
            0,-1,-1,0),
            nrow =4,
            ncol = 5) %>% t()
ci = matrix(c(0,0,0,0,-1))

result <- constrOptim(para, 
                      loglik,
                      ui = ui,
                      ci = ci,
                      method = 'Nelder-Mead' ,
                      data = TurnBAC$V_t, 
                      control = list(trace = 1,
                                     fnscale = -1))

mu_hat <- data_generate(result$par, data = TurnBAC$V_t)
epsilon_hat <- TurnBAC$V_t/mu_hat

for(i in c(1,5,22)) {
  
  Box.test(epsilon_hat %>% pgamma(shape = 1/result$par[4]) %>% qinvgauss(), lag = i) %>% print()
  
}


V_hat <- mu_hat*epsilon_hat
plot(V_hat, type = 'l', col = 1)
lines(TurnBAC$V_t, col = 2)

