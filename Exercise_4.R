rm(list = ls())

library(dplyr)

state_space = function(para) {
  
  list(T = para[1] %>% 
         rbind(1) %>%
         cbind(matrix(0,2,1)),
       Z = 1 %>% cbind(0),
       H = sqrt(para[2] %>% abs()) %>% rbind(0),
       D = sqrt(para[3] %>% abs()))
}

generate_data = function(N = 2500) {
  
  N = 2500
  x = arima.sim(n= N, list(ar = c(0.9)), sd = sqrt(0.5)) %>% 
    as.matrix()
  
  y = x + rnorm(n = N, mean = 0, sd = sqrt(3))
}

Kalman_filter = function(para, p = 1, q = 0, y, alpha) {
  
  N = nrow(y %>% as.matrix())
  m = max(p,q)

  const <- state_space(para)
  
  Z = const$Z
  T = const$T
  H = const$H
  G = matrix(1,1,1)
  
  v = matrix(0, N, 1)
  J = matrix(1, 2, 1)
  P = matrix(1, ncol(Z), ncol(Z))
  F = matrix(0, N, 1)
  
  y <- y %>% as.matrix()
  
  
  for(t in 1:N) {
    
    v[t,] <- y[t] - Z%*%alpha[,t]
    F[t] <- Z%*%P%*%t(Z) + G%*%t(G)
    
    K <- (T%*%P%*%t(Z) + H%*%t(G))%*%(F[t]^-1)
    L <- T - K%*%Z
    
    if(t < N){
      alpha[,(t+1)] <- T%*%alpha[,t]+K*v[t,]}
    
    P <- P*T%*%t(L) + H%*%t(J)
    
    J <- H - K%*%G  
    
  }
  
  list(v = v, F = F, alpha = alpha)
  
}

value_function = function(para, p = 1, q = 0, y = generate_data(N), alpha, epsilon) {
  
  N = y %>% as.matrix() %>% nrow()
  
  val = Kalman_filter(para, p, q, y, alpha = alpha)
  
  F = val$F
  v = val$v
  
  sigma = var(epsilon)
  
  likelihood = matrix(0, N, 1)
  
  options(warn = -1)
  for(t in 1:N){
    
    likelihood[t] = -1/2*(log(pi) + log(sigma*F[t]) + (v[t])^2/(sigma*F[t]))
    
  }
  options(warn = 0)
  
  sum(likelihood, na.rm = TRUE)
}

para = matrix(0.1, 3, 1)
N = 2500 
alpha = matrix(1, nrow = 2, ncol = N)
y = generate_data(N)
Z = state_space(para)$Z

epsilon = (t(-Z%*%alpha)+y)


res <- optim(par = para, 
             fn = value_function,
             method = c('Nelder-Mead'),
             control = list(trace = 1,
                            fnscale = -1),
             y = y,
             p = 1,
             q = 0,
             alpha = alpha,
             epsilon = epsilon)

res$par[2:3] <- res$par[2:3] %>% abs() 

Z = state_space(res$par)$Z
D = state_space(res$par)$D
T = state_space(res$par)$T
H = state_space(res$par)$H
A <- Kalman_filter(para = res$par, p = 1, q = 0, y = y, alpha = alpha)$alpha
C = (t(-B$Z%*%A)+y)

x <- matrix(0, 100, 1)

for(t in 2400:2500) {
  
  x[(t-2399)] <- Z%*%(A[,t] %>% as.matrix()) +D%*%C[t]
  if(t > nrow(A)) {break}
  A[,(t+1)] <- T%*%A[,t] + H%*%rnorm(1,0,1)
  
}


