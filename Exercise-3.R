rm(list = ls())

library(dplyr)



State_space = function(para, p, q) {
  
  m = max(p,q)
  
  phi = para[1:p]
  theta = para[(p+1):(length(para))]
  
  if(length(phi) > length(theta)) {
    
    theta = theta %>% 
      rbind(matrix(1, nrow = length(phi) - length(theta), 1))
    
  } else {
    phi = phi %>% as.matrix() %>% 
      rbind(matrix(1, nrow = length(theta) - length(phi), 1))
  }
  
  
  Z = matrix(0, nrow = 1, ncol = m)
  Z[1] = 1
  
  temp = diag(m-1) %>% 
    rbind(matrix(0, nrow = 1, ncol = m-1))
  T = phi %>% cbind(temp)
  
  H = phi + theta
  
  G = 1
  
  list(Z = Z, T = T, H = H, G = G)
  
}

Kalman_filter = function(para, p, q, y, alpha) {
  
  N = nrow(y %>% as.matrix())
  m = max(p,q)
  
  const <- State_space(para, p, q)
  
  Z = const$Z
  T = const$T
  H = const$H
  G = const$G
  
  v = matrix(0, N, 1)
  J = matrix(1, m, 1)
  P = matrix(1, m, m)
  F = matrix(0, N, 1)
  
  y <- y %>% as.matrix()
  
  
  for(t in 1:N) {
    
    v[t,] <- t(y[t] - Z%*%alpha[,t])
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

value_function = function(para, p, q, y, alpha, epsilon) {
  
  N = y %>% as.matrix() %>% nrow()
  
  val = Kalman_filter(para, p, q, y, alpha)
  
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



N = 1500
y = arima.sim(n= N, list(ar = c(0.9, -0.2), ma = c(0.3)), sd = sqrt(2)) %>% 
  as.matrix()

p = 2
q = 1

para = matrix(0.6/max(p,q),p+q,1)

m = max(p,q)
 
alpha = matrix(1, nrow = m, ncol = N)
const = State_space(para, p, q)
Z = const$Z
G = const$G

epsilon = (G^-1)*(t(-Z%*%alpha)+y)
U.epsilon = (G^-1)*(t(-Z%*%alpha)+y)
temp = matrix(1, N, 1)

i = 1
while(sum(abs(U.epsilon-temp), na.rm = TRUE) > 0.5) {
  temp <- epsilon
  epsilon <- U.epsilon
  
  
res <- optim(par = para, 
      fn = value_function,
      method = c('Nelder-Mead'),
      control = list(trace = 0,
                     fnscale = -1),
      y = y,
      p = p,
      q = q,
      alpha = alpha,
      epsilon = epsilon)

alpha <- Kalman_filter(para = res$par, p = p, q = q, y = y, alpha = alpha)$alpha

U.epsilon = (G^-1)*(t(-Z%*%alpha)+y)


cat('Iteration',i, 'with error', sum(abs(U.epsilon-temp), na.rm = TRUE), 
    'and parameters:')
print(res$par)

i <- i + 1
if(i > 10){break}
}

####### Exercise 2 #########

VIX <- read.csv("~/Desktop/Financial_Econometrics/Exercises/data/VIX.csv", sep = ';')
y <- VIX$Open %>% as.numeric()
N <- y %>% as.matrix %>% nrow()

p = 1
q = 2

para = matrix(0.6/max(p,q),p+q,1)

m = max(p,q)

alpha = matrix(1, nrow = m, ncol = N)
const = State_space(para, p, q)
Z = const$Z
G = const$G

epsilon = (G^-1)*(t(-Z%*%alpha)+y)
U.epsilon = (G^-1)*(t(-Z%*%alpha)+y)
temp = matrix(1, N, 1)

i = 1
while(sum(abs(U.epsilon-temp), na.rm = TRUE) > 0.5) {
  temp <- epsilon
  epsilon <- U.epsilon
  
  
  res <- optim(par = para, 
               fn = value_function,
               method = c('Nelder-Mead'),
               control = list(trace = 0,
                              fnscale = -1),
               y = y,
               p = p,
               q = q,
               alpha = alpha,
               epsilon = epsilon)
  
  alpha <- Kalman_filter(para = res$par, p = p, q = q, y = y, alpha = alpha)$alpha
  
  U.epsilon = (G^-1)*(t(-Z%*%alpha)+y)
  
  
  cat('Iteration',i, 'with error', sum(abs(U.epsilon-temp), na.rm = TRUE), 
      'and parameters:')
  print(res$par)
  
  i <- i + 1
  if(i > 25){break}
}



