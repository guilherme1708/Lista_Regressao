# Lista 5

library(ggplot2)

# Exercio 4

y <- c(17.93,18.82,15.39,14.21,14.16,14.87,9.67,13.41,11.06,8.55,5.12,7.91,
      15.69,17.04,16.21,16.96,14.83,15.07,15.33,14.46,12.70,13.11,11.23,11.32)

x1 <- c(0,0,0,2,2,2,6,6,6,10,10,10,rep(0,12))
x2 <- c(rep(0,12),0,0,0,2,2,2,6,6,6,10,10,10)
x <- matrix(c(x1,x2),24,2)

Y <- cbind(y)
X <- cbind(1,x)

n <- length(Y)  # tamanho da amostra
V <- solve(t(X)%*%X) # (X'X)-¹
D <- V%*%t(X) # (X'X)-¹*X'

beta.MQ <- D%*%Y # matriz Beta0 e beta1 estimados
sigma2.hat <- as.numeric(t(Y)%*%(diag(n) - X%*%D)%*%Y/(n-4)) # variancia estimada não viciada
sigma.hat <- sqrt(sigma2.hat) # desvio padrão estimada não viciada
Matrix.cov <- sigma2.hat*V # matriz de covariancias 

model <- lm(formula = Y ~ x) # modelo de regressão

confint(model) # Intervalos de confiança marginais

# item c

F <- function(C,c,beta,sigma2,n){ # Função para calcular a estatística F e o p-valor
  k <- nrow(C)
  F.obs <- as.numeric(t(C%*%beta - c)%*% solve( C%*% solve(t(X)%*%X) %*% t(C)) %*% (C%*%beta - c)/(k*sigma2))
  p.valor = 1 - pf(F.obs,k,n-2)
  return(c(F.obs,p.valor))
}

#H_0: B1 - B2 = 0
C0 <- matrix(c(0,-1,1),1,3) # Matriz C
c0 <- 0 # Vetor c
f <- F(C0, c0, beta.MQ, sigma2.hat, length(Y)) # estatística F e o p-valor

# item d

x0 <- c(1,10,0)
mu0 <- x0%*%beta.MQ
