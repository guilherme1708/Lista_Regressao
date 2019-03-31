# Lista 3 regressão

library(ggplot2)

## Exercício 2

#item a

x <- c(-13, -11, -7, -5, -3, -2, 2, 3, 5, 7, 11, 13)
y <- c(-26.10, -22.74, -16.10, -11.50, -9.98, -7.61, 0.63, 0.41, 2.03, 8.89, 14.04, 19.20)

Y <- cbind(y) # matriz de variável resposta
X <- cbind(1,x) # matriz de dados (fixados)
n <- length(Y)  # tamanho da amostra
V <- solve(t(X)%*%X) # (X'X)-¹
D <- V%*%t(X) # (X'X)-¹*X'

beta.MQ <- D%*%Y # matriz Beta0 e beta1 estimados
sigma2.hat <- as.numeric(t(Y)%*%(diag(n) - X%*%D)%*%Y/(n-2)) # variancia estimada não viciada
sigma.hat <- sqrt(sigma2.hat) # desvio padrão estimada não viciada

Matrix.cov <- sigma2.hat*V # matriz de covariancias 

reg.lin <- lm(y~x) # modelo de regressão linear
summary(reg.lin) # resumo da regressõa

dados <- data.frame(x,y)

ggplot(dados, aes(x=x, y=y)) + geom_point() + geom_smooth(method=lm, se=F) + labs(y = "Variável Resposta", x = "Variável explicativa") +
  ggtitle("Diagrama de Dispersão") 

# item b

F <- function(C,c,beta,sigma2,n){ # Função para calcular a estatística F e o p-valor
  k <- nrow(C)
  F.obs <- as.numeric(t(C%*%beta - c)%*% solve( C%*% solve(t(X)%*%X) %*% t(C)) %*% (C%*%beta - c)/(k*sigma2))
  p.valor = 1 - pf(F.obs,k,n-2)
  return(c(F.obs,p.valor))
}

#H_0: B0 = -5
C0 <- matrix(c(1,0),1,2) # Matriz C
c0 <- -5 # Vetor c
F(C0, c0, beta.MQ, sigma2.hat, length(Y)) # estatística F e o p-valor

#H_0: B1 = 0
C0 <- matrix(c(0,1),1,2) # Matriz C
c0 <- 0 # Vetor c
F(C0, c0, beta.MQ, sigma2.hat, length(Y)) # estatística F e o p-valor

# item c

confint(lm(y~x),level=0.90) # intevalo de confiança para o vetor beta com confiança de 90%
c(sigma2.hat*(n-2)/qchisq(0.95,10),sigma2.hat*(n-2)/qchisq(0.05,10)) # intevalo de confiança para o sigma^2 com confiança de 90%

#item d

#H_0: 1/2B0 + 5/4B1 = 0
C0 <- matrix(c(1/2,5/4),1,2) # Matriz C
c0 <- 0 # Vetor c
F(C0, c0, beta.MQ, sigma2.hat, length(Y)) # estatística F e o p-valor

#H_0: B0 = B1 = 1
C0 <- matrix(c(1,0,0,1),2,2) # Matriz C
c0 <- c(0,0) # Vetor c
F(C0, c0, beta.MQ, sigma2.hat, length(Y)) # estatística F e o p-valor

#H_0: B0/B1 = 3
C0 <- matrix(c(1,-3),1,2) # Matriz C
c0 <- 0 # Vetor c
F(C0, c0, beta.MQ, sigma2.hat, length(Y)) # estatística F e o p-valor

#H_0: B0 + 5/2B1 = 0 B1-B0=6
C0 <- matrix(c(1,-1,5/2,1),2,2) # Matrz C
c0 <- matrix(c(0,6),2,1) # Matriz c
F(C0, c0, beta.MQ, sigma2.hat, length(Y)) # estatística F e o p-valor

# Exercício 3

N <- 1000 # Número de simulações
n <- length(Y)  # tamanho da amostra
V <- solve(t(X)%*%X) # (X'X)-¹
D <- V%*%t(X) # (X'X)-¹*X'
gamma <- 0.975 # Confiança do intervalo

c <- qt((1+gamma)/2,10) # Calculo do quartil da t-student com gl=10

g1 <- qchisq((1-gamma)/2,10) # Calculo do quartil da X² com gl=10
g2 <- qchisq(((1+gamma)/2),10) # Calculo do quartil da X² com gl=10

# Matrizes para armazenar os IC's
b0 <- matrix(NA,N,2)
b1 <- matrix(NA,N,2)
s2 <- matrix(NA,N,2)

for (i in 1:N){
  erro <- rnorm(n,0,sqrt(2)) # erro simulado com distribuição N(0,2) e n=12
  yy <- cbind(-4 + (8/5)*x + erro) # vetor coluna de Y simulados dado que os o vetor x está fixo
  beta.S <- D%*%yy # vetor coluna de simuações de beta0 e beta1
  sigma2.hatS <- as.numeric(t(yy)%*%(diag(n) - X%*%D)%*%yy/(n-2)) # vetor de simuações de sigma²
  Matrix.cov <- sigma2.hatS*V # Matriz de covariâncias
  b0[i,] <- c(beta.S[1] -c*sqrt(Matrix.cov[1,1]), beta.S[1]+c*sqrt(Matrix.cov[1,1])) # IC p/ B0
  b1[i,] <- c(beta.S[2] -c*sqrt(Matrix.cov[2,2]), beta.S[2]+c*sqrt(Matrix.cov[2,2])) # IC p/ B1
  s2[i,] <- c(sigma2.hatS*(n-2)/g2,sigma2.hatS*(n-2)/g1) # IC p/ sigma²
}

# Verifica se o verdadeiro valor está nos intervalos
aa <- ifelse(b0[,1]<= -4 & b0[,2] >= -4 ,1,0)
bb <- ifelse(b1[,1]<= 8/5 & b1[,2] >= 8/5 ,1,0)
ss <- ifelse(s2[,1]<= 2 & s2[,2] >= 2 ,1,0)

# Junta a variaável binária com a matriz de IC's
b0 <- cbind(b0,aa)
b1 <- cbind(b1,bb)
ss <- cbind(s2,ss)

mean(b0[,3]) # Proporção dos IC's que tem o verdadeiro valor de B0=-4
mean(b1[,3]) # Proporção dos IC's que tem o verdadeiro valor de B1=8/5
mean(ss[,3]) # Proporção dos IC's que tem o verdadeiro valor de sigma^2=2

