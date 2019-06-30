# Lista 7

# Exercício 1
setwd("/home/gui/Área de Trabalho/lista 7")
Dados <- read.csv("data.csv",header = T)

attach(Dados)

n <- length(y)

# item a
plot(y~x, main="Diagrama de dispersão", pch=19)
abline(lm(y~x))

# item d

EF = function(x,theta0,epsilon=10^(-6), it=1000){ # Dados, estimativa inicial, erro, iterações
  
  U <- function(b0,b1,s2){ # Função Escore
    U1 <- 1/s2*sum(x^b1*(y-b0*x^b1))
    U2 <- b0/s2*sum(x^b1*log(x)*(y-b0*x^b1))
    U3 <- -n/2*(1/s2)+sum(2/(2*s2)^2*((y-(b0*x^b1))^2))
    return (c(U1,U2,U3))
  }
  
  IF <- function(b0,b1,s2){ # Matriz informação de Fisher
    IF11 <- -(1/(2*s2)*sum((2*(x^b1*x^b1))))
    IF12 <- 1/(2*s2)*sum((2*(x^b1*log(x)*(y-(b0*x^b1))-x^b1*(b0*(x^b1*log(x))))))
    IF13 <- -(2/(2*s2)^2*sum((2*(x^b1*(y-(b0*x^b1))))))
    IF22 <- 1/(2*s2)*sum((2*(b0*(x^b1*log(x)*log(x))*(y-(b0*x^b1))-b0*(x^b1*log(x))*(b0*(x^b1*log(x))))))
    IF23 <- -(2/(2*s2)^2*sum((2*(b0*(x^b1*log(x))*(y-(b0*x^b1))))))
    IF33 <- -(2*(2*(2*(2*s2)))/((2*s2)^2)^2*sum(((y-(b0*x^b1))^2)) -n/2*(1/s2^2))
    M <- matrix(c(IF11,IF12,IF13,IF12,IF22,IF23,IF13,IF23,IF33),3,3)
    return (M)
  }
  
  erro <- 10
  j <- 0
  b0 <- numeric()
  b1 <- numeric()
  s2 <- numeric()
  b0[1] <-theta0[1]
  b1[1] <-theta0[2]
  s2[1] <-theta0[3]
  
  while(erro > epsilon & j < it){
    j <- j+1
    Aux <- c(b0[j],b1[j],s2[j]) - solve(IF(b0[j],b1[j],s2[j]))%*%U(b0[j],b1[j],s2[j])
    b0[j+1] <- Aux[1]
    b1[j+1] <- Aux[2]
    s2[j+1] <- Aux[3]
    erro <- max(abs(c(b0[j+1]-b0[j],b1[j+1]-b1[j],s2[j+1]-s2[j])))
  }
  
  S <- list()
  S$erro <- erro
  S$Iteracoes <- j
  S$theta <- c(b0[j+1],b1[j+1],s2[j+1])
  S$IF <- IF(b0[j+1],b1[j+1],s2[j+1])
  S$U <- U(b0[j+1],b1[j+1],s2[j+1])
  return(S)
}

# Encontrando as estimativas de MV para theta
theta0 <-  c(10,10,10) # Estimativa inicial
fit <- EF(x,theta0)

# Estimativas:
fit$theta
fit$Iteracoes
# Verificando se os autovalores são negativos
eigen(fit$IF) # Autovalores são estritamente negativos -> ponto de máximo

# E(Y) para cada modelo (modelo não linear) - Escore-Fisher
f1 <- function(x) fit$theta[1]*x^fit$theta[2]

# item e

# Função de log-verossimilhança
log.veross =function(theta){
  b0 = theta[1]
  b1 = theta[2]
  s2 = theta[3]
  return( -(-n/2*log(s2)-1/(2*s2)*sum((y-(b0*x^b1))^2)))
}

# Obtendo estimativas iniciais tentando linearizar o modelo
fit00 = lm(y^1/2 ~ x)
theta00 = c(fit00$coef,summary(fit00)$sigma^2)

# Obtendo as estimativas de MV usando a função optim
fib1 = suppressWarnings(optim(theta00,log.veross, hessian=TRUE))

# E(Y) para cada modelo (modelo não linear) - optim()
f2 = function(x) fib1$par[1]*x^fib1$par[2]

# Exercicio 2

# Função de log-verossimilhança
log.v <- function(theta){
  b0 <- theta[1]
  b1 <- theta[2]
  g0 <- theta[3]
  g1 <- theta[4]
  return( -(-1/2*sum(g0+g1*x)-1/2*sum((y-b0*x^b1)^2/(exp(g0+g1*x)))))
}

# Obtendo estimativas iniciais tentando linearizar o modelo
fit01 = lm(log(y) ~ x)
theta01 = c(37,0.2,-1,0)

# Obtendo as estimativas de MV usando a função optim
fib2 <- suppressWarnings(optim(theta01,log.v, hessian=TRUE))
fib2$par
fib2$convergence
fib2$counts

f3 <- function(x) fib2$par[1]*x^fib2$par[2]

plot(y~x, main="Diagrama de dispersão", pch=19)
curve(f3,add=T, pch=15,col=3)
curve(f1,add=T, pch=15,col=2)

nls(y~b0*x^b1,data=Dados,start=c(b0=3,b1=0.2))

