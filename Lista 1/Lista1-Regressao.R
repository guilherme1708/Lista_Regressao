library(ggplot2)
library(MASS)

# Exercício 2

x <- c(17.33,-90.91,18.49,25.03,19.88,18.51,12.39,11.24,11.78,15.15,
       12.24,-238.41,15.19,24.51,11.82,18.16,15.73,11.26,7.41,19.46,
       25.58,-22.27,16.72,17.37,17.42,13.37,12.15,13.75,24.93,34.49,
       21.66,-7.00,16.21,42.07,-2.70,14.91,16.36,14.57,15.11,12.63,
       2.96,13.24,18.68,19.22,30.43,12.87,11.40,-68.39,14.43,13.28,
       3.24,12.63,11.92,22.08,19.47,13.84,20.68,-10.52,12.78,19.03,
       16.65,16.46,2.65,13.55,15.61,59.20,11.46,17.26,18.17,21.50,
       20.31,9.63,13.84,17.94,16.35,18.40,-16.15,17.27,255.83,20.74,
       -1.08,13.67,-5.82,-0.75,21.05,16.70,23.67,11.82,6.23,17.61,
       5.83,19.83,11.14,18.06,12.11,19.83,16.20,11.03,13.99,8.70)

theta0 <- c(median(x),1) # Estimativa inicial

MV = function(x,theta0,epsilon=10^(-6), int=1000){ # Dados, estimativa inicial, erro, iterações
  n <- length(x)
  
  U <- function(t1,t2){
    U1 <- 2*sum((x-t1)/(t2^2+(x-t1)^2))
    U2 <- n/t2 - 2*sum(t2/(t2^2 + (x-t1)^2))
    return (c(U1,U2))
  }
  
  H <- function(t1,t2){
    H11 <- 2*sum(((x-t1)^2-t2^2)/((x-t1)^2+t2^2)^2)
    H12 <- 4*t2*sum((x-t1)/((x-t1)^2+t2^2)^2)
    H22 <- sum(((x-t1)^2-t2^2)/((x-t1)^2+t2^2)^2) - n/t2^2
    M <- matrix(c(H11,H12,H12,H22),2,2)
    return (M)
  }
  
  erro <- 10
  j <- 0
  t1 <- numeric ()
  t2 <- numeric()
  t1[1]<-theta0[1]
  t2[1]<-theta0[2]
  
  while(erro > epsilon & j < int){
    j <- j+1
    Aux <- c(t1[j],t2[j]) - solve(H(t1[j],t2[j]))%*%U(t1[j],t2[j])
    t1[j+1] <- Aux[1]
    t2[j+1] <- Aux[2]
    erro <- max(abs(c(t1[j+1]-t1[j],t2[j+1]-t2[j])))
    print(erro)
  }
  
  S <- list()
  S$erro <- erro
  S$Iteracoes <- j
  S$theta <- c(t1[j+1],t2[j+1])
  S$H <- H(t1[j+1],t2[j+1])
  S$U <- U(t1[j+1],t2[j+1])
  return(S)
}

#Encontrando as estimativas de MV para o exercício
fit <- MV(x,theta0)

# Estimativas:
fit$theta

# Verificando se os autovalores são negativos
eigen(fit$H) # Autovalores são estritamente negativos -> ponto de máximo

# Estimativas
t1.hat <- fit$theta[1]
t2.hat <- fit$theta[2]

# Gráficos QQ
x <- data.frame(x)
ggplot(x, aes(sample = x)) + stat_qq(distribution = qcauchy) + stat_qq_line(distribution = qcauchy) + ggtitle("QQ-Plot Cauchy")
ggplot(x, aes(sample = x)) + stat_qq(distribution = qnorm) + stat_qq_line(distribution = qnorm) + ggtitle("QQ-Plot Normal")

# Histograma com densidade
ggplot(x, aes(x=x)) + geom_histogram(bins = 30, aes(y=..density..)) + labs(y = "Frequência", x = "Dados Observados") + ggtitle("Histograma da amostra") + stat_function(fun = function (y) dcauchy(y,t1.hat,t2.hat))


# Exercício 3

X <- c(21, 24, 32, 47, 50, 68, 74, 62, 50, 41, 30) # Variável explicativa  
Y <- c(185.79, 214.47, 288.03, 424.84, 539.03, 621.55, 675.06, 562.03, 452.93, 369.95, 273.98) #  Variável resposta
X <- X - median(X) # correção dos dados para a interpretação

df <- data.frame(Y,X)

ggplot(df, aes(x=X, y=Y)) + geom_point() + geom_smooth(method=lm, se=F) + labs(y = "Uso de Vapor (lb)", x = "Temperatura (ºF)") +
   ggtitle("Diagrama de Dispersão")

cor(X,Y) # correlação X,Y

reg_lin <- lm(Y ~ X) # Regressão linear
summary(reg_lin) 
