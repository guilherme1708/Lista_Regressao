# Lista 2 MAE328

library(readxl)
library(ggplot2)

# Exercício 5
x <- c(18.10,	10.71,	13.43,	11.04,	22.81,	18.71,	6.52 ,	15.78,	
       19.59,	21.18,	15.82,	15.43,	20.61,	16.70,	18.52,	12.75,	
       19.61,	14.80,	8.87 ,	21.82,	21.09,	26.47,	19.18,	23.31,	
       15.9 ,	11.51,	7.09 ,	12.95,	21.99,	5.75 , 	15.88,	14.52,	
       17.93,	5.19 ,  12.57,	19.78,	18.96,	14.24,	13.38,	19.38,	
       2.47 ,	12.12,	5.46 ,	22.37,	22.51,	20.90,	9.92 ,	19.97,	
       18.10,	20.50,	12.98,	21.18,	20.34,	30.71,	6.27 ,	15.47,	
       18.08,	24.60,	7.26 ,	23.03,	23.67,	24.29,	11.70,	15.21,	
       10.00,	10.93,	12.65,	15.70,	20.04,	21.97,	17.63,	10.72,	
       19.71,	16.45,	15.04)

y <- c(23.38,	38.07,	29.79,	31.47,	18.44,	32.69,	36.14,	26.18,
       20.58,	18.69,	26.09,	9.27 ,	19.09,	24.20,	20.82,	29.33,	
       19.94,	28.82,	35.16,	21.52,	17.11,	16.17,	22.28,	15.77,	
       26.54,	33.40,	50.76,	30.41,	16.43,	35.83,	29.19,	28.16,	
       23.34,	43.88,	31.64,	18.41,	21.59,	29.18,	31.34,	26.06,	
       46.42,	32.52,	63.93,	15.96,	16.88,	19.55,	35.06,	18.90,	
       24.13,	17.46,	30.12,	18.86,	9.10 ,	2.41 ,  39.80,	24.73,
       19.67,	13.06,	24.18,	14.01,	15.27,	13.23,	31.12,	27.41,
       70.36,	31.15,	32.10,	25.87,	20.60,	17.28,	24.25,	34.70,
       8.55 ,  24.48,	27.06)

n <- length(y)

MV = function(x,theta0,epsilon=10^(-6), it=1000){ # Dados, estimativa inicial, erro, iterações
  
  U <- function(t1,t2){
    U1 <- 2*sum((-t1-t2*x+y)/(((-t2*x-t1+y)^2)+1))
    U2 <- 2*sum((x*(-t2*x-t1+y))/(((-t2*x-t1+y)^2)+1))
    return (c(U1,U2))
  }
  
  H <- function(t1,t2){
    H11 <- 2*sum(((-t1-t2*x+y)^2-1)/((-t1-t2*x+y)^2+1)^2)
    H12 <- 2*sum(x*((-t1-t2*x+y)^2-1)/((-t1-t2*x+y)^2+1)^2)
    H22 <- 2*sum(x^2*((-t1-t2*x+y)^2-1)/((-t1-t2*x+y)^2+1)^2)
    M <- matrix(c(H11,H12,H12,H22),2,2)
    return (M)
  }
  
  erro <- 10
  j <- 0
  t1 <- numeric()
  t2 <- numeric()
  t1[1]<-theta0[1]
  t2[1]<-theta0[2]
  
  while(erro > epsilon & j < it){
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
theta0 <- c(52.45,-1.63)
fit <- MV(x,theta0)

# Estimativas:
fit$theta

# Verificando se os autovalores são negativos
eigen(fit$H) # Autovalores são estritamente negativos -> ponto de máximo

# Estimativas
t1.hat <- fit$theta[1]
t2.hat <- fit$theta[2]

# exercício d
like <- function(b,x,y){ # Função log-verossimilhança
  b0 <- b[1]
  b1 <- b[2]
  f <- -n*log(pi)-sum(log(1+(y-b0-b1*x)^2))
  return(-f)
}

cau.fit <- optim(par = c(52.45,-1.63), # função optim()
                        fn = like, 
                        y=y,
                        x=x,
                        method = "CG",
                        hessian = T)

e.optim <- cau.fit$par # estimativas

dados <- data.frame(x,y)
ggplot(dados, aes(x=x, y=y)) + geom_point() + labs(y ="variável resposta", x ="Variável explicativa") +
  ggtitle("Diagrama de Dispersão") +
  stat_function(fun = function (z) e.optim[1] + e.optim[2]*z)

# [48.274476  , -1.39879787]

# Exercício 6

setwd("~/Downloads") # fixando o local 
dados <- read_excel('slr06.xlsx') # import dos dados

Y <- cbind(dados$Y) # matriz de variável resposta
X <- cbind(1,dados$X) # matriz de dados (fixados)
n <- length(Y)  # tamanho da amostra
V <- solve(t(X)%*%X) # (X'X)-¹
D <- V%*%t(X) # (X'X)-¹*X'

beta.MQ <- D%*%Y # matriz Beta0 e beta1 estimados
sigma2.hat <- as.numeric(t(Y)%*%(diag(n) - X%*%D)%*%Y/(n-2)) # variancia estimada não viciada
sigma.hat <- sqrt(sigma2.hat) # desvio padrão estimada não viciada

Matrix.cov <- sigma2.hat*V # matriz de covariancias 

x<-dados$X
reg.lin <- lm(Y~x) # modelo de regressão linear
summary(reg.lin) # resumo da regressõa

ggplot(dados, aes(x=X, y=Y)) + geom_point() + geom_smooth(method=lm, se=F) + labs(y = "pagamento total em milhares (SEK).", x = "Número de sinistros") +
  ggtitle("Diagrama de Dispersão") 

predict(reg.lin, data.frame(x=c(100)), se.fit = TRUE)$fit # valor predito

# Exercício 7

N <- 10000 # Número de simulações
M <- matrix(NA,N,3) # matrix com as simulções dos parametros

for (i in 1:N){
  erro <- rnorm(nrow(X),0,sigma.hat)
  yy <- cbind(beta.MQ[1] + beta.MQ[2]*x + erro) # vetor coluna de Y simulados dado que os o vetor x está fixo
  beta.S <- D%*%yy # vetor coluna de simuações de beta0 e beta1
  sigma2.hatS <- as.numeric(t(yy)%*%(diag(n) - X%*%D)%*%yy/(n-2)) # vetor coluna de simuações de sigma² (ñ viciado)
  M[i,] <- c(beta.S,sigma2.hatS) # matriz com as simulcões
}


# Histograma com densidade
M.df <- data.frame(M) # data frame para o ggplot
ggplot(M.df, aes(x=X1)) + geom_histogram(bins = 30, aes(y=..density..)) +
  labs(y = "Frequência", x = "Dados Observados") +
  ggtitle("Histograma da amostra") + 
  stat_function(fun = function (y) dnorm(y,beta.MQ[1],sqrt(Matrix.cov[1,1])))

ggplot(M.df, aes(x=X2)) + geom_histogram(bins = 30, aes(y=..density..)) + 
  labs(y = "Frequência", x = "Dados Observados") + 
  ggtitle("Histograma da amostra") +
  stat_function(fun = function (y) dnorm(y,beta.MQ[2],sqrt(Matrix.cov[2,2])))

ggplot(M.df, aes(x=((n-2)*X3)/sigma.hat^2)) + geom_histogram(bins = 40, aes(y=..density..)) + 
  labs(y = "Densidade de frequência", x = "Dados Observados") + 
  ggtitle("Histograma da amostra") +
  stat_function(fun = function (y) dchisq(y,n-2))