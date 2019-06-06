library(MPV)
library(MASS)
library(dplyr)

# exercício 1

data("cement")
attach(cement)

# item a

model <- lm(y~x1+x2+x3+x4)
p <- 4
Influence <- lm.influence(model)
sigma.hat <- summary(model)$sigma
t <- resid(model)/(sigma.hat*sqrt(1-Influence$hat))
t.M <- resid(model)/(Influence$sigma*sqrt(1-Influence$hat))
D <- t^2*Influence$hat/((p+1)*(1- Influence$hat))

# item b

## i. de quantil dos resíduos contras os quantis esperados de uma normal padrão;
qqnorm(resid(model))
qqline(resid(model))

## ii. de resíduos contra os índices das unidades amostrais;
plot(rep(1:13),resid(model), xlab="Unidades Amostrais", ylab="Resíduos", xlim=c(0.5,13.5), main="Resíduos x Unidades amostrais",pch=19)
abline(h=-3,lty=3)
abline(h=3,lty=3)

## iii. dos resíduos contra µ chapeu, i in {1,...,13};
plot(model$fitted.values,resid(model), pch=19, main="Resíduos x Valores Preditos", ylab="Resíduos", xlab="Valores Preditos")
abline(h=0,lty=3)

## iv. de resíduos contra cada uma das covariáveis;
par(mfrow=c(2,2))
plot(x1,resid(model),pch=19, ylab="Resíduos", xlab=expression(x[1]), main=expression(paste("Resíduos x Covariável ",x[1])))
abline(h=-3,lty=3)
abline(h=3,lty=3)
plot(x2,resid(model),pch=19, ylab="Resíduos", main=expression(paste("Resíduos x Covariável ",x[2])))
abline(h=-3,lty=3)
abline(h=3,lty=3)
plot(x3,resid(model),pch=19, ylab="Resíduos", main=expression(paste("Resíduos x Covariável ",x[3])))
abline(h=-3,lty=3)
abline(h=3,lty=3)
plot(x4,resid(model),pch=19, ylab="Resíduos", main=expression(paste("Resíduos x Covariável ",x[4])))
abline(h=-3,lty=3)
abline(h=3,lty=3)

## v. das distâncias de Cook;

plot(D,ylab=expression(D[i]), main="Distância de Cook", pch=19)
points(8,D[8],pch=19,col=2)
abline(h=3*mean(D),lty=3)

## vi. de alavancagem;

plot(Influence$hat, ylab=expression(h[ii]), main="Pontos de alavanca")
abline(h=2*(p+1)/n,lty=3)

## vii. de envelope simulado.

envelope = function(fit0){
  Influence = lm.influence(fit0)
  t.M = resid(fit0)/(Influence$sigma*sqrt(1-Influence$hat))
  n = length(t.M)
  X  = model.matrix(fit0)
  H = X%*%solve(t(X)%*%X)%*%t(X) 
  I = diag(n)
  yy = matrix(rnorm(n*100),n,100)
  A = matrix(0,n,100)
  LI = LS = numeric(n)
  A = (I - H)%*%yy/sqrt(diag(I - H))
  A = apply(A,2,sort)
  B = t(apply(A,1,sort))
  LI = (B[,2]+B[,3])/2
  LS = (B[,97]+B[,98])/2
  med= apply(B,1,mean)
  aux = range(t.M,LS,LI)
  #par(pty="s")
  qqnorm(t.M,xlab="Percentil da N(0,1)", ylab="Ressíduo Studentizado modificado", ylim=aux,pch=16, main="")
  par(new=TRUE)
  qqnorm(LI,axes=F,xlab="",ylab="",type="l",lty=1,ylim=aux)
  par(new=TRUE)
  qqnorm(LS,axes=F,xlab="",ylab="", type="l", lty=1,ylim=aux)
  par(new=TRUE)
  qqnorm(med,axes=F,xlab="",ylab="",type="l",ylim=aux,lty=2)
}

envelope(model)

# item c

model1 <- lm(x1 ~ x2 + x3 + x4, data=cement)
FIV.1 <- 1/(1-summary(model1)$r.squared)

model2 <- lm(x2 ~ x1 + x3 + x4, data=cement)
FIV.2 <- 1/(1-summary(model2)$r.squared)

model3 <- lm(x3 ~ x2 + x1 + x4, data=cement)
FIV.3 <- 1/(1-summary(model3)$r.squared)

model4 <- lm(x4 ~ x1 + x2 + x3, data=cement)
FIV.4 <- 1/(1-summary(model4)$r.squared)

# item d

fit0 <- lm(y ~ x1 + x2 + x3)

model5 <- lm(x1 ~ x2 + x3, data=cement)
FIV.5 <- 1/(1-summary(model5)$r.squared)

model6 <- lm(x2 ~ x1 + x3, data=cement)
FIV.6 <- 1/(1-summary(model6)$r.squared)

model7 <- lm(x3 ~ x2 + x1, data=cement)
FIV.7 <- 1/(1-summary(model7)$r.squared)

# item e

# Exercício 2

Y <- c(20,60,20,60,30,30,10,60,50,80,40,90,140,40,100,
  150,80,130,60,30,150,80,100,200,120,100,110,80,160,150)

X <- c(1,3,4,6,8,9,11,13,15,16,18,20,21,23,25,
       26,28,30,31,33,35,36,38,40,42,43,45,47,48,50)

# item a

fit <- lm(Y~X)

# item b
n <- length(X)
p <- length(X)-1
Influence <- lm.influence(fit)
sigma.hat <- summary(fit)$sigma
d <- resid(fit)/sigma.hat
t <- resid(fit)/(sigma.hat*sqrt(1-Influence$hat))
t.M <- resid(fit)/(Influence$sigma*sqrt(1-Influence$hat))
D <- t^2*Influence$hat/(3*(1- Influence$hat))

plot(Y~X, main="Diagrama de dispersão", xlab="Tempo", ylab="Número de Bactérias")
abline(fit)

par(mfrow=c(2,2))
plot(Influence$hat, ylab=expression(h[ii]), main="Pontos de alavanca")
abline(h=2*(p+1)/n,lty=3)
plot(t.M,ylab="Resíduo studentizado modificado", main="Resíduos studentizados Modificados",ylim=c(-4,4),type="l")
abline(h=-3,lty=3)
abline(h=3,lty=3)

plot(t.M ~ predict(fit),ylab="Resíduo studentizado modificado", main="Resíduos studentizados Modificados",ylim=c(-4,4),xlab=expression(hat(mu)[i]))
abline(h=-3,lty=3)
abline(h=3,lty=3)

plot(D,ylab=expression(D[i]), main="Distância de Cook")
abline(h=3*mean(D),lty=3)

qqnorm(resid(fit))
qqline(resid(fit))

# item c

boxcox(fit) # 1/2

Y <- c(20,60,20,60,30,30,10,60,50,80,40,90,140,40,100,
       150,80,130,60,30,150,80,100,200,120,100,110,80,160,150)
X <- c(1,3,4,6,8,9,11,13,15,16,18,20,21,23,25,
       26,28,30,31,33,35,36,38,40,42,43,45,47,48,50)

Y1 <- log(Y)

fit1 <- lm(Y1~X)

Influence <- lm.influence(fit1)
sigma.hat <- summary(fit1)$sigma
d <- resid(fit1)/sigma.hat
t <- resid(fit1)/(sigma.hat*sqrt(1-Influence$hat))
t.M <- resid(fit1)/(Influence$sigma*sqrt(1-Influence$hat))
D <- t^2*Influence$hat/(2*(1- Influence$hat))

par(mfrow=c(2,2))

plot(Influence$hat, ylab=expression(h[ii]), main="Pontos de alavanca")
abline(h=2*(p)/n,lty=3)
plot(t.M,ylab="Resíduo studentizado modificado", main="Resíduos studentizados Modificados",ylim=c(-4,4),type="l")
abline(h=-3,lty=3)
abline(h=3,lty=3)

plot(t.M ~ predict(fit1),ylab="Resíduo studentizado modificado", main="Resíduos studentizados Modificados",ylim=c(-4,4),xlab=expression(hat(mu)[i]))
abline(h=-3,lty=3)
abline(h=3,lty=3)

plot(D,ylab=expression(D[i]), main="Distância de Cook")
abline(h=3*mean(D),lty=3)

envelope(fit1)

# Exercício 5

y <- c(28.29,45.58,33.49,40.27,49.56,40.00,21.83,33.96,
       57.07,49.90,39.45,36.90,42.07,31.83,55.11,48.29,20.58,43.49,60.27,34.56)

x <- c(6,10,7,7,9,15,5,7,9,9,8,8,8,8,11,10,5,9,11,6)

# item a

Y <- cbind(y) # matriz de variável resposta
X <- cbind(1,x) # matriz de dados (fixados)
n <- length(Y)  # tamanho da amostra
V <- solve(t(X)%*%X) # (X'X)-¹
D <- V%*%t(X) # (X'X)-¹*X'

beta.MQ <- D%*%Y # matriz Beta0 e beta1 estimados
sigma2.hat <- as.numeric(t(Y)%*%(diag(n) - X%*%D)%*%Y/(n-2)) # variancia estimada não viciada
sigma.hat <- sqrt(sigma2.hat) # desvio padrão estimada não viciada

fit0 <- lm(y~x)

F <- function(C,c,beta,sigma2,n){ # Função para calcular a estatística F e o p-valor
  k <- nrow(C)
  F.obs <- as.numeric(t(C%*%beta - c)%*% solve( C%*% solve(t(X)%*%X) %*% t(C)) %*% (C%*%beta - c)/(k*sigma2))
  p.valor = 1 - pf(F.obs,k,n-2)
  return(c(F.obs,p.valor))
}

#H_0: B1 = 5
C <- matrix(c(0,1),1,2) # Matriz C
d <- 5 # Vetor d
F(C, d, beta.MQ, sigma2.hat, length(Y)) # estatística F e o p-valor


# item b

Influence <- lm.influence(fit0)
p <- 2
sigma.hat <- summary(fit0)$sigma
d <- resid(fit0)/sigma.hat
t <- resid(fit0)/(sigma.hat*sqrt(1-Influence$hat))
t.M <- resid(fit0)/(Influence$sigma*sqrt(1-Influence$hat))
D <- t^2*Influence$hat/((p)*(1- Influence$hat))

plot(D,ylab=expression(D[i]), main="Distância de Cook")
points(6,D[6],pch=19,col=2)
abline(h=3*mean(D),lty=3)

# item c

y1 <- c(28.29,45.58,33.49,40.27,49.56,21.83,33.96,
       57.07,49.90,39.45,36.90,42.07,31.83,55.11,48.29,20.58,43.49,60.27,34.56)

x1 <- c(6,10,7,7,9,5,7,9,9,8,8,8,8,11,10,5,9,11,6)

Y <- cbind(y1) # matriz de variável resposta
X <- cbind(1,x1) # matriz de dados (fixados)
n <- length(Y)  # tamanho da amostra
V <- solve(t(X)%*%X) # (X'X)-¹
D <- V%*%t(X) # (X'X)-¹*X'

beta.MQ <- D%*%Y # matriz Beta0 e beta1 estimados
sigma2.hat <- as.numeric(t(Y)%*%(diag(n) - X%*%D)%*%Y/(n-2)) # variancia estimada não viciada
sigma.hat <- sqrt(sigma2.hat) # desvio padrão estimada não viciada

fit1 <- lm(y1~x1)


Influence <- lm.influence(fit1)
p <- 2
sigma.hat <- summary(fit1)$sigma
d <- resid(fit1)/sigma.hat
t <- resid(fit1)/(sigma.hat*sqrt(1-Influence$hat))
t.M <- resid(fit1)/(Influence$sigma*sqrt(1-Influence$hat))
D <- t^2*Influence$hat/((p)*(1- Influence$hat))

plot(D,ylab=expression(D[i]), main="Distância de Cook")
abline(h=3*mean(D),lty=3)

#H_0: B1 = 5
C <- matrix(c(0,1),1,2) # Matriz C
d <- 5 # Vetor d
F1 <- F(C, d, beta.MQ, sigma2.hat, length(Y)) # estatística F e o p-valor

