# Lista 4 - Regressão I MAE0328

# Exercício 5

library(readxl)

# item a
dados <- read_xls("data-table-B6.XLS")
attach(dados)

model <- lm(y~x1+x4,data = dados)

# item b

model$fitted.values
max(model$fitted.values)


# item c
cor(x1,x4)

Y <- cbind(y)
X <- cbind(1,x1,x4)
V <- (t(X)%*%X) # (X'X)
