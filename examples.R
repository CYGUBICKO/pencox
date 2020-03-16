# Using veteran dataset and compare result with glmnet

library(dplyr)
library(glmnet)
library(survival)

source("proxupdate.R")
source("nloglik.R")
source("gradient.R")
source("pencox.R")

df <- veteran

## Sort with time
df <- df[order(df$time),]

## Current, ties are not implemented, remove ties
df <- df[!duplicated(df$time),]

eventvar <- df$status
X <- (df
	%>% select(-c("trt", "time", "status"))
	%>% data.frame()
)
pencox_res <- pencox(eventvar, X, gamma = 0.1, lambda = 0.5, maxiter = 500, standardise = TRUE)
pencox_res

## glmnet

timevar <- df$time
y <- Surv(time = timevar, event = eventvar)
formula <- as.formula(paste0("~", colnames(X), collapse = "+"))
X <- scale(model.matrix(formula, X)[,-1])
glmnet_res <- coef(glmnet(x = X ,y = y, family = 'cox', alpha = 0, lambda = 0.5,standardize = FALSE))
glmnet_res

## Cox-PH
d1 <- (X
	%>% data.frame()
	%>% mutate(timevar = timevar, eventvar = eventvar)
)
coxph_res <- coxph(Surv(time = timevar, event = eventvar) ~., data = d1)
summary(coxph_res)
