library(glmnet)
library(survival)
df <- survival::veteran
df <- df[order(df$time),]
df <- df[!duplicated(df$time),]
delta <- df$status
time <- df$time
So <- Surv(time=time,event=delta)
X <- as.matrix(df[,c('karno','diagtime','age')])
X <- apply(X, 2, function(x){x/sqrt(sum(x^2))})
lambda <- 0


tol = 1e-9          # tolerance
iter = 1000         # number of max iterations
verbose = T

w <- solve(crossprod(X) + diag(lambda, ncol(X))) %*% crossprod(X,delta)
w

tol_curr = 1
J = ncol(X)
a = rep(0, J)
c_ = rep(0, J)
i = 1
ss <- function(j, X, delta, init_beta) {
	N <- NROW(X)
	eta <- as.vector(X[, -j, drop = FALSE] %*% init_beta[-j]) ## BMB: drop()?
	rel_haz <- exp(eta)			# w[i]
	risk_set <- rev(cumsum(rev(rel_haz)))	# W[i]
	P_mat <- outer(rel_haz, risk_set, "/")
	P_mat[upper.tri(P_mat)] <- 0
	w1 <- sum(rel_haz * risk_set - rel_haz^2)
	z1 <- eta + 1/w1 * (delta - sum(P_mat))
	ll <- 1/N * sum(w1 * X[, j, drop = FALSE] * (z1 - X[, -j, drop = FALSE] %*% init_beta[-j]))
	return(ll)
}

c_ <- sapply(1:J, function(j) ss(j, X, delta, w))

soft_thresh <- function(a, b) {
 out = rep(0, length(a))
 out[a >  b] = a[a > b] - b
 out[a < -b] = a[a < -b] + b
 out
}
proxupdate <- function(z, theta){
   prox <- sign(z) * pmax(abs(z) - rep(theta, length(z)), 0)
   return(prox)
}

nloglik <- function(j, X, delta, init_beta) {
	N <- NROW(X)
	eta <- as.vector(X[, -j, drop = FALSE] %*% init_beta[-j]) ## BMB: drop()?
	Xp <- X[, j, drop = FALSE]
	rel_haz <- exp(eta)			# w[i]
	risk_set <- rev(cumsum(rev(rel_haz)))	# W[i]
	P_mat <- outer(rel_haz, risk_set, "/")
	P_mat[upper.tri(P_mat)] <- 0
	nll <- t(Xp) %*% (delta - P_mat %*% delta)
#	nll <- sum(delta * (Xp - P_mat))

#	W <- -P_mat %*% diag(delta) %*% t(P_mat)
#	diag(W) <- diag(P_mat %*% diag(delta) %*% t(1-P_mat))
#	b <- solve(t(Xp)%*%W%*% Xp) %*% nll + init_beta[j]
	return(nll)
#	return(init_beta)
}

nn <- function(j, X, d, b){
	eta <- X[, -j, drop = FALSE] %*% b[-j]
	haz <- as.numeric(exp(eta)) # w[i]
	rsk <- rev(cumsum(rev(haz))) # W[i]
	xx <- X[, j, drop = FALSE]
	P <- outer(haz, rsk, '/')
	P[upper.tri(P)] <- 0
	W <- -P %*% diag(d) %*% t(P)
	diag(W) <- diag(P %*% diag(d) %*% t(1-P))
	b <- solve(t(xx)%*%W%*% xx) %*% t(xx) %*% (d - P%*%d)
	return(b)
}

c_ = sapply(1:J, function(j) nloglik(j, X, delta, w))
c_

#quit()

lambda <- 0
i <- 1
soft <- TRUE
while (tol < tol_curr && i < iter) {
 w_old = w
 a = colSums(X^2)
 l = length(delta) * lambda  # for consistency with glmnet approach
# l = lambda/2  # for consistency with glmnet approach
 c_ = sapply(1:J, function(j) nn(j, X, delta, w_old))
 if (soft) {
	for (j in 1:J) {
	  w[j] = proxupdate(c_[j]/a[j], l/a[j])
	}
 }
 else {
	w = w_old
	w[c_< l & c_ > -l] = 0
 }

 tol_curr = crossprod(w - w_old)
 i = i + 1
 if (verbose && i%%10 == 0) message(i)
}
w
lam <- 0 #1/70
alpha <- 1
mdl.glmnet <- glmnet(x=X,y=So,family='cox',alpha=alpha,lambda=lam,standardize = F)
coef(mdl.glmnet)[,1]



set.seed(8675309)
N = 500
p = 3
X = scale(matrix(rnorm(N*p), ncol=p))
b = c(.5, -.5, .25)
y = scale(X %*% b + rnorm(N, sd=.5))
lambda = .1

w = solve(crossprod(X) + diag(lambda, ncol(X))) %*% crossprod(X,y)
tol_curr = 1
J = ncol(X)
a = rep(0, J)
c_ = rep(0, J)
i = 1

w_old = w
a = colSums(X^2)
l = length(y)*lambda  # for consistency with glmnet approach
c_ = sapply(1:J, function(j)  sum( X[,j] * (y - X[,-j] %*% w_old[-j]) ))
c_



quit()
w <- solve(crossprod(X) + diag(lambda, ncol(X))) %*% crossprod(X,delta)
w[,1] <- 0
lambda <- 1/80
for(step in 1:50){for(j in 1:ncol(X)){
  w_old <- w
  # partial residuals
  eta <- X[, -j, drop = FALSE] %*% w_old[-j]
  xx <- X[, j, drop = FALSE]
  haz <- as.numeric(exp(eta))
  rsk <- rev(cumsum(rev(haz)))
  P <- outer(haz, rsk, '/')
  P[upper.tri(P)] <- 0
  W <- -P %*% diag(d) %*% t(P)
  diag(W) <- diag(P %*% diag(d) %*% t(1-P))
  b <- t(xx) %*% (d - P%*%d)

  # soft-threshold solution
  x2 <- sum(xx^2)
  w[j] <- (abs(b)-lambda * length(d))/(x2 * t(xx)%*%W%*% xx)
  w[j] <- sign(b)*ifelse(w[j]>0,w[j],0)
  # print(w)
  #  print(w)
}}
w
