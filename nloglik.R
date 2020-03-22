# Define negative log likelihood

nloglik <- function(X, delta, init_beta) {
	N <- NROW(X)
	eta <- as.vector(X %*% init_beta)
	rel_haz <- as.numeric(exp(eta))			# w[i]
	risk_set <- rev(cumsum(rev(rel_haz)))	# W[i]
	P_mat <- outer(rel_haz, risk_set, "/")
	P_mat[upper.tri(P_mat)] <- 0
	nll <- 1/N * (t(X) %*% (delta - P_mat %*% delta))
	return(nll)
}
 
