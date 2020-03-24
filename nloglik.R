# Define negative log likelihood

nloglik <- function(X, delta, init_beta) {
	N <- NROW(X)
	eta <- as.vector(X %*% init_beta) ## BMB: drop()?
	rel_haz <- exp(eta)			# w[i]
	risk_set <- rev(cumsum(rev(rel_haz)))	# W[i]
        ## BMB: this is inefficient since it
        ##      computes for all {i,j} (not i<j)
        ##   I don't know of any easy, terse method
        ##    that computes only the upper triangle
        ##    if this is a bottleneck, could rewrite in Rcpp ...
		  ## Steve: For loop over all delta but I guess it would be much slower?
	P_mat <- outer(rel_haz, risk_set, "/")
	P_mat[upper.tri(P_mat)] <- 0
	nll <- 1/N * (t(X) %*% (delta - P_mat %*% delta))
	return(nll)
}
