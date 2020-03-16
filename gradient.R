# Compute the penalisation term of proximal GD

gradient <- function(nll, init_beta, gamma, lambda, alpha) {
	res <- gamma*nll - gamma*lambda*(1 - alpha)*init_beta
	return(res)
}

