# Sets up the simple time-invariant survival model
# Used the soft-threshold operator and the negative loglikelihood 
# already defined

pencox <- function(eventvar, X, alpha = 0, lambda = 1
	, gamma = 0.1, standardise = FALSE, maxiter = 5000, tol = 1e-8) {
#	y <- parseFormula(formula = formula, data = data)
#	ynames <- colnames(y)
#	if (!grepl("^Surv", ynames) || NCOL(y) > 1) {
#		stop("formula: must be a survival formula. ?survival")
#	}
	
	formula <- as.formula(paste0("~", colnames(X), collapse = "+"))
	X <- model.matrix(formula, X)[,-1]
	p <- NCOL(X)

	if (standardise) {
		X <- scale(X)
		scale_mu <- attr(X, "scaled:center")
		scale_sd <- attr(X, "scaled:scale")
	} else {
		scale_mu <- rep(0, p)
		scale_sd <- rep(0, p)
		
	}

	# Initialise beta values
	init_beta <- rep(0, p)

	# Coefficient matrix
	beta <- matrix(rep(NA, p*(maxiter+1)), nrow = p, dimnames = list(colnames(X)))

	# Initialise likelihood function
	nll0 <- nloglik(X = X, delta = eventvar, init_beta = init_beta)
	grad0 <- gradient(nll = nll0, init_beta = init_beta, gamma = gamma, lambda = lambda, alpha = alpha)
	beta[, 1] <- proxupdate(init_beta + grad0, lambda * alpha * gamma)

	# Convergence message
	message <- sprintf("Model did not converge after %i iterations, consider increasing maxiter...", maxiter)

	# Proximal update
	for (k in 1:maxiter) {
		nll <- nloglik(X = X, delta = eventvar, init_beta = beta[,k])
		grad <- gradient(nll = nll
			, init_beta = beta[, k]
			, gamma = gamma
			, lambda = lambda
			, alpha = alpha
		)
		beta[, k+1] <- proxupdate(beta[,k] + grad, lambda * alpha * gamma) 

		# Check for convergence
		concheck <- sqrt(sum((beta[,k+1] - beta[, k])^2))/gamma
		if (concheck < tol) {
     		beta <- beta[, -((k+2):ncol(beta))]

			message <- sprintf("Model converged after %i iterations", (k+1))
			break
		} 
	}
	print(message)
	beta_hat <- as.matrix(beta[,ncol(beta)])
	colnames(beta_hat) <- "s0"
	return("beta_hat" = beta_hat)
}
