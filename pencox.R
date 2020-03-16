# Sets up the simple time-invariant survival model
# Used the soft-threshold operator and the negative loglikelihood 
# already defined

pencox <- function(eventvar, X, alpha = 0, lambda = 1
	, gamma = 0.1, standardise = FALSE, maxiter = 500) {
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
		scale_sd <- attr(X,'scaled:scale')
	} else {
		scale_mu <- rep(0, p)
		scale_sd <- rep(0, p)
		
	}

	# Estimate betas
	init_beta <- as.matrix(rep(0, p))
	
	# Proximal update
	for (k in 1:maxiter) {
		nll <- nloglik(X = X, delta = eventvar, init_beta = init_beta)
		grad <- gradient(nll = nll
			, init_beta = init_beta
			, gamma = gamma
			, lambda = lambda
			, alpha = alpha
		)
		init_beta <- proxupdate(init_beta + grad, lambda * alpha * gamma) 
	}

	return(init_beta)
}
