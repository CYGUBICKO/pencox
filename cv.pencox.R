# Setup for performing CV. Currently, minimizing sum(b_k+1 - b_k)^(0.5)

cv.pencox <- function(eventvar, X, alpha = 0, lambda = 1
	, gamma = 0.1, standardise = FALSE, maxiter = 5000, tol = 1e-8) {
	
	cv.result <- foreach(l = 1:length(lambda)) %dopar% {
		model <- pencox(eventvar, X, alpha, lambda[l], gamma, standardise, maxiter, tol)
		result <- data.frame(beta.hat = model[["beta.hat"]]
			, lambda = model[["lambda"]]
			, alpha = model[["alpha"]]
			, min.deviance = model[["min.deviance"]]
		)
	}
	# yet to figure out how to return min lambda and the rest
	return(list(cv.list = cv.result))
}
