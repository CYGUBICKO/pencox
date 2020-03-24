## Based on Breslow-estimator

cumhaz <- function(fit, centered = TRUE){
	beta.hat <- fit$beta.hat
	delta <- fit$delta
	X <- fit$Xmatrix
	eta <- drop(X %*% beta.hat)
	rel_haz <- exp(eta)
	risk_set <- rev(cumsum(rev(rel_haz)))
	h0 <- delta/risk_set
	hazard <- cumsum(h0)
	surv.est <- exp(-hazard)
	if (centered){
		beta.mean <- apply(X, 2, mean)
		surv.est <- surv.est^(exp(drop(X %*% beta.mean)))
		hazard <- -log(surv.est)
	}
	## Index is a place holder for time, future revisions will add time var
	H0 <- data.frame(index = 1:length(hazard), cumhazard = hazard, survestimate = surv.est)
	return(H0 = H0)
}
