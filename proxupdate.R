# Defining soft-thresholding operator
# for proximal gradient descent update
set.seed(7778)
proxupdate <- function(z, theta){
	prox <- sign(z) * pmax(abs(z) - rep(theta, length(z)), 0)
	return(prox)
}

