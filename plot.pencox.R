# Plot coefficient estimates as a function of lambda

plot.pencox <- function(cv.object){
	plot.df <- data.frame(do.call(cbind, cv.object$cv.list))
	beta.df <- plot.df[, grepl("s0", names(plot.df))]
	lambda.df <- plot.df[, grepl("lambda", names(plot.df))][1,]
	plot(0, 0, type = "n", ylim = c(min(beta.df), max(beta.df)), cex.main = 0.8
		, ylab=expression(hat(beta)[lambda]), xlim = log(c(min(lambda.df), max(lambda.df)))
		, xlab=expression(paste('log(',lambda,')'))
		, main = 'Proximal Gradient Cox PH coefficients'
	)
	
	for (i in 1:nrow(beta.df)) { lines(log(lambda.df[1,]), beta.df[i, ], col=i) }
}
