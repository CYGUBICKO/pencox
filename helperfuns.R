# Get formula response variable
getResponse <- function(formula) {
	tt <- terms(formula)
	vars <- as.character(attr(tt, "variables"))[-1]
	response <- attr(tt, "response")
	vars[response]
}

# Parse formula and return response variable
parseFormula <- function(formula, data, env = parent.frame()) {
	f <- as.formula(formula)
	t <- terms(f, data = data)

	## Get dependent var(s)
	response <- data.frame(eval(f[[2]], envir = data, enclos = env))
	colnames(response) <- deparse(f[[2]])
	return(response)
}
