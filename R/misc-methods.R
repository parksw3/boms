##' @export
print.bomsmodel <- function(x, ...) {
	h <- paste0("d", x$state, "/dt = ", sapply(x$grad, function(x) deparse(x[[3]])))
	cat("\nModel:\n")
	for(i in 1:length(x$grad)) 
		cat(h[i], "\n")
	
	cat("\nObservation: ")
	cat(deparse(x$observation), "\n")
	
	cat("\nInitial values:\n")
	g <- paste0(x$state, "(0) = ", sapply(x$initial, function(x) deparse(x[[3]])))
	for(i in 1:length(g))
		cat(g[i], "\n")
	
	cat("\nFamily:", x$family$family)
	cat("\nParameters:", x$par)
}

print.bomssummary <- function(x, digit = 2, ...) {
	cat(paste0(
		"  Model: ", x$model.name,
		" (Observation: ", deparse(x$model$observation), ") \n"
	))
	
	cat(" Family: ")
	cat(x$family$family, "\n")
	cat("  Links: ")
	cat(paste(paste(names(x$link), c(x$link), sep=":"), collapse = " "), "\n")
	
	cat(paste0(
		"   Data: ", x$data.name, 
		" (Number of observations: ", x$nobs, ") \n"
	))
	
	final_samples <- ceiling((x$iter - x$warmup) / x$thin * x$chains)
	cat(paste0(
		"Samples: ", x$chains, " chains, each with iter = ", x$iter, 
		"; warmup = ", x$warmup, "; thin = ", x$thin, ";\n",
		"         total post-warmup samples = ", final_samples, "\n\n"
	))
	
	if (nrow(x$prior)) {
		cat("Priors: \n")
		print(x$prior, show_df = FALSE)
		cat("\n")
	}
	
	if (nrow(x$cor_pars)) {
		cat("Correlation Structures:\n")
		print_format(x$cor_pars, digits)
		cat("\n")
	}
	if (length(x$random)) {
		cat("Group-Level Effects: \n")
		for (i in seq_along(x$random)) {
			g <- names(x$random)[i]
			cat(paste0("~", g, " (Number of levels: ", x$ngrps[[g]], ") \n"))
			print_format(x$random[[g]], digits)
			cat("\n")
		}
	}
	if (nrow(x$fixed)) {
		cat("Population-Level Effects: \n")
		print_format(x$fixed, digits)
		cat("\n")
	}
	
	if (nrow(x$spec_pars)) {
		cat("Family Specific Parameters: \n")
		print_format(x$spec_pars, digits)
		cat("\n")
	}
	if (length(x$rescor_pars)) {
		cat("Residual Correlations: \n")
		print_format(x$rescor, digits)
		cat("\n")
	}
	cat(paste0("Samples were drawn using ", x$sampler, ". "))
	if (x$algorithm == "sampling") {
		cat(paste0(
			"For each parameter, Eff.Sample \n",
			"is a crude measure of effective sample size, ", 
			"and Rhat is the potential \n",
			"scale reduction factor on split chains ",
			"(at convergence, Rhat = 1)."
		))
	}
	cat("\n")
	
	invisible(x)
}

