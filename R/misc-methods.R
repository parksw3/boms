##' @export
print.bomsmodel <- function(x, ...) {
	h <- paste0("d", x$state, "/dt = ", sapply(x$grad, function(x) deparse(x[[3]])))
	cat("\nModel:\n")
	for(i in 1:length(x$grad)) 
		cat(h[i], "\n")
	
	cat("\nObservations:\n")
	cat(deparse(x$observation), "\n")
	
	cat("\nInitial values:\n")
	g <- paste0(x$state, "(0) = ", sapply(x$initial, function(x) deparse(x[[3]])))
	for(i in 1:length(g))
		cat(g[i], "\n")
	
	cat("\nFamily:", x$family$family)
	cat("\nParameters:", x$par)
}

print.bomssummary <- function(x, digit = 2, ...) {
	cat(" Family: ")
	cat(summarise_families(x$formula), "\n")
	cat("  Links: ")
	cat(summarise_links(x$formula, wsp = 9), "\n")
	cat("Formula: ")
	
	invisible(x)
}