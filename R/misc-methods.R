##' @export
print.bomsmodel <- function(x, ...) {
	h <- paste0("d", x$state, "/dt = ", sapply(x$model, function(x) deparse(x[[3]])))
	cat("\nModel:\n")
	for(i in 1:length(x$model)) 
		cat(h[i], "\n")
	
	cat("\nObservations:\n")
	cat(deparse(x$observation), "\n")
	
	cat("\nInitial values:\n")
	g <- paste0(x$state, "(0) = ", sapply(x$initial, function(x) deparse(x[[3]])))
	for(i in 1:length(g))
		cat(g[i], "\n")
	
	cat("\nFamily:", x$family)
	cat("\nParameters:", x$par)
}
