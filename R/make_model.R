##' Set up link functions for ode/loglik parameters
##'
##' @param link named list or vector of strings specifying link functions
##' @param model model.ode object
##' @param loglik loglik.ode object
##' @seealso \code{\link{make.link}}
##' @return list of strings specifying link functions
set_link <- function(link, modelpar, data) {
	link_default <- c(rep("log", length(modelpar)))
	names(link_default) <- modelpar
	
	link_default[names(link_default) %in% colnames(data)] <- "identity"
	
	if (!missing(link)) link_default[names(link)] <- link
	
	link_default
}

##' Apply link functions to parameters
##'
##' @param par vector of parameter values
##' @param linklist list containing \code{linkfun}, \code{linkinv}, and \code{mu.eta} for each link
##' @param type string specifying which function should be applied
##' @seealso \code{\link{make.link}}
apply_link <- function(par, linklist, type=c("linkfun", "linkinv", "mu.eta")) {
	type <- match.arg(type)
	
	if (type=="linkinv" || type=="mu.eta") {
		filter <- match(names(par), names(linklist$linkfun))
	} else {
		filter <- match(names(par), names(linklist$linkinv))
	}
	
	ff <- linklist[[type]][filter]
	pp <- unlist(Map(function(x, fun) fun(x), x=par, fun=ff))
	names(pp) <- names(ff)
	
	pp
}

stan_link <- function(link) {
	switch (link,
			logit="inv_logit",
			log="exp"
	)
}

stan_link_deriv <- function(link) {
	switch (link,
			logit="logistic_lpdf",
			log=""
			)
}

##' @export
make_model <- function(grad,
					   observation,
					   initial,
					   par,
					   data, tcol,
					   family,
					   link,
					   solver=c("rk45", "bdf")) {
	## TODO: add more warnings and errors
	if (any(sapply(grad, class) != "formula"))
		stop("grad must be a list of formulas")
	
	if (any(sapply(initial, class) != "formula"))
		stop("initial must be a list of formulas")
	
	if (class(observation) != "formula")
		stop("observation must be a formulas")
	
	solver <- match.arg(solver)
	
	state <- sapply(grad, function(x) deparse(x[[2]]))
	nstate <- length(state)
	
	## substitute state in gradient and initial
	tstate <- paste0("y[", 1:length(state),"]")
	
	## substitute states in observation
	tobs <- paste0("y_hat[", 1:length(state),"]")
	
	## substitute parameters
	tpar <- paste0("params[", 1:length(par), "]")
	
	state_transforms <- sfun(state, tstate)
	
	obs_transforms <- sfun(state, tobs)
	
	param_transforms <- sfun(par, tpar)
	
	gtextlist <- itextlist <- vector('list', length(grad))
	
	for (i in 1:length(grad)) {
		tt <- subst(subst(grad[[i]][[3]], state_transforms), param_transforms)
		
		gtextlist[[i]] <- paste0("    dydt[", i, "] = ", do.call(paste0, as.list(deparse(tt))), "; \n")
	}
	
	for (i in 1:length(initial)) {
		ii <- subst(initial[[i]][[3]], param_transforms)
		
		itextlist[[i]] <- paste0("    y_ini[", i, "] = ", do.call(paste0, as.list(deparse(ii))), "; \n")
	}
	
	oo <- subst(observation[[3]], obs_transforms)
	
	link <- set_link(link, par, data)
	
	link_data <- lapply(link, make.link)
	
	linklist <- lapply(c("linkfun", "linkinv", "mu.eta"),
					   function(x) lapply(link_data, "[[", x))
	
	names(linklist) <- c("linkfun", "linkinv", "mu.eta")
	
	newpar <- Map(function(x, y) ifelse(x=="identity", y, paste(x, y, sep="_")), x=link, y=par)
	newpar <- unname(unlist(newpar))
	
	names(linklist$linkfun) <- names(linklist$mu.eta) <- newpar
	linklist$newpar <- newpar
	
	attr(link, "linklist") <- linklist
	
	ptextlist <- as.list(paste0("    params[", 1:length(par), "] = ", par, "; \n"))
	
	parg <- c("t", "t0", par)
	
	farg <- paste(c(paste0("real ", parg), "real[] y0"), collapse=", ")
	
	scode_functions <- paste0(
		"  real[] gfun(real t, \n",
		"              real[] y, \n",
		"              real[] params, \n",
		"              real[] x_r, \n",
		"              int[] x_i) { \n",
		"    real dydt[", length(grad), "]; \n\n",
		do.call(paste0, gtextlist), "\n",
		"    return dydt;\n  }\n\n",
		"  real[] ifun(", farg, ") { \n",
		"    real y_ini[", length(initial), "]; \n",
		"    real params[", length(par) , "];\n",
		paste(ptextlist, collapse=""),
		do.call(paste0, itextlist), "\n",
		"    return y_ini;\n",
		"  }\n\n",
		"  real[] simfun(", farg, ", real[] x_r, int[] x_i){\n",
		"    real y_hat[1,", nstate, "];\n",
		"    real params[", length(par) , "];\n",
		paste(ptextlist, collapse=""),
			"    y_hat = integrate_ode_", solver, "(gfun, y0, t0, rep_array(t, 1), params, x_r, x_i); \n\n",
		"    return y_hat[1,];\n",
		"  }\n\n",
		"  real mufun(", farg, ", real[] y_hat){\n",
		"    return (", deparse(oo), ");\n",
		"  }"
	)
	
	mm <- list(
		grad=grad,
		observation=observation,
		initial=initial,
		par=par,
		state=state,
		scode_functions=scode_functions,
		family=family,
		link=link,
		data=data,
		tcol=tcol
	)
	
	class(mm) <- c("list", "bomsmodel")
	
	mm
}
