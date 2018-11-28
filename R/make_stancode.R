##' Set up link functions for ode/loglik parameters
##'
##' @param link named list or vector of strings specifying link functions
##' @param model model.ode object
##' @param loglik loglik.ode object
##' @seealso \code{\link{make.link}}
##' @return list of strings specifying link functions
set_link <- function(link, modelpar) {
	link_default <- as.list(rep("log", length(modelpar)))
	names(link_default) <- modelpar
	
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
		probit="Phi",
		cloglog="inv_cloglog",
		cauchit="atan",
		log="exp",
		loglog="inv_logit"
	)
}

input_type <- function(dd) {
	ifelse(
		dd %in% c("neg_binomial_2", "binomial", "beta_binomial", "poisson"),
		"int",
		"real"
	)
}

## took the name from brms

##' write a stan code
make_stancode <- function(model, 
						  link,
						  fixed,
						  solver=c("rk45", "bdf")) {
	solver <- match.arg(solver)
	
	par <- model@par

	observation <- model@observation
	
	if (!missing(link)) {
		link <- link[names(link) %in% par]
		
		if (any(is.na(match(names(link), par)))) stop("Some link functions do not correspond to the model parameters.")
	}
	
	if (!missing(fixed)) {
		which_fix <- match(fixed, par)
		which_est <- (1:length(par))[-which_fix]
	} else {
		which_est <- 1:length(par)
	}
	
	if (!missing(link)) {
		link <- link[names(link) %in% par]
		
		if (any(is.na(match(names(link), par)))) 
			stop("Some link functions do not correspond to the model parameters.")
	}
	
	link <- set_link(link, par)[which_est]
	
	link_data <- lapply(link, make.link)
	
	linklist <- lapply(c("linkfun", "linkinv", "mu.eta"),
					   function(x) lapply(link_data, "[[", x))
	
	names(linklist) <- c("linkfun", "linkinv", "mu.eta")
	
	newpar <- Map(function(x, y) ifelse(x=="identity", y, paste(x, y, sep="_")), x=link, y=par[which_est])
	newpar <- unname(unlist(newpar))
	
	names(linklist$linkfun) <- names(linklist$mu.eta) <- newpar
	linklist$newpar <- newpar
	
	etextlist <- vector('list', length(which_est))
	for (i in 1:length(which_est)) {
		ee <- paste0("real ", par[which_est][i], " = ", stan_link(link[[i]]), "(", newpar[i], ")")
		
		if (link[[i]]=="cauchit") ee <- paste0(ee, "/pi() + 0.5")
		
		etextlist[[i]] <- paste0(ee, "; \n")
	}
	
	if (!missing(fixed)) {
		ftextlist <- as.list(paste0("real ", par[which_fix], "; \n"))
	} else {
		ftextlist <- list()
		fixed <- character(0)
	}
	
	ptextlist <- as.list(paste0("params[", 1:length(par), "] = ", par, "; \n"))
	
	sname <- sapply(observation, function(x) deparse(x[[2]]))
	dtype <- sapply(observation, function(x) deparse(x[[3]][[1]]))
	
	dtextlist <- as.list(paste0(input_type(dtype), " ", sname, "[n_obs]; \n"))
	
	## TODO: do observation thing ...
	stancode_data <- paste0(
		"data {\n",
		"int<lower = 1> n_obs; // Number of days sampled \n",
		"int<lower = 1> n_params; // Number of model parameters \n",
		"int<lower = 1> n_state; \n",		
		"real t0; \n",
		"real ts[n_obs]; \n\n",
		do.call(paste0, ftextlist), "\n",
		do.call(paste0, dtextlist),
		"} \n"
	)
	
	stancode_tdata <- paste0(
		"transformed data { \n",
		"real x_r[0]; \n",
		"int x_i[0]; \n",
		"} \n"
	)
	
	stancode_par <- paste0(
		"parameters { \n",
		do.call(paste0, as.list(paste0("real ", newpar, "; \n"))),
		"} \n"
	)
	
	stancode_tpar <- paste0(
		"transformed parameters{ \n",
		"real y_hat[n_obs, n_state]; \n",
		"real y0[n_state]; \n",
		"real params[n_params]; \n\n",
		do.call(paste0, etextlist), "\n",
		do.call(paste0, ptextlist), "\n",
		"y0 = ifun(params); \n\n",
		"y_hat = integrate_ode_", solver, "(gfun, y0, t0, ts, params, x_r, x_i); \n\n",
		"} \n"
	)
	
	## TODO: prior
	stancode_model <- paste0(
		"model{ \n",
		model@stancode_observation,
		"} \n"
	)
	
	complete_model <- paste0(
		model@stancode_functions,
		stancode_data, 
		stancode_tdata, 
		stancode_par,
		stancode_tpar,
		stancode_model
	)
	
	new("stancode.ode",
		stancode=complete_model,
		model=model,
		fixed=fixed,
		linklist=linklist
	)
}
