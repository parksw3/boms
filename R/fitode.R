##' fit ode
##' @rdname fitode
##' @name fitode
##' @param model ode model
##' @param data data frame with time column and observation column
##' @param init named vector of initial parameter values
##' @param tcol time column
##' @param method fit method
##' @param link named vector or list of link functions for ode/log-likelihood parameters
##' @param debug print debugging output?
##' @param ... rstan arguments
##' @import rstan
##' @seealso \code{\link{mle2}}
##' @export fitode
fitode <- function(model, data,
                   init, tcol="times", t0,
				   method=c("optimizing", "sampling"),
                   link,
                   fixed,
                   ...) {
	method <- match.arg(method)
	
    modelpar <- model@par
    
    if (!missing(fixed)) {
    	if (is.list(fixed)) fixed <- unlist(fixed)
    	if (any(!(names(fixed) %in% modelpar)))
    		stop("`fixed must be a named vector whose names correspond to model parameters")
    }
    
    if ("t" %in% modelpar) {
        stop("`t` is reserved for time variable. Try a different parameterization?")
    }
    
    stancode <- make_stancode(model,
    					fixed=fixed,
    					link=link)

    linklist <- stancode@linklist
    
    if (missing(init)) {
    	stop("Provide initial parameters")
    } else {
    	if (any(!(names(init) %in% modelpar))) stop("Initial parameters do not match model parameters")
    	
    	init <- init[match(names(init), modelpar)]
    	
    	newinit <- as.list(apply_link(init, linklist, "linkfun"))
    }
    
    ts <- data[[tcol]]
    
    if (missing(t0)) t0 <- min(ts) - min(diff(ts))
    
    standata <- list(
    	n_obs=nrow(data),
    	n_params=length(modelpar),
    	n_state=length(model@state),
    	t0=t0,
    	ts=ts
    )
    
    standata <- c(standata, data[sapply(model@observation, function(x) as.character(x[[2]]))])
    
    if (!missing(fixed)) standata <- c(standata, fixed)
    
    message("Compiling the C++ model")
    
    mm <- stan_model(model_code=stancode@stancode)
    
    message("Fitting the Stan model")
    
    stanarg <- list(
    	object=mm,
    	data=standata,
    	init=newinit,
    	...
    )
    
    if (method=="optimizing") stanarg <- c(stanarg, as_vector=FALSE)
    
    fitfun <- switch(method, optimizing=optimizing, sampling=sampling)
    
    stanfit <- do.call(fitfun, stanarg)
    
    stanfit
}
