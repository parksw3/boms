##' the initializer for model.ode
##'
##' @param .Object object
##' @slot name name of the model
##' @slot model ode model
##' @slot observation observation model
##' @slot initial initial values
##' @slot par parameters
##' @examples
##' SI_model <- new("model.ode",
##'     name = "SI",
##'     model = list(
##'         S ~ - beta*S*I/N,
##'         I ~ beta*S*I/N - gamma*I
##'     ),
##'     observation = list(
##'         susceptible ~ normal(mean=S, sd=sigma1),
##'         infected ~ normal(mean=I, sd=sigma2)
##'     ),
##'     initial = list(
##'         S ~ N * (1 - i0),
##'         I ~ N * i0
##'     ),
##'     par= c("beta", "gamma", "N", "i0", "sigma1", "sigma2")
##' )
##' @docType methods
##' @exportMethod initialize
setMethod(
    "initialize",
    "model.ode",
    function(.Object, name,
             model,
             observation,
             initial,
             par) {
    	## TODO: add more warnings and errors
    	if (any(sapply(model, class) != "formula"))
    		stop("model must be a list of formulas")
    	
        if (any(sapply(initial, class) != "formula"))
            stop("initial must be a list of formulas")
    	
    	if (any(sapply(initial, class) != "formula"))
    		stop("observation must be a list of formulas")

        state <- sapply(model, function(x) deparse(x[[2]]))
        nstate <- length(state)

        ## substitute state in gradient and initial
        tstate <- paste0("y[", 1:length(state),"]")
        
        ## substitute states in observation
        tobs <- paste0("y_hat[,", 1:length(state),"]")
        
        ## substitute parameters
        tpar <- paste0("params[", 1:length(par), "]")
        
        state_transforms <- sfun(state, tstate)
        
        obs_transforms <- sfun(state, tobs)
        
        param_transforms <- sfun(par, tpar)
        
        gtextlist <- itextlist <- vector('list', length(model))
        otextlist <- vector('list', length(observation))
        
        for (i in 1:length(model)) {
        	tt <- subst(subst(model[[i]][[3]], state_transforms), param_transforms)
        	
        	gtextlist[[i]] <- paste0("dydt[", i, "] = ", do.call(paste0, as.list(deparse(tt))), "; \n")
        }
        
        for (i in 1:length(initial)) {
        	ii <- subst(initial[[i]][[3]], param_transforms)
        	
        	itextlist[[i]] <- paste0("y0[", i, "] = ", do.call(paste0, as.list(deparse(ii))), "; \n")
        }
        
        for (i in 1:length(observation)) {
        	oo <- subst(subst(observation[[i]][[3]], obs_transforms), param_transforms)
        	
        	otextlist[[i]] <- paste0(deparse(observation[[i]][[2]])," ~ ", do.call(paste0, as.list(deparse(oo))), "; \n")
        }
        
        stancode_functions <- paste0(
        	"// generated with fitode2", "\n",
        	"functions { \n",
        	"real[] gfun(real t, \n",
        	"real[] y, \n",
        	"real[] params, \n",
        	"real[] x_r, \n",
        	"int[] x_i) { \n\n",
        	"real dydt[", length(model), "]; \n\n",
        	do.call(paste0, gtextlist), "\n",
        	"return dydt;\n}\n\n",
        	"real[] ifun(real[] params) { \n\n",
        	"real y0[", length(initial), "]; \n\n",
        	do.call(paste0, itextlist), "\n",
        	"return y0;\n",
			"}\n\n} \n"
        )
        
		stancode_observation <- do.call(paste0, otextlist)
		
		.Object@name <- name
		.Object@model <- model
		.Object@observation <- observation
		.Object@initial <- initial
		.Object@par <- par
		.Object@state <- state
		.Object@stancode_functions <- stancode_functions
		.Object@stancode_observation <- stancode_observation

        .Object
    }
)

setMethod("show", "model.ode",
    function(object){
        cat("Name:", object@name, "\n")

    	h <- paste0("d", object@state, "/dt = ", sapply(object@model, function(x) deparse(x[[3]])))
    	cat("\nModel:\n")
    	for(i in 1:length(object@model)) 
    		cat(h[i], "\n")
    	
    	cat("\nObservations:\n")
        for(i in 1:length(object@observation)) {
            cat(deparse(object@observation[[i]]), "\n")
        }

        cat("\nInitial values:\n")
        g <- paste0(object@state, "(0) = ", sapply(object@initial, function(x) deparse(x[[3]])))
        for(i in 1:length(g))
            cat(g[i], "\n")

        cat("\nParameters:", object@par, "\n")
    }
)
