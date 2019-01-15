##' Fit Bayesian Ordinary Differential Equations (ODE) Multilevel Models
##' @import rstan
##' @importFrom rstan stan_model sampling
##' @export
bom <- function(model, 
				effect=NULL,
				prior,
				data, family, 
				# autocor = NULL, cov_ranef = NULL, 
				sample_prior = c("no", "yes", "only"), 
				# sparse = FALSE, knots = NULL, stanvars = NULL,
				# stan_funs = NULL, fit = NA, save_ranef = TRUE, 
				# save_mevars = FALSE, save_all_pars = FALSE, 
				# inits = "random", 
				chains = 1, iter = 2000, 
				warmup = floor(iter / 2), thin = 1,
				cores = getOption("mc.cores", 1L), 
				control = NULL,
				# algorithm = c("sampling", "meanfield", "fullrank"),
				# future = getOption("future", FALSE), 
				seed = NA, 
				# save_model = NULL, stan_model_args = list(),
				# save_dso = TRUE, file = NULL, 
				...) {
	if (missing(data)) {
		bdata <- model$data
		data.name <- model$data.name
	} else {
		bdata <- data
		data.name <- substr(collapse(deparse(substitute(data))), 1, 50)
	}
	
	tcol <- model$tcol
	
	bdata$t <- bdata[[tcol]]
	dt <- diff(bdata$t)
	## temporary
	bdata$t0 <- min(bdata$t) - min(dt[dt > 0])
	
	if (missing(family)) family <- model$family
	
	if (is.null(effect)) {
		enames <- NULL
	} else {
		enames <- lapply(effect, function(x) as.character(x[[2]]))
	}
	
	estpar <- c(model$par[is.na(match(model$par, colnames(bdata)))], family$dpars[-1])
	
	if (sum(!estpar %in% enames) > 0) {
		effect0 <- sapply(paste0(estpar[!estpar %in% enames], " ~ 1"), as.formula)
		names(effect0) <- NULL
		
		effect <- c(effect, effect0)
		effect <- effect[match(estpar, sapply(effect, function(x) as.character(x[[2]])))]
	}
	
	stancode_boms <- make_stancode_boms(
		model=model,
		effect=effect,
		prior=prior,
		data=bdata
	)
	
	standata_boms <- make_standata_boms(
		model=model,
		effect=effect,
		prior=prior,
		data=bdata
	)

	message("Compiling the C++ model")
	
	stanmodel_boms <- rstan::stan_model(model_code=stancode_boms)
	
	ss <- sampling(stanmodel_boms,
			 data=standata_boms,
			 chains=chains, iter=iter,
			 warmup=warmup, thin=thin,
			 control=control, seed=seed,
			 cores=cores,
			 ...)
	
	oformula <- as.formula(paste0(as.character(model$observation[[2]]), "~tmpfun(t, t0, ", paste(model$par, collapse=", "),  ")"))
	
	bf_arg <- c(oformula, effect, nl=TRUE)
	
	formula <- do.call(bf, bf_arg)
	formula$family <- model$family
	
	bterms <- parse_bf(formula)
	
	bdata <- brms:::update_data(bdata, bterms = bterms)
	
	allvars <- all.vars(bterms$allvars)
	allvars <- allvars[!allvars %in% c(deparse(model$observation[[2]]), "t", "t0")]
	
	bdata <- arrange_t0(bdata, allvars)
	
	x <- list(
		model=model,
		effect=effect,
		code=stancode_boms,
		data=bdata,
		standata=standata_boms,
		fit=ss,
		prior=prior,
		allvars=allvars
	)
	
	x$ranef <- brms:::tidy_ranef(bterms, data = bdata)
	
	## TODO: do I really need this?
	x$exclude <- brms:::exclude_pars(
		bterms, data = x$data, ranef = x$ranef
	)
	
	x$model.name <- substr(collapse(deparse(substitute(model))), 1, 50)
	x$data.name <- data.name
	
	class(x) <- c("bomsfit")
	
	x
}
