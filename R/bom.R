##' Fit Bayesian Ordinary Differential Equations (ODE) Multilevel Models
##' @import rstan
##' @importFrom rstan stan_model
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
	tcol <- model$tcol
	
	if (missing(data)) {
		bdata <- model$data
	} else {
		bdata <- data
	}
	
	bdata$t <- bdata[[tcol]]
	bdata <- bdata[order(bdata$t),]
	
	if (missing(family))
		family <- model$family
	
	if (is.null(effect)) {
		enames <- NULL
	} else {
		enames <- lapply(effect, function(x) as.character(x[[2]]))
	}
	
	estpar <- c(model$par[is.na(match(model$par, colnames(bdata)))], family$dpars[2])
	
	effect0 <- sapply(paste0(estpar[!estpar %in% enames], " ~ 1"), as.formula)
	names(effect0) <- NULL
	
	effect <- c(effect, effect0)
	effect <- effect[match(estpar, sapply(effect, function(x) as.character(x[[2]])))]
	
	stancode_boms <- make_stancode_boms(
		model=model,
		effect=effect,
		prior=prior
	)
	
	standata_boms <- make_standata_boms(
		model=model,
		effect=effect,
		prior=prior
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
	
	x <- list(
		model=model,
		effect=effect,
		code=stancode_boms,
		data=bdata,
		standata=standata_boms,
		fit=ss
	)
	
	class(x) <- c("bom")
	
	x
}
