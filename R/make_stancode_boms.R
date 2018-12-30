##' @import brms
##' @importFrom brms bf parse_bf prior
##' @export
make_stancode_boms <- function(model,
							   effect,
							   prior,
							   data,
							   sample_prior = c("no", "yes", "only"),
							   ...) {
	tcol <- model$tcol
	
	if (missing(data)) {
		bdata <- model$data
	} else {
		bdata <- data
	}
	
	bdata$t <- bdata[[tcol]]
	bdata <- bdata[order(bdata$t),]
	dt <- diff(bdata$t)
	## TODO: fix this
	bdata$t0 <- min(bdata$t) - min(dt[dt > 0])
	
	oformula <- as.formula(paste0(as.character(model$observation[[2]]), "~tmpfun(t, t0, ", paste(model$par, collapse=", "),  ")"))
	
	bf_arg <- c(oformula, effect, nl=TRUE)
	
	formula <- do.call(bf, bf_arg)
	formula$family <- model$family
	
	bterms <- parse_bf(formula)

	bprior <- do.call(c, lapply(prior, function(x) {
		arg <- c(x[[3]], nlpar=as.character(x[[2]]), lb=0)
		do.call(brms::prior, arg)
	}))
	
	## TODO check for other options in sample_prior
	sample_prior <- match.arg(sample_prior)
	prior <- check_prior(bprior, formula = formula, data = bdata, 
						 sample_prior = sample_prior, warn = TRUE)
	
	bdata <- brms:::update_data(bdata, bterms = bterms)

	ranef <- brms:::tidy_ranef(bterms, data = bdata)
	meef <- brms:::tidy_meef(bterms, data = bdata)
	
	scode_predictor <- brms:::stan_predictor(bterms, data = bdata, prior = prior, 
									  ranef = ranef, meef = meef, sparse = FALSE, stanvars = NULL)
	
	scode_ranef <- brms:::stan_re(ranef, prior = prior, cov_ranef = NULL)

	scode_llh <- brms:::stan_llh(bterms, data = bdata)
	scode_global_defs <- brms:::stan_global_defs(bterms, prior = prior, 
										  ranef = ranef, cov_ranef = NULL)
	scode_Xme <- brms:::stan_Xme(meef, prior = prior)
	
	par_prior <- scode_predictor$prior
	
	scode_functions <- paste0(
		"// generated with boms ", utils::packageVersion("boms"), "\n",
		"// generated with brms ", utils::packageVersion("brms"), "\n",
		"functions { \n", model$scode_functions, "\n} \n"
	)
	
	scode_data <- paste0(
		"data { \n", "  int<lower=1> N;  // total number of observations \n", 
		scode_predictor$data, scode_ranef$data, scode_Xme$data, 
		"  int which_t0[N]; \n",
		"  int prior_only;  // should the likelihood be ignored? \n", 
		"} \n"
	)
	
	scode_transformed_data <- paste0(
		"transformed data { \n", 
		scode_global_defs$tdataD, scode_predictor$tdataD,
		"    real x_r[0];\n",
		"    int x_i[0];\n",
		"} \n"
	)

	scode_parameters <- paste0(gsub("<lower=0>", "", scode_predictor$par), scode_ranef$par, 
							   scode_Xme$par)
	
	scode_parameters <- paste0(
		"parameters { \n", 
		"  // parameters are defined on unconstrained scales\n",
		scode_parameters, 
	#	scode_rngprior$par, 
		"} \n"
	)
	
	estpar <- model$par[is.na(match(model$par, colnames(bdata)))]
	link <- model$link
	
	etextlist <- vector('list', length(estpar))
	for (i in 1:length(estpar)) {
		j <- match(estpar[i], model$par)
		ee <- paste0("    nlp_", estpar[i], "[n] = ", stan_link(link[j]), "(nlp_", estpar[i] ,"[n])")
		
		etextlist[[i]] <- paste0(ee, "; \n")
	}
	
	yhat_frame <- gsub(
		"    // compute non-linear predictor \n",
		"",
		gsub("mu\\[n\\]", "y_hat[n,]", scode_predictor$modelC4)
	)
	
	scode_model_loop <- paste0(
		scode_predictor$modelC2, 
		scode_predictor$modelC3, 
		paste(etextlist, collapse=""), ## rep_array(0.1, 0) actually doesn't do anything inside ifun
		"    if (which_t0[n]==n) {\n",
		gsub(")", ", rep_array(0.1, 0))", gsub("tmpfun", "ifun", gsub("y_hat", "  y0", yhat_frame))),
		"    } else {\n",
		"      y0[n,] = y_hat[which_t0[n],];\n",
		"    }\n",
		gsub(")", ", y0[n,], x_r, x_i)", gsub("tmpfun", "simfun", yhat_frame)),
		gsub(")", ", y0[n,], y_hat[n,])", gsub("tmpfun", "mufun", gsub("y_hat\\[n,]", "mu[n]", yhat_frame)))
	)
	
	if (isTRUE(nzchar(scode_model_loop))) {
		scode_model_loop <- paste0("  for (n in 1:N) { \n", scode_model_loop, 
								   "  } \n")
	}
	
	prior0 <- scode_predictor$prior
	prior_deriv <- vector('list', length(estpar))
	
	for (i in 1:length(estpar)) {
		j <- match(estpar[i], model$par)
		pp <- paste0("b_", estpar[i])
		tpp <- paste0(stan_link(link[j]), "(", pp, ")")
		dpp <- paste0(stan_link_deriv(link[j]), "(", pp, ")")
		
		prior0 <- gsub(pp, tpp, prior0)
		
		if (link[j] != "identity") {
			prior_deriv[i] <- paste0("  target += ", dpp, ";\n")
			
			if (link[j] == "logit")
				prior_deriv[i] <- gsub(")", "|0, 1)", prior_deriv[i]) ## TODO
		}
	}
	
	scode_prior <- paste0(
		prior0, 
		paste(prior_deriv, collapse = ""),
		scode_ranef$prior, 
		scode_Xme$prior, brms:::stan_prior(class = "", prior = prior))
	
	scode_transformed_parameters <- paste0(
		"transformed parameters { \n", 
		scode_predictor$tparD, scode_ranef$tparD, scode_Xme$tparD, 
		scode_ranef$tparC1, 
		"  real y_hat[N,", length(model$state), "];\n",
		"  real y0[N,", length(model$state), "];\n",
		scode_predictor$modelD, 
		scode_predictor$modelC1, 
		scode_predictor$modelCgp1, 
		scode_predictor$modelC5, 
		scode_model_loop, 
		"} \n"
	)
	
	scode_model <- paste0(
		"model { \n", 
		"  // priors including all constants \n", 
		scode_prior, 
		"  // likelihood including all constants \n", 
		"  if (!prior_only) { \n", 
		scode_llh,
		"  } \n", 
		"} \n"
	)
	
	## TODO: work on initial condition
	
	scode_generated_quantities <- paste0("generated quantities { \n", 
										 scode_predictor$genD, scode_ranef$genD, scode_Xme$genD, 
										 scode_predictor$genC, scode_ranef$genC, 
										 scode_Xme$genC, "} \n")
	
	complete_model <- paste0(scode_functions, scode_data, scode_transformed_data, 
							 scode_parameters, scode_transformed_parameters, scode_model, 
							 scode_generated_quantities)
	
	class(complete_model) <- c("character", "brmsmodel")
	complete_model
}
