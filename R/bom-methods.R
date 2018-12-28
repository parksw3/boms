##' @export
print.bomsfit <- function(x, digits=2, ...) {
	print(summary(x, ...), digits=digits, ...)
}

##' @importMethodsFrom rstan summary
##' @export
summary.bomsfit <- function(object, priors = FALSE, prob = 0.95,
							mc_se = FALSE, use_cache = TRUE, ...) {
	model <- object$model
	effect <- object$effect
	
	out <- list(
		group=unique(object$ranef$group), 
		nobs=nrow(object$data),
		family=object$model$family,
		model=model,
		effect=effect
	)
	
	class(out) <- "bomssummary"
	
	oformula <- as.formula(paste0(as.character(model$observation[[2]]), "~tmpfun(t, t0, ", paste(model$par, collapse=", "),  ")"))
	
	bf_arg <- c(oformula, effect, nl=TRUE)
	
	formula <- do.call(bf, bf_arg)
	formula$family <- model$family
	
	bterms <- parse_bf(formula)
	
	pars <- dimnames(object$fit)$parameters
	
	meta_pars <- object$fit@sim$pars_oi
	meta_pars <- meta_pars[!grepl("^(r|s|zgp|Xme|prior|lp)_", meta_pars)]
	probs <- c((1 - prob) / 2, 1 - (1 - prob) / 2)
	
	fit_summary <- rstan::summary(
		object$fit, pars=meta_pars,
		probs=probs, use_cache=use_cache
	)
	fit_summary <- fit_summary$summary
	if (!mc_se) {
		fit_summary <- fit_summary[, -2, drop = FALSE] 
	}
	CIs <- paste0(c("l-", "u-"), prob * 100, "% CI")
	colnames(fit_summary) <- c(
		"Estimate", if (mc_se) "MC.Error", 
		"Est.Error", CIs, "Eff.Sample", "Rhat"
	)
	
	## TODO: change scale to response
	
	fe_pars <- pars[grepl("^b(()|(s)|(cs)|(sp)|(mo)|(me)|(mi)|(m))_", pars)]
	out$fixed <- fit_summary[fe_pars, , drop = FALSE]
	rownames(out$fixed) <- gsub("^b(()|(s)|(cs)|(sp)|(mo)|(me)|(mi)|(m))_", "", fe_pars)
	
	## TODO ??
	# summary of family specific parameters
	# spec_pars <- c(valid_dpars(object), "delta")
	# spec_pars <- paste0(spec_pars, collapse = "|")
	# spec_pars <- paste0("^(", spec_pars, ")($|_)")
	# spec_pars <- pars[grepl(spec_pars, pars)]
	# out$spec_pars <- fit_summary[spec_pars, , drop = FALSE]
	
	# summary of residual correlations
	rescor_pars <- pars[grepl("^rescor_", pars)]
	if (length(rescor_pars)) {
		out$rescor_pars <- fit_summary[rescor_pars, , drop = FALSE]
		rescor_pars <- sub("__", ",", sub("__", "(", rescor_pars))
		rownames(out$rescor_pars) <- paste0(rescor_pars, ")")
	}
	
	# summary of autocorrelation effects
	# cor_pars <- pars[grepl(regex_cor_pars(), pars)]
	# out$cor_pars <- fit_summary[cor_pars, , drop = FALSE]
	# rownames(out$cor_pars) <- cor_pars
	
	# summary of group-level effects
	for (g in out$group) {
		gregex <- gsub(".", "\\.", g, fixed = TRUE)
		sd_prefix <- paste0("^sd_", gregex, "__")
		sd_pars <- pars[grepl(sd_prefix, pars)]
		cor_prefix <- paste0("^cor_", gregex, "__")
		cor_pars <- pars[grepl(cor_prefix, pars)]
		df_prefix <- paste0("^df_", gregex, "$")
		df_pars <- pars[grepl(df_prefix, pars)]
		gpars <- c(df_pars, sd_pars, cor_pars)
		out$random[[g]] <- fit_summary[gpars, , drop = FALSE]
		if (isTRUE(nrow(out$random[[g]]) > 0L)) {
			sd_names <- sub(sd_prefix, "sd(", sd_pars)
			cor_names <- sub(cor_prefix, "cor(", cor_pars)
			cor_names <- sub("__", ",", cor_names)
			df_names <- sub(df_prefix, "df", df_pars)
			gnames <- c(df_names, paste0(c(sd_names, cor_names), ")"))
			rownames(out$random[[g]]) <- gnames
		}
	}
	
	out
}
