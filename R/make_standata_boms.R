arrange_t0 <- function(x) {
	x <- x[order(x$t),]
	tmp <- rep(NA, nrow(x))
	for (i in 1:nrow(x)) {
		if (x$t[i] == min(x$t)) {
			tmp[i] <- x$order[i]
		} else {
			y <- x[x$t[i] > x$t,]
			tmp[i] <- y$order[which.max(y$t)]
			x$t0[i] <- y$t[which.max(y$t)]
		}
	}
	x$which_t0 <- tmp
	x
}

##' @export
make_standata_boms <- function(model,
							   effect,
							   prior,
							   ...) {
	dots <- list(...)
	
	tcol <- model$tcol
	
	bdata <- as.data.frame(model$data)
	bdata$t <- bdata[[tcol]]
	dt <- diff(bdata$t)
	## temporary
	bdata$t0 <- min(bdata$t) - min(dt[dt > 0])
	
	linklist <- attr(model$link, "linklist")
	oformula <- as.formula(paste0(as.character(model$observation[[2]]), "~tmpfun(t, t0, ", paste(model$par, collapse=", "),  ")"))
	
	bf_arg <- c(oformula, effect, nl=TRUE)
	
	bformula <- do.call(bf, bf_arg)
	
	formula <- brms:::validate_formula(bformula, data = data, family = model$family)
	
	bterms <- parse_bf(formula)

	bprior <- do.call(c, lapply(prior, function(x) {
		arg <- c(x[[3]], nlpar=as.character(x[[2]]), lb=0)
		do.call(brms::prior, arg)
	}))
	
	sample_prior <- brms:::check_sample_prior("no")
	prior <- brms:::check_prior(bprior, formula = formula, data = bdata, 
								sample_prior = sample_prior, warn = TRUE)
	
	bdata <- brms:::update_data(bdata, bterms = bterms)
	bdata$order <- 1:nrow(bdata)
	
	allvars <- all.vars(bterms$allvars)
	allvars <- allvars[!allvars %in% c(deparse(model$observation[[2]]), "t", "t0")]
	
	if (length(allvars) > 0) {
		bdata <- do.call("rbind", lapply(split(bdata, bdata[allvars]), function(x){
			if (nrow(x) == 1) {
				x$which_t0 <- x$order
			} else if (nrow(x) > 1) {
				x <- arrange_t0(x)
			}
			x
		}))
	} else {
		bdata <- arrange_t0(bdata)
	}
	
	rownames(bdata) <- NULL
	bdata <- bdata[order(bdata$order),]
	
	out <- c(
		list(N = nrow(bdata)), 
		data_response(
			bterms, bdata
		)
	)
	
	ranef <- brms:::tidy_ranef(bterms, data = bdata)
	meef <- brms:::tidy_meef(bterms, data = bdata)
	
	out <- c(out, data_predictor(
		bterms, data = bdata, prior = prior, ranef = ranef, meef = meef
	))
	
	out$which_t0 <- bdata$which_t0
	
	out$prior_only <- as.integer(identical(sample_prior, "only"))
	
	structure(out, class = "standata")
}
