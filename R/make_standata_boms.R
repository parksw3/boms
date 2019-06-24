##' @export
make_standata_boms <- function(model,
							   effect=NULL,
							   prior,
							   data,
							   sample_prior = c("no", "yes", "only"),
							   family,
							   ...) {
	tcol <- model$tcol
	
	if (missing(data)) {
		bdata <- model$data
	} else {
		bdata <- data
	}
	
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
	
	bdata$t <- bdata[[tcol]]
	dt <- diff(bdata$t)
	## temporary
	bdata$t0 <- min(bdata$t) - min(dt[dt > 0])
	
	oformula <- as.formula(paste0(as.character(model$observation[[2]]), "~tmpfun(t, t0, ", paste(model$par, collapse=", "),  ")"))
	
	bf_arg <- c(oformula, effect, nl=TRUE)
	
	formula <- do.call(bf, bf_arg)
	formula$family <- family
	
	bterms <- parse_bf(formula)

	bprior <- do.call(c, lapply(prior, function(x) {
		arg <- c(x[[3]], nlpar=as.character(x[[2]]), lb=0)
		do.call(brms::prior, arg)
	}))
	
	sample_prior <- match.arg(sample_prior)
	prior <- check_prior(bprior, formula = formula, data = bdata, 
						 sample_prior = sample_prior, warn = TRUE)
	
	bdata <- brms:::update_data(bdata, bterms = bterms)
	
	allvars <- all.vars(bterms$allvars)
	allvars <- allvars[!allvars %in% c(deparse(model$observation[[2]]), "t", "t0")]
	
	bdata <- arrange_t0(bdata, allvars)
	
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
	
	out <- c(out, brms:::data_gr_global(ranef, cov_ranef = NULL))
	meef <- brms:::tidy_meef(bterms, data)
	out <- c(out, brms:::data_Xme(meef, data = data))
	
	out$which_t0 <- bdata$which_t0
	
	out$prior_only <- as.integer(identical(sample_prior, "only"))
	
	structure(out, class = "standata")
}

arrange_t0_internal <- function(x) {
	# x <- x[order(x$t),]
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

arrange_t0 <- function(data, allvars=NULL) {
	if (length(data$t0) == 0) {
		dt <- diff(data$t)
		data$t0 <- min(data$t) - min(dt[dt > 0])
	}
	
	if (length(allvars) > 0) {
		data <- do.call("rbind", lapply(split(data, data[allvars]), function(x){
			x[order(x$t),]
		}))
		
		data$order <- 1:nrow(data)
		
		data <- do.call("rbind", lapply(split(data, data[allvars]), function(x){
			if (nrow(x) == 1) {
				x$which_t0 <- x$order
			} else if (nrow(x) > 1) {
				x <- arrange_t0_internal(x)
			}
			x
		}))
	} else {
		data <- data[order(data$t),]
		data$order <- 1:nrow(data)
		data <- arrange_t0_internal(data)
	}
	## probably unnecessary
	data <- data[order(data$order),]
	rownames(data) <- NULL
	
	data
}
