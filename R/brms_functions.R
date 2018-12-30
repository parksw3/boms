## using unexported functions from brms 2.7.0
## copied from 
## https://github.com/paul-buerkner/brms/tree/56dd64f08295ced87167fad79a9ef502199775c3
## slightly modified to skip some steps

rename <- function(x, pattern = NULL, replacement = NULL, 
				   fixed = TRUE, check_dup = FALSE, ...) {
	# rename certain patterns in a character vector
	# Args:
	#   x: a character vector to be renamed
	#   pattern: the regular expressions in x to be replaced
	#   replacement: the replacements
	#   fixed: same as for sub, grepl etc
	#   check_dup: logical; check for duplications in x after renaming
	#   ...: passed to gsub
	# Returns: 
	#   renamed character vector of the same length as x
	pattern <- as.character(pattern)
	replacement <- as.character(replacement)
	if (!length(pattern) && !length(replacement)) {
		# default renaming to avoid special characters in coeffcient names 
		pattern <- c(
			" ", "(", ")", "[", "]", ",", "\"", "'", 
			"?", "+", "-", "*", "/", "^", "="
		)
		replacement <- c(rep("", 9), "P", "M", "MU", "D", "E", "EQ")
	}
	if (length(replacement) == 1L) {
		replacement <- rep(replacement, length(pattern))
	}
	stopifnot(length(pattern) == length(replacement))
	# avoid zero-length pattern error
	has_chars <- nzchar(pattern)
	pattern <- pattern[has_chars]
	replacement <- replacement[has_chars]
	out <- x
	for (i in seq_along(pattern)) {
		out <- gsub(pattern[i], replacement[i], out, fixed = fixed, ...)
	}
	dup <- duplicated(out)
	if (check_dup && any(dup)) {
		dup <- x[out %in% out[dup]]
		stop2("Internal renaming led to duplicated names. \n",
			  "Occured for: ", collapse_comma(dup))
	}
	out
}

seq_rows <- function(x) {
	seq_len(NROW(x))
}

##' @importFrom brms get_prior
check_prior <- function(prior, formula, data, 
						sample_prior = c("no", "yes", "only"),
						warn = FALSE) {
	# check prior input and amend it if needed
	# Args:
	#   same as the respective parameters in brm
	#   warn: passed to check_prior_content
	# Returns:
	#   a data.frame of prior specifications to be used in stan_prior
	if (isTRUE(attr(prior, "checked"))) {
		# prior has already been checked; no need to do it twice
		# attributes may still need to be updated
		attr(prior, "sample_prior") <- sample_prior
		return(prior)
	}
	bterms <- parse_bf(formula)
	all_priors <- get_prior(formula, data, internal = TRUE)
	if (is.null(prior)) {
		prior <- all_priors  
	} else if (!is.brmsprior(prior)) {
		stop2("Argument 'prior' must be a 'brmsprior' object.")
	}
	# temporarily exclude priors that should not be checked
	no_checks <- !nzchar(prior$class)
	prior_no_checks <- prior[no_checks, ]
	prior <- prior[!no_checks, ]
	# check for duplicated priors
	prior$class <- rename(
		prior$class, c("^cor$", "^rescor$", "^corme$"), 
		c("L", "Lrescor", "Lme"), fixed = FALSE
	)
	rcols <- c("class", "coef", "group", "resp", "dpar", "nlpar")
	duplicated_input <- duplicated(prior[, rcols])
	if (any(duplicated_input)) {
		stop2("Duplicated prior specifications are not allowed.")
	}
	## check if parameters in prior are valid
	#if (nrow(prior)) {
	#	valid <- which(duplicated(rbind(all_priors[, rcols], prior[, rcols])))
	#	invalid <- which(!seq_rows(prior) %in% (valid - nrow(all_priors)))
	#	if (length(invalid)) {
	#		msg_priors <- .print_prior(prior[invalid, ])
	#		stop2(
	#			"The following priors do not correspond ", 
	#			"to any model parameter: \n",
	#			collapse(.print_prior(prior[invalid, ]), "\n"),
	#			"Function 'get_prior' might be helpful to you."
	#		)
	#	}
	#}
	prior$prior <- sub("^(lkj|lkj_corr)\\(", "lkj_corr_cholesky(", prior$prior)
	# check_prior_content(prior, warn = warn)
	# merge user-specified priors with default priors
	prior$new <- rep(TRUE, nrow(prior))
	all_priors$new <- rep(FALSE, nrow(all_priors))
	prior <- prior + all_priors
	prior <- prior[!duplicated(prior[, rcols]), ]
	# prior <- check_prior_special(prior, bterms = bterms, data = data)
	prior <- prior[with(prior, order(class, group, resp, dpar, nlpar, coef)), ]
	prior <- prior + prior_no_checks
	rownames(prior) <- NULL
	attr(prior, "sample_prior") <- sample_prior
	attr(prior, "checked") <- TRUE
	prior
}

update_data <- function(data, bterms, na.action = na.omit2,
						drop.unused.levels = TRUE,
						terms_attr = NULL, knots = NULL) {
	# Update data for use in brms functions
	# Args:
	#   data: the original data.frame
	#   bterms: object of class brmsterms
	#   na.action: function defining how to treat NAs
	#   drop.unused.levels: indicates if unused factor levels
	#     should be removed
	#   terms_attr: a list of attributes of the terms object of 
	#     the original model.frame; only used with newdata;
	#     this ensures that (1) calls to 'poly' work correctly
	#     and (2) that the number of variables matches the number 
	#     of variable names; fixes issue #73
	#   knots: a list of knot values for GAMMs
	# Returns:
	#   model.frame for use in brms functions
	if (missing(data)) {
		stop2("Argument 'data' is missing.")
	}
	if (isTRUE(attr(data, "brmsframe"))) {
		return(data)
	}
	if (is.null(knots)) {
		knots <- attr(data, "knots", TRUE)
	}
	data <- try(as.data.frame(data), silent = TRUE)
	if (is(data, "try-error")) {
		stop2("Argument 'data' must be coercible to a data.frame.")
	}
	if (!isTRUE(nrow(data) > 0L)) {
		stop2("Argument 'data' does not contain observations.")
	}
	bterms$allvars <- terms(bterms$allvars)
	attributes(bterms$allvars)[names(terms_attr)] <- terms_attr
	data <- data_rsv_intercept(data, bterms = bterms)
	missing_vars <- setdiff(all.vars(bterms$allvars), names(data))
	if (length(missing_vars)) {
		stop2("The following variables are missing in 'data':\n",
			  collapse_comma(missing_vars))
	}
	for (v in intersect(vars_keep_na(bterms), names(data))) {
		attr(data[[v]], "keep_na") <- TRUE
	}
	data <- model.frame(
		bterms$allvars, data, na.action = na.action,
		drop.unused.levels = drop.unused.levels
	)
	if (any(grepl("__|_$", colnames(data)))) {
		stop2("Variable names may not contain double underscores ",
			  "or underscores at the end.")
	}
	groups <- get_group_vars(bterms)
	data <- combine_groups(data, groups)
	data <- fix_factor_contrasts(data, ignore = groups)
	attr(data, "knots") <- knots
	attr(data, "brmsframe") <- TRUE
	data
}

collapse <- function(..., sep = "") {
	# wrapper for paste with collapse = ""
	paste(..., sep = sep, collapse = "")
}

print_format <- function(x, digits = 2, no_digits = "Eff.Sample") {
	# helper function to print summary matrices
	# in nice format, also showing -0.00 (#263)
	# Args:
	#   x: object to be printed; coerced to matrix
	#   digits: number of digits to show
	#   no_digits: names of columns for which no digits should be shown
	x <- as.matrix(x)
	digits <- as.numeric(digits)
	if (length(digits) != 1L) {
		stop2("'digits' should be a single numeric value.")
	}
	out <- x
	fmt <- paste0("%.", digits, "f")
	for (i in seq_cols(x)) {
		if (isTRUE(colnames(x)[i] %in% no_digits)) {
			out[, i] <- sprintf("%.0f", x[, i])
		} else {
			out[, i] <- sprintf(fmt, x[, i])
		}
	}
	print(out, quote = FALSE, right = TRUE)
	invisible(x)
}

seq_cols <- function(x) {
	seq_len(NCOL(x))
}
