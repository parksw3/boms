predict_internal <- function(x, newdata) {
	
	model <- x$model
	fit <- x$fit
	stanmodel <- fit@stanmodel
	
	bdata <- newdata
	
	tcol <- model$tcol
	
	bdata$t <- bdata[[tcol]]
	bdata <- bdata[order(bdata$t),]
	dt <- diff(bdata$t)
	## temporary
	bdata$t0 <- min(bdata$t) - min(dt[dt > 0])
	
	standata <- make_standata_boms(
		model,
		effect=x$effect,
		prior=x$prior,
		data=bdata,
		sample_prior="no"
	)
	
	ee <- rstan::extract(fit)
	
	mu_mat <- matrix(0, nrow=standata$N, ncol=length(ee$lp__))
	
	for (i in 1:length(ee$lp__)) {
		par <- list()
		
		for (j in 1:length(fit@model_pars)) {
			pp <- fit@model_pars[j]
			par[[pp]] <- select_indice(ee[[pp]], i) 
			if (grepl("temp", pp)) par[[pp]] <- c(par[[pp]])
			
		}
		
		output <- capture.output(
			new_fit <- rstan::sampling(
				stanmodel, data=standata,
				chains=1, iter=1,
				pars="mu",
				init=list(par),
				algorithm="Fixed_param", show_messages=FALSE,
				refresh=0,
				cores=1
			)
		)
		
		new_ee <- rstan::extract(new_fit)
		
		mu_mat[,i] <- c(new_ee$mu)
	}
	
	t(mu_mat)
}

select_indice <- function(x, n) {
	switch(sum(dim(x) != 1),
		   "1"=as.array(x[n]),
		   "2"=x[n,],
		   "3"=x[n,,])
}
