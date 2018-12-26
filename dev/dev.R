library(boms)
library(rstan)

harbin <- fitsir::harbin

SI_model <- make_model(
	grad = list(
        S ~ - beta * S * I/N,
        I ~ beta * S * I/N - gamma * I
    ),
    observation = Deaths ~ I,
    initial = list(
        S ~ N * (1 - i0),
        I ~ N * i0
    ),
    par= c("beta", "gamma", "N", "i0"),
	link=c(i0="logit"),
	data=harbin, tcol="week",
	family="negbinomial"
)

mm <- make_stancode_boms(
	model=SI_model,
	effect=list(
		beta~1,
		gamma~1,
		N~1,
		i0~1,
		shape~1
	),
	prior = list(
		beta ~ gamma(2, 1),
		gamma ~ gamma(10, 10),
		N ~ gamma(2, .001),
		i0 ~ beta(1, 100)
	)
)

dd <- make_standata_boms(
	model=SI_model,
	effect=list(
		beta~1,
		gamma~1,
		N~1,
		i0~1,
		shape~1
	),
	prior = list(
		beta ~ gamma(2, 1),
		gamma ~ gamma(10, 10),
		N ~ gamma(2, .001),
		i0 ~ beta(1, 100)
	)
)

sm <- stan_model(model_code=mm)

oo <- sampling(sm, data=dd,
			   chains=1,
			   iter=2000)

save("oo", file="dev.rda")

ee <- extract(oo, include=TRUE)

hist(exp(ee$b_beta)) 
hist(exp(ee$b_gamma))
hist(exp(ee$b_N))
hist(plogis(ee$b_i0))
hist(exp(ee$b_shape_Intercept), breaks=50)

expose_stan_functions(oo)

bom(
	model=SI_model,
	prior=list(
		beta ~ gamma(2, 1),
		gamma ~ gamma(10, 10),
		N ~ gamma(2, .001),
		i0 ~ beta(1, 100)
	)
)
