library(boms)

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

oo <- bom(
	model=SI_model,
	prior=list(
		beta ~ gamma(2, 1),
		gamma ~ gamma(10, 10),
		N ~ gamma(250, .01),
		i0 ~ beta(1, 100)
	)
)

brms:::extract_draws.btl

brms:::extract_draws.brmsterms

brms:::prepare_family

brms:::predict.brmsfit

brms:::predict_internal.brmsdraws

