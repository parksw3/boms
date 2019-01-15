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
		gamma ~ gamma(50, 50),
		N ~ gamma(2, .001),
		i0 ~ beta(1, 100)
	),
	seed=101
)
