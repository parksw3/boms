library(fitode2)

harbin <- fitsir::harbin

SI_model <- fitode2:::make_model(
	name = "SI",
    model = list(
        S ~ - beta * S * I/N,
        I ~ beta * S * I/N - gamma * I
    ),
    observation = list(
        Deaths ~ neg_binomial_2(I, phi)
    ),
    initial = list(
        S ~ N * (1 - i0),
        I ~ N * i0
    ),
    par= c("beta", "gamma", "N", "i0", "phi")
)

init <- c(
	gamma=1,
	beta=2,
	N=2000,
	i0=1e-3,
	phi=10
)

ff <- fitode(SI_model, 
			 harbin, 
			 tcol="week", 
			 init=init,
			 link=c(i0="logit"))



