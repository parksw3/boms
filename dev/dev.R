library(rstan)
library(tsiR)

london <- twentymeas$London[1:100,]

SI_model <- make_model(
	name = "SI",
    model = list(
        S ~ mu * (N-S) - (c0+c1*(1+cos(2*3.141592*t)))*S*I/N,
        E ~ (c0+c1*(1+cos(2*3.141592*t)))*S*I/N - (mu + sigma) * E,
        I ~ sigma * E - (mu + gamma)*I
    ),
    observation = list(
        cases ~ neg_binomial_2(I, phi)
    ),
    initial = list(
        S ~ N * s0,
        E ~ N * e0,
        I ~ N * i0
    ),
    par= c("c0", "c1", "mu", "sigma", "gamma", "phi", "N", "s0", "e0", "i0")
)

sm <- make_stancode(SI_model,
					fixed=c("mu", "sigma", "gamma", "N"),
					link=c(s0="logit", e0="logit", i0="logit"))

mm <- stan_model(model_code=sm)

standata <- list(
	n_obs=nrow(london),
	n_params=10,
	n_state=3,
	t0=min(london$time)-1/26,
	ts=london$time,
	mu=1/50,
	sigma=35.6,
	gamma=33,
	N=london$pop[1],
	cases=london$cases
)

oo <- optimizing(mm, data=standata, as_vector=FALSE,
				 init=list(
				 	log_c0=log(300), log_c1=log(100), log_phi=5,
				 	logit_s0=qlogis(0.04), logit_e0=qlogis(1e-7), logit_i0=qlogis(1e-6))
				 )

oo$value

plot(london$cases)
lines(oo$par$y_hat[,3], type="l")


ss <- sampling(mm, data=standata, 
			   chains=1,
			   init=list(list(
			   	log_c0=log(20), log_c1=log(3), log_phi=1,
			   	logit_s0=qlogis(0.04), logit_e0=qlogis(1e-10), logit_i0=qlogis(1e-7))))
 
traceplot(ss, pars=c("beta", "gamma", "i0"))

ext <- extract(ss)

plot(standata$infected)

lines(apply(ext$y_hat[,,2], 2, mean))
matlines(t(apply(ext$y_hat[,,2], 2, quantile, c(0.025, 0.975))), col=1, lty=2)


