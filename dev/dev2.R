library(boms)
library(data.table)

load("lossData0.rda")

model <- make_model(
	grad = list(
		EX ~ - Ber * t * EX,
		OS ~ Ber * t * RLR * EX - kp * OS,
		PD ~ kp * RRF * OS
	),
	observation = loss_train ~ OS * (1 - delta) + PD * delta,
	initial = list(
		EX ~ premium,
		OS ~ 0,
		PD ~ 0
	),
	par= c("Ber", "RLR", "RRF", "kp", "premium", "delta"),
	data=lossData0[cal <= max(accident_year) & dev > 0],
	tcol="dev",
	link=c(Ber="identity", RLR="identity", RRF="identity", kp="identity"),
	family="gaussian"
)

mm <- make_stancode_boms(
	model=model,
	effect=list(
		RLR ~ 1 + (1 | p | accident_year), # 'p' allow for correlation with RRF 
		RRF ~ 1 + (1 | p | accident_year), # 'p' allow for correlation with RLR
		sigma ~ 0 + deltaf
	),
	prior=list(
		RLR ~ gamma(4, 5),
		RRF ~ gamma(4, 5),
		Ber ~ gamma(12, 3),
		kp ~ gamma(3, 4)
	)
)

dd <- make_standata_boms(
	model=model,
	effect=list(
		RLR ~ 1 + (1 | p | accident_year), # 'p' allow for correlation with RRF 
		RRF ~ 1 + (1 | p | accident_year), # 'p' allow for correlation with RLR
		sigma ~ 0 + deltaf
	),
	prior=list(
		RLR ~ gamma(4, 5),
		RRF ~ gamma(4, 5),
		Ber ~ gamma(12, 3),
		kp ~ gamma(3, 4)
	)
)

rt <- rstan::stanc(model_code=mm)
sm <- rstan::stan_model(stanc_ret=rt)


## similar estimates as 
## https://magesblog.com/post/2018-01-30-pkpd-reserving-models/
oo <- rstan::optimizing(sm, data=dd, seed=101, as_vector=FALSE,
						iter=4000)
## don't run; takes too long
ss <- rstan::sampling(sm, data=dd, seed=101, iter=500, chains=1)

