## these will be model arguments?

model = list(
	S ~ - beta*S*I/N,
	I ~ beta*S*I/N - gamma*I
)

observation = list(
	susceptible ~ dnorm(mean=S, sd=sigma1),
	infected ~ dnorm(mean=I, sd=sigma2)
)

initial = list(
	S ~ N * (1 - i0),
	I ~ N * i0
)

par = c("beta", "gamma", "N", "i0", "sigma1", "sigma2")

tvar = "t"

## code starts here

## actual names
state <- sapply(model, function(x) deparse(x[[2]]))

## substitute names
tstate <- paste0("y[", 1:length(state),"]")

## substitute 

tpar <- paste0("params[", 1:length(par), "]")

tfun <- function(a, b) {
	fixed <- c(as.symbol("~"), parse(text=a))
	as.formula(as.call(c(fixed, parse(text=b))))
}

state_transforms <- mapply(tfun, a=state, b=tstate, SIMPLIFY=FALSE)

state_transforms2 <- trans(state_transforms, state)

param_transforms <- mapply(tfun, a=par, b=tpar, SIMPLIFY=FALSE)

param_transforms2 <- trans(param_transforms, par)



gtextlist <- list()

for (i in 1:length(model)) {
	s1 <- subst(model[[i]][[3]], state_transforms2)
	s2 <- subst(s1, param_transforms2)
	
	gtextlist[[i]] <- paste0("dydt[", i, "] = ", deparse(s2), "; \n")
}

scode_functions <- paste0(
	"// generated with fitode2", "\n",
	"functions { \n",
	"real[] gfun(real t, \n",
	"real[] y, \n",
	"real[] params, \n",
	"real[] x_r, \n",
	"real[] x_i, { \n\n",
	"real dydt[", length(model), "]; \n\n",
	do.call(paste0, gtextlist), "\n\n",
	"return dydt;\n\n}\n\n}"
)

cat(scode_functions)



