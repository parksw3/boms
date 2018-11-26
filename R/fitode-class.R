##' Class representing ode models
##'
##' @slot name name of the model
##' @slot gfun gradient function
##' @slot grad list of gradients
##' @slot observation list of observation models
##' @slot initial list of expressions representing the initial values
##' @slot state state variables
##' @slot par parameters
##' @slot diffnames state variables that needs to be diff-ed
##' @slot jacobian.initial Jacobian of initial values with respect to its parameters
##' @slot jacobian.state Jacobian with respect to its states
##' @slot jacobian.par Jacobian with repsect to its parameters
##' @slot loglik list of log-likelihood functions
##' @slot expr expressions for true trajectories
##' @slot expr.sensitivity sensitivity of the expressions with respect to state variables and parameters
##' @slot keep_sensitivity (logical) keep sensitivity equations?
##' @exportClass model.ode
make_model <- setClass(
    "model.ode",
    slots = c(
        name = "character",
        model = "list",
        observation = "list",
        initial= "list",
        par = "character",
        state = "character",
        stancode_functions = "character",
        stancode_observation = "character"
    )
)

setClass(
	"stancode.ode",
	slots = c(
		stancode="character",
		model="model.ode",
		fixed="character",
		linklist="list"
	)
)

##' Class "fitode".
##' Result of ode fitting based on Maximum Likelihood Estimation
##'
##' @name fitode
##' @rdname fitode
##' @seealso \code{\link{mle2-class}}
##' @exportClass fitode
##' @export
setClass("fitode",
    slots = c(
        model="model.ode",
        data="data.frame",
        coef="numeric",
        vcov="matrix",
        min="numeric",
        mle2="mle2",
        link="list",
        fixed="list"
    )
)
