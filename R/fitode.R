##' fit ode
##' @rdname fitode
##' @name fitode
##' @param model ode model
##' @param data data frame with time column and observation column
##' @param start named vector of starting parameter values
##' @param tcol time column
##' @param method optimization method
##' @param optimizer optimizer
##' @param link named vector or list of link functions for ode/log-likelihood parameters
##' @param control see optim
##' @param solver.opts options for ode integration. See ode
##' @param skip.hessian skip hessian calculation
##' @param use.ginv use generalized inverse (\code{\link{ginv}}) to compute approximate vcov
##' @param debug print debugging output?
##' @param ... mle2 arguments
##' @import bbmle
##' @importFrom numDeriv jacobian hessian
##' @importFrom MASS ginv
##' @seealso \code{\link{mle2}}
##' @export fitode
fitode <- function(model, data,
                   start, tcol="times",
                   method="BFGS",
                   optimizer="optim",
                   link,
                   fixed=list(),
                   control=list(maxit=1e5),
                   solver.opts=list(method="rk4"),
                   solver=ode,
                   skip.hessian=FALSE,
                   use.ginv=TRUE,
                   debug=FALSE,
                   ...) {

    modelpar <- model@par

    if ("t" %in% modelpar) {
        stop("`t` is reserved for time variable. Try a different parameterization?")
    }
    
    stancode <- make_stancode(model,
    					fixed=names(fixed),
    					link=link)

    
}
