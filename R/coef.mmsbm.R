#' Extract Regression Coefficients for a Fitted \code{mmsbm} Object
#'
#' 
#'
#' @param object An object of class \code{mmsbm}, a result of a call to \code{mmsbm}
#' @param param Character string, which set of parameters should the vcov be extracted for? One
#'              of \code{"MonadCoef"}, \code{"DyadCoef"} or \code{"All"} (the default).  
#' @param ... Currently ignored
#' @return For \code{param="DyadCoef"}, a numeric vector. For \code{param="MonadCoef"}, an array
#'         with HMM states along the third dimension. For \code{param="All"}, named list of individual return components.   
#'
#' @method coef mmsbm
#' 
#' @author Santiago Olivella (olivella@@unc.edu), Adeline Lo (adelinel@@princeton.edu), Tyler Pratt (tyler.pratt@@yale.edu), Kosuke Imai (imai@@harvard.edu)
#' 
#' 
#' @examples 
#' library(NetMix)
#' ## Load datasets
#' data("lazega_dyadic")
#' data("lazega_monadic")
#' ## Estimate model with 2 groups
#' set.seed(123)
#' lazega_mmsbm <- mmsbm(SocializeWith ~ Coworkers,
#'                       ~  School + Practice + Status,
#'                       senderID = "Lawyer1",
#'                       receiverID = "Lawyer2",
#'                       nodeID = "Lawyer",
#'                       data.dyad = lazega_dyadic,
#'                       data.monad = lazega_monadic,
#'                       n.blocks = 2,
#'                       mmsbm.control = list(hessian = FALSE))
#' 
#' coef(lazega_mmsbm, "MonadCoef")
#' 

coef.mmsbm <- function(object,
                       param = "All",
                       ...)
{
  switch(param,
         MonadCoef = object$MonadCoef,
         DyadCoef = object$DyadCoef,
         All = list(MonadCoef = object$MonadCoef,
                    DyadCoef = object$DyadCoef))
}