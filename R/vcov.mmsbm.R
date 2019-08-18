#' Extract Variance-Covariance Matrix for a Fitted \code{mmsbm} Object
#'
#' 
#'
#' @param object An object of class \code{mmsbm}, a result of a call to \code{mmsbm}
#' @param param Character string, which set of parameters should the vcov be extracted for? One
#'              of \code{"MonadCoef"}, \code{"DyadCoef"}, \code{"BlockModel"} or \code{"All"} (the default).  
#' @param ... Currently ignored
#' @return For \code{param="DyadCoef"} and \code{param="BlockModel"}, a numeric matrix. For \code{param="MonadCoef"} a list of numeric matrices (one per HMM state) with rows/columns
#'         ordered first by predictor and then by group. For For \code{param="All"}, named list of individual return components.   
#'
#' @method vcov mmsbm
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
#'                       n.blocks = 2)
#' 
#' vcov(lazega_mmsbm, "MonadCoef")
#' 

vcov.mmsbm <- function(object,
                       param = "All",
                       ...)
{
  switch(param,
         MonadCoef = object$vcov_monad,
         DyadCoef = object$vcov_dyad,
         BlockModel = object$vcov_blockmodel,
         All = list(MonadCoef = object$vcov_monad,
                    DyadCoef = object$vcov_dyad,
                    BlockModel = object$vcov_blockmodel))
}