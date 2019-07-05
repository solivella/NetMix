#' Dyadic predictors in the Lazega friendship network (Lazega 2001).
#'
#' A dataset containing edges and dyad-level predictors in the network
#' of friendships among lawyers in a New England law firm. More 
#' details are available in Lazega (2001).
#'
#' @format A data frame with 5041 rows and 4 variables:
#' \describe{
#'   \item{Lawyer1, Lawyer2}{lawyer ID, corresponding to identifiers common
#'          to those in \code{lazega_monadic}; numeric}
#'   \item{SocializeWith}{value of edge in network; binary}
#'   \item{Coworkers}{are the corresponding lawyers in the same office? boolean}
#' }
#' @usage data(lazega_dyadic)
#' @source \url{https://github.com/Z-co/networkdata/blob/master/networkdata/data/lazega.rda}
#' @references Emmanuel Lazega, \emph{The Collegial Phenomenon: The Social Mechanisms of Cooperation Among Peers in a Corporate Law Partnership}, Oxford University Press (2001).
"lazega_dyadic"

#' Monadic predictors in the Lazega friendship network (Lazega 2001).
#'
#' A dataset containing vertex-level predictors in the network
#' of sought-after advise among lawyers in a New England law firm. More 
#' details are available in Lazega (2001).
#'
#' @format A data frame with 71 rows and 7 variables:
#' \describe{
#'   \item{Lawyer}{lawyer ID,corresponding to identifiers common
#'          to those in \code{lazega_dyadic}; numeric}
#'   \item{Age}{age, in years; numeric}
#'   \item{Gender}{1=man; 2=woman; factor}
#'   \item{School}{1=harvard, yale; 2=ucon; 3= other; factor}
#'   \item{Practice}{1=litigation; 2=corporate; factor}
#'   \item{Seniority}{time in the firm, in years; numeric}
#'   \item{Status}{1=partner; 2=associate; factor}
#' }
#' @usage data(lazega_monadic)
#' @source Emmanuel Lazega, \emph{The Collegial Phenomenon: The Social Mechanisms of Cooperation Among Peers in a Corporate Law Partnership}, Oxford University Press (2001).
#' @source \url{https://github.com/Z-co/networkdata/blob/master/networkdata/data/lazega.rda}
"lazega_monadic"

