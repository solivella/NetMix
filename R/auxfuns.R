#' @name auxfuns
#' @title Internal functions and generics for \code{mmsbm} package
#' 
#' @description  These are various utilities and generic methods used by 
#' the main package function. 
#' 
#' @details These functions are meant for internal use only.
#' 
#' @param mat Numeric matrix
#' @param p Numeric scalar; power to raise matrix to.
#' @param block.list List of matrices; each element is a square, numeric matrix 
#'                   that defines a blockmodel,
#' @param target.mat Numeric matrix; reference blockmodel that those in block.list should 
#'                   be aligned to. Optional, defaults to \code{NULL}.
#' @param use.perms Boolean; should all row/column permutations be explored when
#'                  realigning matrices? defaults to \code{TRUE}.
#' @param X Numeric matrix; design matrix of monadic predictors.
#' @param beta Numeric array; array of coefficients associated with monadic predictors. 
#'             It of dimensions Nr. Predictors by Nr. of Blocks by Nr. of HMM states.
#' @param pi_l List of mixed-membership matrices. 
#' @param kappa Numeric matrix; matrix of marginal HMM state probabilities.
#' @param y Numeric vector; vector of edge values.
#' @param d_id Integer matrix; two-column matrix with nr. dyads rows, containing zero-based
#'             sender (first column) and receiver (second column) node id's for each dyad. 
#' @param pi_mat Numeric matrix; row-stochastic matrix of mixed-memberships. 
#' @param ... Numeric vectors; vectors of potentially different length to be cbind-ed.
#' 
#' @author Kosuke Imai (imai@@harvard.edu), Tyler Pratt (tyler.pratt@@yale.edu), Santiago Olivella (olivella@@unc.edu)
#'
#' @return See individual return section for each function:
#' \describe{
#'       \item{.cbind.fill}{Matrix of \code{cbind}'ed elements in \code{...}, with missing values in each vector filled with \code{NA}.}
#'       \item{.mpower}{Matrix; the result of raising \code{mat} to the \code{p} power.}
#'       \item{.findPerm}{List of permuted blockmodel matrices}
#'       \item{.transf}{Matrix with transformed mixed-membership vectors along its rows, s.t. no element is equal to 0.0 or 1.0.}
#'       \item{.pi.hat}{List of predicted mixed-membership matrices, one element per HMM state.}
#'       \item{.e.pi }{Matrix of expected mixed-membership vectors along its rows, with expectation computed over marginal 
#'                     distribution over HMM states for each time period.}
#'     }
#' 
#' @rdname auxfuns
.cbind.fill<-function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  rnames <- unique(unlist(lapply(nm, rownames)))
  mat <- matrix(NA, nrow = length(rnames), ncol = length(nm), 
                dimnames = list(rnames, 1:length(nm)))
  for (x in seq_along(nm)) {
    mat[rownames(nm[[x]]), x] <- nm[[x]]
  }
  return(mat)
}

#' @rdname auxfuns
.mpower <- function(mat, p){
  orig <- mat
  while(p > 1){
    mat <- mat %*% orig
    p <- p - 1 
  }
  return(mat)
}

#' @rdname auxfuns
.findPerm <- function(block.list, target.mat=NULL, use.perms=TRUE){
  tar_nprov <- is.null(target.mat) 
  if(tar_nprov){
    target.mat <- block.list[[1]] 
  }
  n <- ncol(target.mat)
  if(use.perms){
    all_perms <- gtools::permutations(n, n)
  } 
  res <- lapply(seq_along(block.list),
                function(x, t_mat, tar_nprov){
                  if(use.perms){
                    norms <- apply(all_perms, 1,
                                   function(p){
                                     base::norm(block.list[[x]][p, p]-t_mat, type="f")
                                   })
                    P <- as.matrix(as(all_perms[which.min(norms),], "pMatrix"))
                  } else {
                    P <- as.matrix(igraph::match_vertices(plogis(block.list[[x]]),
                                                          plogis(t_mat),
                                                          m = 0,
                                                          start = diag(ncol(block.list[[x]])),
                                                          iteration = 10)$P)
                  }
                  if(tar_nprov){
                    t_mat <- t(P) %*% block.list[[x]] %*% P
                  }
                  return(P)
                }, t_mat = target.mat, tar_nprov = tar_nprov)
  return(res)
}

#' @rdname auxfuns
.transf <- function(mat){
  (mat * (nrow(mat) - 1) + 1/ncol(mat))/nrow(mat)
}

#' @rdname auxfuns
.pi.hat <- function(X, beta){
  if(dim(beta)[3]==1){
    mu <- exp(X %*% beta[,,1])
    pi.states <- list(t(mu))
  } else {
    pi.states <- lapply(1:dim(beta)[3], function(m){
      mu <- exp(X %*% beta[,,m])
      return(t(mu))
    })
  }
  return(pi.states)
}

#' @rdname auxfuns
.e.pi <- function(pi_l, kappa){
  if(is.null(dim(kappa))){
    return(pi_l[[1]])
  } else {
    n_states <- nrow(kappa)
    pi.states <- lapply(1:n_states,
                        function(m){
                          pi_l[[m]] * rep(kappa[m,], each=nrow(pi_l[[m]])) 
                        })
    return(Reduce("+", pi.states))
  }
}

