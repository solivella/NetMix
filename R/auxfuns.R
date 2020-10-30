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
#' @param block_list List of matrices; each element is a square, numeric matrix 
#'                   that defines a blockmodel,
#' @param target_mat Numeric matrix; reference blockmodel that those in block_list should 
#'                   be aligned to. Optional, defaults to \code{NULL}.
#' @param use_perms Boolean; should all row/column permutations be explored when
#'                  realigning matrices? defaults to \code{TRUE}.
#' @param X Numeric matrix; design matrix of monadic predictors.
#' @param beta Numeric array; array of coefficients associated with monadic predictors. 
#'             It of dimensions Nr. Predictors by Nr. of Blocks by Nr. of HMM states.
#' @param alpha_list List of mixed-membership parameter matrices. 
#' @param kappa Numeric matrix; matrix of marginal HMM state probabilities.
#' @param C_mat Numeric matrix; matrix of posterior counts of block instantiations per node. 
#' @param y Numeric vector; vector of edge values.
#' @param d_id Integer matrix; two-column matrix with nr. dyads rows, containing zero-based
#'             sender (first column) and receiver (second column) node id's for each dyad. 
#' @param pi_mat Numeric matrix; row-stochastic matrix of mixed-memberships. 
#' @param colPalette A function produced by \code{colorRamp}.
#' @param range The range of values to label the legend.
#' @param par Vector of parameter values.
#' @param tot_nodes Integer vector; total number of nodes each node interacts with.
#' @param c_t Integer matrix; samples from Poisson-Binomial counts of a node instantiating a group.
#' @param x_t Numeric matrix; transposed monadic design matrices.
#' @param s_mat Integer matrix; Samples of HMM states by time period. 
#' @param t_id Integer vector; for each node, what time-period is it observed in? zero-indexed.
#' @param mu_beta,var_beta Numeric arrays; prior mean and variances of monadic coefficients.
#' @param orig Object to be transformed.
#' @param is_var Boolean. Is the object to be transformed a variance term?
#' @param is_array Boolean. Is the object to be transformed an array?
#' @param des.mat Numeric matrix. Design matrix corresponding to transformed object.
#' @param nblock Number of groups in model, defaults to \code{NULL}.
#' @param nstate Number of hidden Markov states in model, defaults to \code{NULL}.
#' @param devs Vector of standard deviations to use in transformation of variances. Defaults to \code{NULL}.
#' @param soc_mats,Y,dyads,edges,t_id_d,t_id_n,nodes_pp,dyads_pp,nt_id,node_id_period,mu_b,var_b,n_dyads,n.blocks,periods,directed,ctrl Internal arguments for initialization.
#' @param ... Numeric vectors; vectors of potentially different length to be cbind-ed.
#' 
#' @author Santiago Olivella (olivella@@unc.edu), Adeline Lo (aylo@@wisc.edu), Tyler Pratt (tyler.pratt@@yale.edu), Kosuke Imai (imai@@harvard.edu)
#'
#' @return See individual return section for each function:
#' \describe{
#'       \item{.cbind.fill}{Matrix of \code{cbind}'ed elements in \code{...}, with missing values in each vector filled with \code{NA}.}
#'       \item{.mpower}{Matrix; the result of raising \code{mat} to the \code{p} power.}
#'       \item{.findPerm}{List of permuted blockmodel matrices}
#'       \item{.transf}{Matrix with transformed mixed-membership vectors along its rows, s.t. no element is equal to 0.0 or 1.0.}
#'       \item{.compute.alpha}{List of predicted alpha matrices, one element per HMM state.}
#'       \item{.e.pi}{Matrix of expected mixed-membership vectors along its rows, with expectation computed over marginal 
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
.transf_muvar <- function(orig, is_var, is_array, des.mat, nblock=NULL, nstate=NULL,devs=NULL){
  if(is_array){
    tmp <- array(ifelse(is_var,1.0, 0.0), c(ncol(des.mat), nblock, nstate))
    rownames(tmp) <- colnames(des.mat)
  } else {
    tmp <- array(ifelse(is_var,1.0, 0.0), ncol(des.mat))
    names(tmp) <- colnames(des.mat)
  }
  if(length(orig) == 0){
    if(length(dim(tmp))==1){
      non_miss <- !grepl("_missing", colnames(des.mat))
      tmp[non_miss] <- orig
    } else {
      tmp[non_miss,,] <- orig
    }
  } else {
    if(length(dim(tmp))==1){
      tmp[names(orig)] <- orig
    } else {
      tmp[rownames(orig),,] <- orig
    }
  }
  if(is_var){
    if(length(dim(tmp))==1){
      tmp <- sqrt(tmp)/devs
    } else {
      tmp <- sqrt(tmp)/rep(devs, nblock*nstate)
    }
    tmp <- tmp^2
  }
  return(tmp)
}

## Adapted from `fields`` package under GPL
#' @rdname auxfuns
.bar.legend <- function(colPalette, range){
  col <- rgb(colPalette(seq(0,1, length.out = 100)), maxColorValue = 255)
  zlim <- c(0, 1)
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar), add = TRUE)
  stick <- TRUE
  nlevel <- length(col)
  midpoints <- seq(zlim[1], zlim[2], length.out=nlevel)
  delta <- (midpoints[2] - midpoints[1])/2
  breaks <- c(midpoints[1] - delta, midpoints + delta)
  legend.mar <- 5.1
  char.size <- par()$cin[1]/par()$din[1]
  offset <- char.size * par()$mar[4]
  legend.width <- char.size * 1.2
  legend.mar <- legend.mar * char.size
  smallplot <- opar$plt
  smallplot[2] <- 1 - legend.mar
  smallplot[1] <- smallplot[2] - legend.width
  pr <- (smallplot[4] - smallplot[3]) * ((1 - 0.9)/2)
  smallplot[4] <- smallplot[4] - pr
  smallplot[3] <- smallplot[3] + pr
  bigplot <- opar$plt
  bigplot[2] <- min(bigplot[2], smallplot[1] - offset)
  dp <- smallplot[2] - smallplot[1]
  smallplot[1] <- min(bigplot[2] + offset, smallplot[1])
  smallplot[2] <- smallplot[1] + dp
  
  if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
    stop("plot region too small to add legend\n")
  }
  ix <- 1:2
  iy <- breaks
  nBreaks <- length(breaks)
  midpoints <- (breaks[1:(nBreaks - 1)] + breaks[2:nBreaks])/2
  iz <- matrix(midpoints, nrow = 1, ncol = length(midpoints))
  
  par(new = TRUE, pty = "m", plt = smallplot, err = -1)
  graphics::image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
        ylab = "", col = col, breaks = breaks)
  axis.args <- c(list(side =  4, at = seq(0,1, length.out = 10),
                      labels = round(seq(range[1], range[2], length.out=10),2), 
                      cex.axis=0.75,
                      mgp = c(3, 1, 0), las = 2))
  
  do.call(graphics::axis, axis.args)
  graphics::box()
  legend.args <- list(text = "Edge\nProbability", side = 3, line = 0.5, cex = 1, adj = 0)
  do.call(mtext, legend.args)
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
.findPerm <- function(block_list, target_mat=NULL, use_perms=TRUE){
  tar_nprov <- is.null(target_mat) 
  if(tar_nprov){
    target_mat <- block_list[[1]] 
  }
  n <- ncol(target_mat)
  if(use_perms){
    all_perms <- gtools::permutations(n, n)
  } 
  res <- lapply(seq_along(block_list),
                function(x, t_mat, tar_nprov){
                  if(use_perms){
                    norms <- apply(all_perms, 1,
                                   function(p){
                                     base::norm(block_list[[x]][p, p]-t_mat, type="f")
                                   })
                    P <- as.matrix(as(all_perms[which.min(norms),], "pMatrix"))
                  } else {
                    P <- as.matrix(igraph::match_vertices(plogis(block_list[[x]]),
                                                          plogis(t_mat),
                                                          m = 0,
                                                          start = diag(ncol(block_list[[x]])),
                                                          iteration = 10)$P)
                  }
                  if(tar_nprov){
                    t_mat <- t(P) %*% block_list[[x]] %*% P
                  }
                  return(P)
                }, t_mat = target_mat, tar_nprov = tar_nprov)
  return(res)
}

#' @rdname auxfuns
.transf <- function(mat){
  (mat * (nrow(mat) - 1) + 1/ncol(mat))/nrow(mat)
}

#' @rdname auxfuns
.compute.alpha <- function(X, beta){
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
.e.pi <- function(alpha_list, kappa, C_mat = NULL){
  if(is.null(C_mat)){
    pi_l <- alpha_list
  } else {
    pi_l <- lapply(alpha_list, function(x) x + t(C_mat))
  }
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

#' @rdname auxfuns
.initPi <- function(soc_mats,
                     Y,
                     dyads,
                     edges,
                     t_id_d,
                     t_id_n,
                     nodes_pp,
                     dyads_pp,
                     nt_id,
                     node_id_period,
                     mu_b,
                     var_b,
                     n_dyads, n.blocks, periods, directed, ctrl){
  max_ll <- -Inf
  best_pi <- matrix(as.double(0.0), ncol=n_dyads, nrow=n.blocks)
  rownames(best_pi) <- 1:n.blocks
  for(n in 1:ctrl$nstarts){
    temp_res <- vector("list", periods)
    for(i in 1:periods){
      if(!ctrl$init_gibbs) {
        mn <- ncol(soc_mats[[i]]) 
        if(directed){
          D_o <- 1/sqrt(.rowSums(soc_mats[[i]], mn, mn) + 1)
          D_i <- 1/sqrt(.colSums(soc_mats[[i]], mn, mn) + 1)
          C_o <- t(D_o * soc_mats[[i]])
          C_i <- t(D_i * t(soc_mats[[i]]))
          U <- t(C_o * D_i) %*% C_o +
            t(C_i * D_o) %*% C_i
        } else {
          D <- 1/sqrt(.rowSums(soc_mats[[i]], mn, mn) + 1)
          U <- t(D * soc_mats[[i]]) * D
        }
        if(ctrl$spectral) {
          n_elem <- n.blocks + 1
          res <- RSpectra::eigs_sym(U, n_elem)
          eta <- res$vectors[,1:n_elem] %*% diag(res$values[1:n_elem])
          target <- eta[,2:n_elem] / (eta[,1] + 1e-8)
          sig <- 1 - res$values[n_elem] / (res$values[n.blocks])
          sig <- ifelse(is.finite(sig), sig, 0)
          if(sig > 0.1){
            target <- target[,1:(n_elem - 2), drop = FALSE]
          }
        } else {
          target <- U
        }
        if(nrow(unique(target)) > n.blocks){
          clust_internal <- fitted(kmeans(x = target,
                                          centers = n.blocks,
                                          iter.max = 15,
                                          nstart = 1), "classes")
          
        } else {
          init_c <- sample(1:nrow(target), n.blocks, replace = FALSE)
          cents <- jitter(target[init_c, ])
          clust_internal <- fitted(suppressWarnings(kmeans(x = target,
                                                           centers = cents,
                                                           iter.max = 15,
                                                           algorithm = "Lloyd",
                                                           nstart = 1)), "classes")
        }
        
        pi_internal <- model.matrix(~ factor(clust_internal, 1:n.blocks) - 1)
        pi_internal <- .transf(pi_internal)
        rownames(pi_internal) <- rownames(soc_mats[[i]])
        colnames(pi_internal) <- 1:n.blocks
        MixedMembership <- t(pi_internal)
        int_dyad_id <- apply(dyads[[i]], 2, function(x)match(x, colnames(MixedMembership)) - 1)
        BlockModel <- approxB(edges[[i]], int_dyad_id, MixedMembership)
        temp_res[[i]] <- list(BlockModel = BlockModel,
                              MixedMembership = MixedMembership)                         
        
      } else {
        n_prior <- (dyads_pp[i] - nodes_pp[i]) * .05
        a <- plogis(ctrl$mu_b) * n_prior
        b <- n_prior - a
        lda_beta_prior <- lapply(list(b,a),
                                 function(prior){
                                   mat <- matrix(prior[2], n.blocks, n.blocks)
                                   diag(mat) <- prior[1]
                                   return(mat)
                                 })
        ret <- lda::mmsb.collapsed.gibbs.sampler(network = soc_mats[[i]],
                                                 K = n.blocks,
                                                 num.iterations = 100L,
                                                 burnin = 50L,
                                                 alpha = ctrl$alpha,
                                                 beta.prior = lda_beta_prior)
        MixedMembership <- prop.table(ret$document_expects, 2)
        colnames(MixedMembership) <- colnames(soc_mats[[i]])
        int_dyad_id <- apply(dyads[[i]], 2, function(x)match(x, colnames(MixedMembership)) - 1)
        BlockModel <- approxB(edges[[i]], int_dyad_id, MixedMembership)
        if(any(is.nan(BlockModel))){
          BlockModel[is.nan(BlockModel)] <- 0.0
        }
        temp_res[[i]] <- list(BlockModel = BlockModel,
                              MixedMembership = MixedMembership)
        
      }
    }
    block_models <- lapply(temp_res, function(x)x$BlockModel)
    target_ind <- which.max(sapply(soc_mats, ncol))
    perms_temp <- .findPerm(block_models, target_mat = block_models[[target_ind]], use_perms = ctrl$permute)
    pis_temp <- lapply(temp_res, function(x)x$MixedMembership) 
    pi.ord <- as.numeric(lapply(pis_temp, function(x)strsplit(colnames(x), "@")[[1]][2])) # to get correct temporal order
    pi_init_t <- do.call(cbind,mapply(function(phi_tmp,perm){perm %*% phi_tmp},
                                       pis_temp[order(pi.ord)], perms_temp, SIMPLIFY = FALSE)) 
    rownames(pi_init_t) <- 1:n.blocks
    Z_tmp <- matrix(0, ncol = n_dyads, nrow = 1)
    X_tmp <- matrix(1, ncol = ncol(pi_init_t), nrow = 1)
    ctrl$vi_iter <- 2
    ctrl$verbose <- FALSE
    ##sparsity <- mean(Y >= 0.5)
    fit_tmp <- mmsbm_fit(Z_tmp,
                         ##Z_tmp,
                         X_tmp,
                         Y,
                         ##Y,
                         t_id_d,
                         t_id_n,
                         nodes_pp,
                         ##nt_id,
                         nt_id,
                         node_id_period,
                         mu_b,
                         var_b,
                         array(0.0, c(1, n.blocks, ctrl$states)),
                         array(1.0, c(1, n.blocks, ctrl$states)),
                         0.0,
                         1.0,
                         pi_init_t,
                         ctrl$kappa_init_t,
                         qlogis(block_models[[target_ind]]),
                         array(ctrl$alpha, c(1,n.blocks,ctrl$states)),
                         0.0,
                         ##sparsity,
                         ctrl)
    if(fit_tmp$LowerBound >= max_ll){
      best_pi <- pi_init_t
    }
  }
  return(best_pi)
}


