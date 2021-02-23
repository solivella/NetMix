#' @name auxfuns
#' @title Internal functions and generics for \code{mmsbm} package
#' 
#' @description  These are various utilities and generic methods used by 
#' the main package function. 
#' 
#' @details These functions are meant for internal use only.
#' 
#' @param fomula,formula.dyad,formula.monad1,formula.monad2 A formula object.
#' @param data,data.dyad,data.monad1,data.monad2 A data-frame object.
#' @param method String, indicating type of missing data handling procesure. 
#' @param mat Numeric matrix.
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
#' @param pi1_mat,pi2_mat_tmp Numeric matrices (or NULL for pi2_mat_tmp); row-stochastic matrices of mixed-memberships. 
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
#' @param dyads,all.nodes1,all.nodes2 Arguments to internal function creating adjacency/affiliation matrices from data.frames
#' @param .soc_mats,Y,dyads,edges,t_id_d,t_id_n,nodes_pp,dyads_pp,nt_id,node_id_period,mu_b,var_b,n_dyads,n.blocks,periods,ctrl Internal arguments for initialization.
#' @param m1,bootrep,blist1,blist2 Arguments to internal fiunctions for bootstrapping 
#' @param ... Numeric vectors; vectors of potentially different length to be cbind-ed.
#' 
#' @author Santiago Olivella (olivella@@unc.edu), Adeline Lo (aylo@@wisc.edu), Tyler Pratt (tyler.pratt@@yale.edu), Kosuke Imai (imai@@harvard.edu)
#'
#' @return See individual return section for each function:
#' \describe{
#'       \item{.cbind.fill}{Matrix of \code{cbind}'ed elements in \code{...}, with missing values in each vector filled with \code{NA}.}
#'       \item{.mpower}{Matrix; the result of raising \code{mat} to the \code{p} power.}
#'       \item{.findPerm}{List of permuted blockmodel matrices.}
#'       \item{.transf}{Matrix with transformed mixed-membership vectors along its rows, s.t. no element is equal to 0.0 or 1.0.}
#'       \item{.compute.alpha}{List of predicted alpha matrices, one element per HMM state.}
#'       \item{.pi.hat}{Matrix of predicted mixed-membership vectors along its rows, with expectation computed over marginal 
#'                     distribution over HMM states for each time period.}
#'       \item{.missing}{Transformed data.frame with missing values list-wise deleted, or expanded
#'                       with missing indicator variables.}
#'       \item{.createSocioB}{List of sociomatrices.}
#'       \item{.vertboot2}{List of bootstrapped sociomatrices.}
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
.scaleVars <- function(x, keep_const = TRUE){
  A <- model.matrix(terms(x), x) #assumes intercept comes first
  which_cont <- which(apply(A, 2, function(a)length(unique(a))>2))
  A_sd <- rep(1, ncol(A))
  A_sd[which_cont] <- apply(A[,which_cont, drop = FALSE], 2, sd)*2
  A <- scale(A, TRUE, A_sd)
  constx <- which(colnames(A)=="(Intercept)")
  if(keep_const){
    A[,constx] <- 1
    attr(A, "scaled:scale") <- c(1, attr(A, "scaled:scale")[-constx])
    attr(A, "scaled:center") <- c(0, attr(A, "scaled:center")[-constx])
  } else {
    attr_tmp <- list(attr(A, "scaled:center")[-constx],
                     attr(A, "scaled:scale")[-constx])
    A <- A[,-constx, drop = FALSE]
    attr(A, "scaled:center") <- attr_tmp[[1]]
    attr(A, "scaled:scale") <- attr_tmp[[2]]
  }
  return(A)
}

#' @rdname auxfuns
.transf_muvar <- function(orig, is_var, is_array, des.mat, nblock=NULL, nstate=NULL){
  if(is_array){
    tmp <- array(ifelse(is_var, 1.0, 0.0), c(ncol(des.mat), nblock, nstate))
    rownames(tmp) <- colnames(des.mat)
  } else {
    tmp <- array(ifelse(is_var, 1.0, 0.0), ncol(des.mat))
    names(tmp) <- colnames(des.mat)
  }
  if(length(orig) > 1){
    if(is_array){
      tmp[rownames(orig),,] <- orig
    } else {
      tmp[names(orig)] <- orig
    }
  } else {
    non_miss <- !grepl("_missing", colnames(des.mat))
    if(is_array){
      tmp[non_miss,,] <- orig
    } else {
      tmp[non_miss] <- orig
    }
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
.vcovBeta <- function(all_phi, beta_coef, n.sim, n.blk, n.hmm, n.nodes, n.periods,
                      mu.beta, var.beta, est_kappa, t_id_n, X, fit){
  sampleC_perm <- do.call(rbind,
                          lapply(all_phi,
                         function(mat){
                           apply(mat, 2, function(vec)poisbinom::rpoisbinom(n.sim, vec))
                         })) 
  C_samples <- split.data.frame(sampleC_perm, rep(1:n.sim, times = length(all_phi)))
  S_samples <- replicate(n.sim, apply(est_kappa, 2, function(x)sample(1:n.hmm, 1, prob = x)), simplify = FALSE)
  hessBeta_list <- mapply(
    function(C_samp, S_samp, tidn, X_i, Nvec, beta_vec, vbeta, mbeta, periods)
    {
      if(n.hmm > 1) {
        s_matrix <- t(model.matrix(~factor(S_samp, 1:n.hmm) - 1))
      } else {
        s_matrix <- matrix(1, ncol=periods)
      }
      tot_in_state <- rowSums(s_matrix)
      if(any(tot_in_state == 0.0)){
        which_empty_s <- which(tot_in_state < 1.0)
        
        warning("Some HMM states are empty; no standard errors will be returned for coefficients associated with them.")
      }  
      hess_tmp <- optimHess(c(beta_vec),alphaLBound,alphaGrad,
                            tot_nodes = Nvec,
                            c_t = t(C_samp),
                            x_t = X_i,
                            s_mat = s_matrix,
                            t_id = tidn,
                            var_beta = vbeta,
                            mu_beta = mbeta)
      vc_tmp <- Matrix::forceSymmetric(solve(hess_tmp))
      ev <- eigen(vc_tmp)$value
      if(any(ev<0)){
        vc_tmp <- vc_tmp - diag(min(ev)-1e-4, ncol(vc_tmp))
      }
      ch_vc <- chol(vc_tmp)
      return(t(ch_vc) %*% ch_vc)
    },
    C_samples, S_samples,
    MoreArgs = list(tidn = t_id_n,
                    X_i = X,
                    Nvec = n.nodes,
                    beta_vec = beta_coef, 
                    vbeta = var.beta, 
                    mbeta = mu.beta,
                    periods = n.periods),
    SIMPLIFY=FALSE)
  vcov_monad <- Reduce("+", hessBeta_list)/n.sim
  
  colnames(vcov_monad) <- rownames(vcov_monad) <- paste(rep(paste("State",1:n.hmm), each = prod(dim(beta_coef)[1:2])), #beta_coef used to be fbeta_coef??
                                                        rep(colnames(beta_coef), each = nrow(beta_coef), times = n.hmm),#beta_coef used to be fbeta_coef??
                                                        rep(rownames(beta_coef), times = n.blk*n.hmm),
                                                        sep=":")
  return(vcov_monad)
}



#' @rdname auxfuns
.e.pi <- function(alpha_list, kappa, C_mat = NULL){
  if(is.null(C_mat)){
    pi_l <- alpha_list
  } else {
    pi_l <- lapply(alpha_list, function(x) x + t(C_mat))
  }
  if(is.null(dim(kappa))){
    return(proportions(pi_l[[1]], 2))
  } else {
    n_states <- nrow(kappa)
    pi.states <- lapply(1:n_states,
                        function(m){
                          pi_l[[m]] * rep(kappa[m,], each=nrow(pi_l[[m]])) 
                        })
    return(proportions(Reduce("+", pi.states), 2))
  }
}

#' @rdname auxfuns
.initPi <- function(soc_mats,
                    dyads,
                    edges,
                    nodes_pp,
                    dyads_pp,
                    n.blocks, periods, directed, ctrl){
  res <- vector("list", 2L)
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
        n_elem <- n.blocks[1] + 1
        #res <- RSpectra::eigs_sym(U, n_elem)
        res <- eigen(U)
        sel_val <- order(abs(res$values), decreasing = TRUE)[1:n_elem] 
        eta <- res$vectors[,sel_val] %*% diag(res$values[sel_val])
        target <- eta[,2:n_elem] / (eta[,1] + 1e-8)
        sig <- 1 - (res$values[n_elem] / (res$values[n.blocks[1]]))
        sig <- ifelse(is.finite(sig), sig, 0)
        if(abs(sig) > 0.1){
          target <- target[,1:(n_elem - 2), drop = FALSE]
        }
      } else {
        target <- U
      }
      if(nrow(unique(target)) > n.blocks[1]){
        clust_internal <- fitted(kmeans(x = target,
                                        centers = n.blocks[1],
                                        iter.max = 15,
                                        nstart = 10), "classes")
        
      } else {
        init_c <- sample(1:nrow(target), n.blocks[1], replace = FALSE)
        cents <- jitter(target[init_c, ])
        clust_internal <- fitted(suppressWarnings(kmeans(x = target,
                                                         centers = cents,
                                                         iter.max = 15,
                                                         algorithm = "Lloyd",
                                                         nstart = 1)), "classes")
      }
      
      phi_internal <- model.matrix(~ factor(clust_internal, 1:n.blocks[1]) - 1)
      phi_internal <- .transf(phi_internal)
      rownames(phi_internal) <- rownames(soc_mats[[i]])
      colnames(phi_internal) <- 1:n.blocks[1]
      MixedMembership <- t(phi_internal)
      int_dyad_id <- apply(dyads[[i]][,c("(sid)","(rid)")],
                           2,
                           function(x)match(x, colnames(MixedMembership)) - 1)
      BlockModel <- approxB(edges[[i]], int_dyad_id, MixedMembership)
      temp_res[[i]] <- list(BlockModel = BlockModel,
                            MixedMembership = MixedMembership)
      
    } else {
      n_prior <- (dyads_pp[i] - nodes_pp[i]) * .05
      a <- plogis(ctrl$mu_block) * n_prior
      b <- n_prior - a
      lda_beta_prior <- lapply(list(b,a),
                               function(prior){
                                 mat <- matrix(prior[2], n.blocks[1], n.blocks[1])
                                 diag(mat) <- prior[1]
                                 return(mat)
                               })
      ret <- lda::mmsb.collapsed.gibbs.sampler(network = soc_mats[[i]],
                                               K = n.blocks[1],
                                               num.iterations = 75L,
                                               burnin = 25L,
                                               alpha = ctrl$alpha,
                                               beta.prior = lda_beta_prior)
      MixedMembership <- prop.table(ret$document_expects, 2)
      colnames(MixedMembership) <- colnames(soc_mats[[i]])
      int_dyad_id <- apply(dyads[[i]][,c("(sid)","(rid)")], 2,
                           function(x)match(x, colnames(MixedMembership)) - 1)
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
  phis_temp <- lapply(temp_res, function(x)x$MixedMembership)
  phi.ord <- as.numeric(lapply(phis_temp, function(x)strsplit(colnames(x), "@")[[1]][2])) # to get correct temporal order
  mm_init_t <- do.call(cbind,mapply(function(phi,perm){perm %*% phi},
                                    phis_temp[order(phi.ord)], perms_temp, SIMPLIFY = FALSE))
  rownames(mm_init_t) <- 1:n.blocks[1]
  res[[1]] <- mm_init_t 
return(res)
}


