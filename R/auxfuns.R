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
#' @param pi_l List of mixed-membership matrices. 
#' @param kappa Numeric matrix; matrix of marginal HMM state probabilities.
#' @param y Numeric vector; vector of edge values.
#' @param d_id Integer matrix; two-column matrix with nr. dyads rows, containing zero-based
#'             sender (first column) and receiver (second column) node id's for each dyad. 
#' @param pi_mat Numeric matrix; row-stochastic matrix of mixed-memberships. 
#' @param colPalette A function produced by \code{colorRamp}.
#' @param range The range of values to label the legend.
#' @param par Vector of parameter values.
#' @param tot_nodes Integer vector; total number of nodes each node interacts with.
#' @param c_t Integer matrix; samples from Poisson-Binomial counts of a node instantiating a group.
#' @param x_t,z_t Numeric matrices; transposed monadic and dyadic design matrices.
#' @param s_mat Integer matrix; Samples of HMM states by time period. 
#' @param t_id Integer vector; for each node, what time-period is it observed in? zero-indexed.
#' @param mu_beta,var_beta Numeric arrays; prior mean and variances of monadic coefficients.
#' @param mu_gamma,var_gamma Numeric vectors; prior mean and variances of dyadic coefficients.
#' @param send_phi,rec_phi Numeric matrices; for each dyad, sampled group instantiated by sender and reciver in pair.
#' @param mu_b_t,var_b_t Numeirc matrices; prior mean and variance for blockmodel parameters.
#' @param directed Boolean; is the network directed?
#' @param orig Object to be transformed.
#' @param is_var Boolean. Is the object to be transformed a variance term?
#' @param is_array Boolean. Is the object to be transformed an array?
#' @param des.mat Numeric matrix. Design matrix corresponding to transformed object.
#' @param nblock Number of groups in model, defaults to \code{NULL}.
#' @param nstate Number of hidden Markov states in model, defaults to \code{NULL}.
#' @param devs Vector of standard deviations to use in transformation of variances. Defaults to \code{NULL}.
#' @param ... Numeric vectors; vectors of potentially different length to be cbind-ed.
#' 
#' @author Santiago Olivella (olivella@@unc.edu), Adeline Lo (adelinel@@princeton.edu), Tyler Pratt (tyler.pratt@@yale.edu), Kosuke Imai (imai@@harvard.edu)
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
.pi.hat <- function(X, beta){
  if(length(dim(beta))==2){
    mu <- exp(X %*% beta)
    pi.states <- list(t(mu))
    return(pi.states)
  }else{
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


## modularized fns from mmsbm.R

#Missing data
.missing <- function(formula.dyad, formula.monad1, formula.monad2=NULL, 
                     data.dyad, data.monad1, data.monad2=NULL, missing){
  
  if(missing=="indicator method"){
    # dyadic dataset
    if(length(all.vars(formula.dyad[[3]]))){
      miss.d <- apply(data.dyad[,all.vars(formula.dyad[[3]]), drop = FALSE], 2, function(x){length(na.omit(x))}) < nrow(data.dyad)
      md <- names(miss.d[miss.d])
      if(length(md)>0){
        m.ind <- apply(as.data.frame(data.dyad[,md]), 2, function(x){
          ifelse(is.na(x), 1, 0)
        })
        colnames(m.ind) <- paste(md, "_missing", sep="")
        data.dyad[,md] <- apply(as.data.frame(data.dyad[,md]), 2, function(x){
          x[is.na(x)] <- 0
          return(x)
        })
        data.dyad <- cbind(data.dyad, m.ind)
        fc <- paste(as.character(formula.dyad[[2]]), as.character(formula.dyad[[1]]),
                    paste(c(all.vars(formula.dyad)[-1], colnames(m.ind)), collapse=" + "))
        formula.dyad <- eval(parse(text=fc))
      }
    }
    # monadic1 dataset
    if(length(all.vars(formula.monad1[[2]]))){
      miss.m <- apply(data.monad1[,all.vars(formula.monad1[[2]]), drop = FALSE], 2, function(x){length(na.omit(x))}) < nrow(data.monad1)
      mm <- names(miss.m[miss.m])
      if(length(mm)>0){
        m.ind <- apply(as.data.frame(data.monad1[,mm]), 2, function(x){
          ifelse(is.na(x), 1, 0)
        })
        colnames(m.ind) <- paste(mm, "_missing", sep="")
        data.monad1[,mm] <- as.vector(apply(as.data.frame(data.monad1[,mm]), 2, function(x){
          x[is.na(x)] <- 0
          return(x)
        }))
        data.monad1 <- cbind(data.monad1, m.ind)
        fc <- paste("~", paste(c(all.vars(formula.monad1), colnames(m.ind)),  collapse=" + "))
        formula.monad1 <- eval(parse(text=fc))
      }
    }
    # monadic2 dataset
    if(length(all.vars(formula.monad2[[2]]))){
      miss.m <- apply(data.monad2[,all.vars(formula.monad2[[2]]), drop = FALSE], 2, function(x){length(na.omit(x))}) < nrow(data.monad2)
      mm <- names(miss.m[miss.m])
      if(length(mm)>0){
        m.ind <- apply(as.data.frame(data.monad2[,mm]), 2, function(x){
          ifelse(is.na(x), 1, 0)
        })
        colnames(m.ind) <- paste(mm, "_missing", sep="")
        data.monad2[,mm] <- as.vector(apply(as.data.frame(data.monad2[,mm]), 2, function(x){
          x[is.na(x)] <- 0
          return(x)
        }))
        data.monad2 <- cbind(data.monad2, m.ind)
        fc <- paste("~", paste(c(all.vars(formula.monad2), colnames(m.ind)),  collapse=" + "))
        formula.monad2 <- eval(parse(text=fc))
      }
    }
  }
  if(missing=="listwise deletion"){
    if(length(all.vars(formula.dyad[[3]]))){
      mdyad <- apply(data.dyad[,all.vars(formula.dyad[[3]])], 1, function(x){!any(is.na(x))})
    } else {
      mdyad <- TRUE
    }
    #monad1
    if(length(all.vars(formula.monad1[[2]]))){
      mmonad1 <- apply(data.monad1[,all.vars(formula.monad1[[2]])], 1, function(x){!any(is.na(x))})
    } else {
      mmonad1 <- TRUE
    }
    #monad2
    if(length(all.vars(formula.monad2[[2]]))){
      mmonad2 <- apply(data.monad2[,all.vars(formula.monad2[[2]])], 1, function(x){!any(is.na(x))})
    } else {
      mmonad2 <- TRUE
    }
    
    data.dyad <- data.dyad[mdyad,]
    data.monad1 <- data.monad1[mmonad1,]
    data.monad2 <- data.monad2[mmonad2,]
    d.keep <- lapply(unique(data.dyad[,timeID]), function(x){
      nts <- data.monad1[data.monad1[,timeID]==x,nodeID] #### How to change this to also include data.monad2?
      dd <- data.dyad[data.dyad[,timeID]==x,]
      dd <- dd[dd[,senderID] %in% nts & dd[,receiverID] %in% nts,]
      return(dd)
    })
    data.dyad <- do.call("rbind", d.keep)
  }
}

#Create sociomatrix for bipartite
.createSocioB <- function(dyads, all.nodes1, all.nodes2, directed){
  res <- lapply(dyads,
                function(dyad_df,
                         nnode1 = length(all.nodes1),
                         nnode2 = length(all.nodes2),
                         nodes1 = all.nodes1,
                         nodes2 = all.nodes2)
                  #if bipartite:
                  #need adj_mat to be composed of nnode1=length(all.nodes1), nnode2=length(all.nodes2), and dimnames are nodes1, nodes2
                  # then create two sets of node_names
                  #create function here that takes 1 and 2 always; then double up in input ifit is bipartite (outside this fn, so remove bipartite arg)
                {
                  indeces <- as.matrix(dyad_df[,c("(sid)","(rid)")])
                  time_ind <- unique(dyad_df[,"(tid)"])
                  adj_mat <-  matrix(NA,
                                     nnode1,
                                     nnode2,
                                     dimnames = list(nodes1,
                                                     nodes2))
                  adj_mat[indeces] <- dyad_df[,4] # out of bounds
                  if(!directed){
                    adj_mat[indeces[,c(2,1)]] <- dyad_df[,4]
                  }
                  diag(adj_mat) <- 0
                  adj_mat[is.na(adj_mat)] <- sample(0:1, sum(is.na(adj_mat)), replace = TRUE)
                  if(!directed){
                    mat_ind <- which(upper.tri(adj_mat), arr.ind = TRUE)
                    adj_mat[mat_ind[,c(2,1)]] <- adj_mat[upper.tri(adj_mat)]
                  }
                  node_names1 <- paste(nodes1, "@", time_ind, sep="")
                  node_names2 <- paste(nodes2, "@", time_ind, sep="")
                  dimnames(adj_mat) <- list(node_names1,
                                            node_names2)
                  return(adj_mat)
                })
  return(res)
}

#' @rdname auxfuns
.vertboot2 <- function(m1, boot_rep){#m1 <- igraph::as_adjacency_matrix(graph_ex) or m1 <- as.matrix(m1)
  res <- list()
  for (i in 1:boot_rep) {
    blist1 <- sample(0:(dim(m1)[1]-1), replace = TRUE) #blist is sampling or rows in C indexing
    blist2 <- sample(0:(dim(m1)[2]-1), replace = TRUE) 
    res <- c(res, list(vertboot_matrix_rcpp2(m1,blist1, blist2)))
  }
  res
}
