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
#' @param prob Numeric vector; predicted probabilities of edges.
#' @param pen Numeric vector; prior coefficient variaces.
#' 
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
#' @param alpha Numeric matrix; Mixed-membership parameter matrix.
#' @param alpha_sum Numeric vector; sums of rows of alpha (Nr. of nodes x 1). 
#' @param kappa Numeric matrix; matrix of marginal HMM state probabilities.
#' @param C_mat Numeric matrix; matrix of posterior counts of block instantiations per node. 
#' @param y Numeric vector; vector of edge values.
#' @param d_id Integer matrix; two-column matrix with nr. dyads rows, containing zero-based
#'             sender (first column) and receiver (second column) node id's for each dyad. 
#' @param pi_mat,pi1_mat,pi2_mat_tmp Numeric matrices (or NULL for pi2_mat_tmp); row-stochastic matrices of mixed-memberships. 
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
#' @param moretimes Run 5 times for each period's initialization?
#' @param realign Use BM realignment method?
#' @param fp5times Run the first period's initialization 5 times?
#' @param gibbs_iter,burn_in,n_cores Arguments for collapsed Gibb's sampler.
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
#'       \item{.collapsedGibbs}{Collasped Gibbs' sampler for MAP estimates}
#'     }
#' 
#' @rdname auxfuns
.cbind.fill<-function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  rnames <- unique(unlist(lapply(nm, rownames)))
  print("start line 80")
  mat <- matrix(NA, nrow = length(rnames), ncol = length(nm), 
                dimnames = list(rnames, 1:length(nm)))
  print("end line 80")
  for (x in seq_along(nm)) {
    mat[rownames(nm[[x]]), x] <- nm[[x]]
  }
  return(mat)
}

#' @rdname auxfuns
.monadData <- function(form, d, timeID, nodeID, dntid, verbose){
  mfm <- do.call(model.frame, list(formula = as.formula(form),
                                   data = d,
                                   drop.unused.levels = TRUE,
                                   tid = as.name(timeID),
                                   nid = as.name(nodeID)))
  mfm[,"(nid)"] <- as.character(mfm[,"(nid)"]) 
  if(anyDuplicated(mfm[,c("(tid)","(nid)")])){
    stop("timeID and nodeID do not uniquely identify observations in one of the monadic datasets.")
  }
  
  ntid <- do.call(paste, c(mfm[c("(nid)","(tid)")], sep="@"))
  if(!all(dntid %in% ntid))
    stop("Nodes in dyadic dataset missing from monadic dataset.
          Are node and time identifiers identical in data.dyad and the monadic datasets?")
  match_ids <- ntid %in% dntid
  if(any(!match_ids)){
    if(verbose){
      cat("\tSome nodes in one of the monadic datsets are not present in data.dyad; dropping them.\n")
    }
    mfm <- mfm[match_ids, ]
    ntid <- do.call(paste, c(mfm[c("(nid)","(tid)")], sep="@"))
  }
  mfm[,"(nid)"] <- as.character(mfm[,"(nid)"]) 
  
  return(list(mf = mfm,
              id = ntid))
}

#' @rdname auxfuns
.missHandle <- function(form, dat, method){
  data_terms <- terms.formula(form)
  var_names <- attr(data_terms, "term.labels")
  if(length(var_names)){
    ret <- list()
    if(identical(method, "indicator method")){
      miss.d <- apply(dat[,var_names, drop = FALSE], 2, anyNA)
      md <- names(miss.d[miss.d])
      if(length(md)>0){
        m.ind <- apply(dat[,md, drop = FALSE], 2, function(x){
          ifelse(is.na(x), 1, 0)
        })
        colnames(m.ind) <- paste(md, "_missing", sep="")
        dat[,md] <- apply(dat[,md, drop = FALSE], 2, function(x){
          x[is.na(x)] <- 0
          return(x)
        })
        ret$data <- cbind(dat, m.ind)
        fc <- as.formula(paste0("~.+",
                                paste(colnames(m.ind), collapse=" + ")))
        ret$form <- update.formula(form, fc)
      } else {
        ret <- list(form = form,
                    dat = dat)
      }
    } else { ## method = "listwise delete"
      keep_ind <- apply(as.matrix(dat[,var_names, drop = FALSE]), 1, function(x){!any(is.na(x))})
      ret$dat <- data[keep_ind,,drop = FALSE]
      ret$form <- form
    }
  } else {
    ret <- list(form = form,
                dat = dat)
  }
  return(ret)
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
    tmp <- array(ifelse(is_var, 5, 0.0), c(ncol(des.mat), nblock, nstate))
    rownames(tmp) <- colnames(des.mat)
    #tmp["(Intercept)",,] <- 5.0
  } else {
    tmp <- array(ifelse(is_var, 5, 0.0), ncol(des.mat))
    names(tmp) <- colnames(des.mat)
  }
  if(length(orig) > 1){
    if(is_array){
      tmp[rownames(tmp),,] <- orig
    } else {
      tmp[rownames(tmp)] <- orig
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
                                                          iteration = 100)$P)
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
#.vcovBeta <- function(beta_coef,tot_nodes, c_t, x_t,s_mat,t_id,var_beta,mu_beta){
# print("start tmp in line 326 auxfuns")
# tmp <- alphaLBound(c(beta_coef),
#                    tot_nodes,
#                    c_t,
#                    x_t,
#                    s_mat,
#                    t_id,
#                    var_beta,
#                     mu_beta)
#  print("end tmp in line 326 auxfuns")
# cat("tmp:", attr(tmp, "hessian"), "\n")
# cat("dim tmp:", dim(attr(tmp, "hessian")), "\n")
# vcov_monad <- as.matrix(Matrix::nearPD(solve(-attr(tmp, "hessian")))$mat)
# print("start colnames(vcov_monad)")
#  cat("ncol (cov_monad)", ncol(vcov_monad), "\n")
#  cat("nrow (cov_monad)", nrow(vcov_monad), "\n")
#  cat("cov_monad", vcov_monad, "\n")
#  cat("dim(beta_coef)",dim(beta_coef),"\n")
#  cat("(beta_coef)",beta_coef,"\n")
#  names_344<-paste(rep(paste("State",1:dim(beta_coef)[3]), each = prod(dim(beta_coef)[1:2])), #beta_coef used to be fbeta_coef??
#                                                       rep(colnames(beta_coef), each = nrow(beta_coef), times = dim(beta_coef)[3]),#beta_coef used to be fbeta_coef??
#                                                       rep(rownames(beta_coef), times = prod(dim(beta_coef)[2:3])),
#                                                       sep=":")
# cat("names_344",names_344,"\n")
# colnames(vcov_monad) <- rownames(vcov_monad) <- paste(rep(paste("State",1:dim(beta_coef)[3]), each = prod(dim(beta_coef)[1:2])), #beta_coef used to be fbeta_coef??
#                                                        rep(colnames(beta_coef), each = nrow(beta_coef), times = dim(beta_coef)[3]),#beta_coef used to be fbeta_coef??
#                                                       rep(rownames(beta_coef), times = prod(dim(beta_coef)[2:3])),
#                                                       sep=":")
# print("end colnames(vcov_monad)")
# return(as.matrix(vcov_monad))
#}

.vcovBeta <- function(c_mat, beta_coef, n.sim, n.blk, n.hmm, n.nodes, n.periods,
                      mu.beta, var.beta, est_kappa, t_id_n, X, tot_nodeid){
  #print(paste0("n.nodes is: ",n.nodes))
  hess_tmp <- optimHess(c(beta_coef),alphaLBound,alphaGrad,
                        tot_nodes = n.nodes,
                        c_t = t(c_mat),
                        x_t = t(X),
                        s_mat = est_kappa,
                        t_id = t_id_n,
                        var_beta = var.beta,
                        mu_beta = mu.beta)
 # print(paste0("hess_tmp is: ",hess_tmp))     
  #print(paste0("tot_nodeid is: ",tot_nodeid))    
  hess_tmp<-hess_tmp/tot_nodeid
 # print(paste0("nrow(c_mat): ",nrow(c_mat) ))
  vcov_monad <- Matrix::forceSymmetric(solve(hess_tmp))
  
  ev <- eigen(vcov_monad)$value
  if(any(ev < 0)){
    vcov_monad <- vcov_monad - diag(min(ev) - 1e-4, ncol(vcov_monad))
  }
  
  colnames(vcov_monad) <- rownames(vcov_monad) <- paste(rep(paste("State",1:n.hmm), each = prod(dim(beta_coef)[1:2])),
                                                        rep(colnames(beta_coef), each = nrow(beta_coef), times = n.hmm),
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
    return(proportions(pi_l[[1]],2))
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
.transfBeta <- function(coefs, n.hmmstates, mean_vec,  sd_vec, n.blk, cnames){
  res <- vapply(1:n.hmmstates,
                function(ind, coefs, sd_vec, mean_vec){
                  mat <- coefs[,,ind, drop=FALSE]
                  constx <- 1
                  mat[-constx, , 1] <- mat[-constx, , 1] / sd_vec[-constx]
                  if(length(constx)!=0){
                    mat[constx, ,1] <- mat[constx, ,1] - mean_vec[-constx] %*% mat[-constx, , 1]
                  }
                  return(mat)
                },
                array(0.0, c(nrow(coefs), n.blk)),
                coefs = coefs,
                sd_vec = sd_vec,
                mean_vec = mean_vec)
  rownames(res) <- cnames
  colnames(res) <- paste("Group", 1:n.blk)
  return(res)
}

#' @rdname auxfuns
.transfHess <- function(vcovmat, n.hmmstates, sd_vec, n.blk){
  # res <- #vapply(1:n.hmmstates,
  #      function(ind, vcovmat, sd_vec){
  #mat <- vcovmat[(ind-1)*(dim(vcovmat)[1]/n.hmmstates)+1:ind*(dim(vcovmat)[1]/n.hmmstates),(ind-1)*(dim(vcovmat)[1]/n.hmmstates)+1:ind*(dim(vcovmat)[1]/n.hmmstates)]
  constx <- 1
  sd_vec[constx]<-1
  sd_vec<-rep(rep(sd_vec, times = n.blk),times=n.hmmstates)
  #  mat[-constx, , 1] <- mat[-constx, , 1] / sd_vec[-constx]
  scaling_mat <- outer(sd_vec, sd_vec, `*`) 
  res <- vcovmat / scaling_mat  # Element-wise division
  #  if(length(constx)!=0){
  #  mat[constx, ,1] <- mat[constx, ,1] #- mean_vec[-constx] %*% mat[-constx, , 1]
  #  }
  #  return(mat)
  #},
  #array(0.0, c(nrow(coefs), n.blk)),
  #vcovmat = vcovmat,
  #sd_vec1 = sd_vec1,
  #sd_vec2 = sd_vec2)
  # rownames(res) <- cnames
  # colnames(res) <- paste("Group", 1:n.blk)
  return(res)
}

#' @rdname auxfuns
.initPi <- function(formula.monad,formula.dyad,soc_mats,
                    bipartite,
                    dyads,
                    edges,
                    nodes_pp,
                    realign,
                    moretimes,
                    fp5times,
                    dyads_pp,
                    n.blocks, periods, directed, ctrl,data.dyad,data.monad){
  res <- vector("list", 2L)
  init_lb<-list()
  init_niter<-list()
  init_lb_best<-list()
  init_niter_best<-list()
  init_bm<-list()
  init_bm_best<-list()
  init_seed<-list()
  init_seed_best<-list()
  init_distance<-list()
  init_distance_best<-list()

  realign<-TRUE #manual
  moretimes<-TRUE
  fp5times<-FALSE
  fp5times<-TRUE
  if(bipartite){
    if (periods==1){
      phi_init_temp <- lapply(soc_mats, function(mat){
        msg <- capture.output(clust.o<-blockcluster::coclusterBinary(mat,nbcocluster=c(n.blocks[1],n.blocks[2])))
        phi1_init_temp<-matrix(0,nrow=nrow(mat),ncol=n.blocks[1])
        phi2_init_temp<-matrix(0,nrow=ncol(mat),ncol=n.blocks[2])
        for(i in 1:nrow(mat)){#node1
          val<-clust.o@rowclass[i]+1
          phi1_init_temp[i,val]<-1
        }
        for(j in 1:ncol(mat)){#node2
          val<-clust.o@colclass[j]+1
          phi2_init_temp[j,val]<-1
        }
        res<-vector("list",length=2)
        res[[1]]<-t(phi1_init_temp)
        res[[2]]<-t(phi2_init_temp)
        return(res)
      }
      )
      res[[1]] <- do.call(cbind, lapply(phi_init_temp, `[[`, 1))#1st matrix of each element of big list
      res[[2]] <- do.call(cbind, lapply(phi_init_temp, `[[`, 2))#2nd matrix of each element of big list
    }else{
      #  init_lb<-c()
      # res[[1]]<-t(netSim$piS)
      # res[[2]]<-t(netSim$piB)
      out<-vector("list",length=periods)
      out2<-vector("list",length=periods)
      for (i in c(1)){
        cat("Now running year:", i, "\n")
        # dy<-netSim[["df_dyad_1"]]%>%filter(year==i)
        dy<-data.dyad%>%filter(year==i)
        sdf<-data.monad[[1]]%>%filter(year==i)
        bdf<-data.monad[[2]]%>%filter(year==i)

        seeds<-c(sample(100:9999, 1)) #run 5 times
        if(fp5times){
          seeds<-c(sample(100:9999,5))
        }
        else{
          seeds<-c(sample(100:9999, 1))
        }
        # seeds<-c(02138) #only run once

        best_model <- NULL
        best_lower_bound <- -Inf

        init_lb_i<-c()
        init_niter_i<-c()
        init_bm_i<-list()
        init_seed_i<-seeds 

        for (s in seeds){
          m_s<-mmsbm(formula.dyad = formula.dyad,
                    formula.monad = list(formula.monad[[1]], 
                                         formula.monad[[2]]),
                     timeID="year",
                     senderID = "id1",
                     receiverID = "id2",
                     nodeID = list("id","id"),
                     bipartite= TRUE,
                     data.dyad = dy,
                     data.monad = list(sdf,bdf),
                     n.blocks = c(n.blocks[1],n.blocks[2]), n.hmmstates = 1,
                     mmsbm.control = list(verbose = TRUE,
                                          threads=1,
                                          svi = TRUE,
                                          vi_iter = 10000,
                                          #   batch_size = 1.0,
                                          conv_tol = 1e-3,
                                          mu_gamma = ctrl[["mu_gamma"]],
                                          var_gamma = ctrl[["var_gamma"]],
                                          var_beta=list(ctrl[["var_beta"]][[1]][,,1],
                                                        ctrl[["var_beta"]][[2]][,,1]),
                                          #    mu_beta=list(ctrl[["mu_beta"]][[1]][,,1],
                                          #         ctrl[["mu_beta"]][[2]][,,1]),
                                          hessian = FALSE,
                                          seed=s))
          cat("Seed:", s, "\n")
          init_lb_i<-c(init_lb_i,m_s$LowerBound)
          init_niter_i<-c(init_niter_i,m_s$niter)
          init_bm_i<-append(init_bm_i,list(m_s$BlockModel))


          if (m_s$LowerBound > best_lower_bound) {
            best_lower_bound <- m_s$LowerBound
            best_model <- m_s
          }
        }
        m<-best_model
        bm_year1<-m$BlockModel
        cat("BM1:", m$BlockModel, "\n")
        init_lb_best<- append(init_lb,list(best_lower_bound))
        init_niter_best<-append(init_niter,list(m$niter))
        init_bm_best<-append(init_bm_best,list(m$BlockModel))
        init_seed_best<-append(init_seed_best,list(m$seed))

        init_lb<- append(init_lb,list(init_lb_i))
        init_niter<-append(init_niter,list(init_niter_i))
        init_bm<-append(init_bm,list(init_bm_i))
        init_seed<-append(init_seed,list(init_seed_i))



        #PredS = matrix(c(t(m$MixedMembership1)),nrow=2,byrow=T)
        #PredB = matrix(c(t(m$MixedMembership2)),nrow=2,byrow = T)
        PredS =m$MixedMembership1
        PredB=m$MixedMembership2

        out[[i]][[1]]<-PredS
        out[[i]][[2]]<-PredB

        out2[[i]]<-m$BlockModel
      }


      for (i in 2:periods){
        cat("Now running year:", i, "\n")
        dy<-data.dyad%>%filter(year==i)
        sdf<-data.monad[[1]]%>%filter(year==i)
        bdf<-data.monad[[2]]%>%filter(year==i)
        if (moretimes){
          seeds<-c(sample(100:9999, 5))} #run 5 times
        else{
          seeds<-c(sample(100:9999, 1))
        }
        # seeds<-c(02138) #only run once

        best_model <- NULL
        best_lower_bound <- -Inf

        init_lb_i<-c()
        init_niter_i<-c()
        init_bm_i<-list()
        init_seed_i<-seeds 
        m_s_list<-list()
        state_current<-ifelse(i<=25,1,2)

        for (s in seeds){
          m_s<-mmsbm(formula.dyad = formula.dyad,
                    formula.monad = list(formula.monad[[1]], 
                                         formula.monad[[2]]),
                     timeID="year",
                     senderID = "id1",
                     receiverID = "id2",
                     nodeID = list("id","id"),
                     bipartite= TRUE,
                     data.dyad = dy,
                     data.monad = list(sdf,bdf),
                     n.blocks = c(n.blocks[1],n.blocks[2]),  n.hmmstates = 1,
                     mmsbm.control = list(verbose = TRUE,
                                          threads=1,
                                          svi = TRUE,
                                          vi_iter = 10000,
                                          #  batch_size = 1.0,
                                          conv_tol = 1e-3,
                                          mu_gamma = ctrl[["mu_gamma"]],
                                          var_gamma = ctrl[["var_gamma"]],
                                          var_beta=list(ctrl[["var_beta"]][[1]][,,1],
                                                        ctrl[["var_beta"]][[2]][,,1]),
                                          #  mu_beta=list(ctrl[["mu_beta"]][[1]][,,state_current],
                                          #          ctrl[["mu_beta"]][[2]][,,state_current]),
                                          hessian = FALSE,
                                          seed=s))
          cat("Seed:", s, "\n")
          cat("BM original",i, m_s$BlockModel, "\n")
          cat("Current LB",i, m_s$LowerBound, "\n")


          init_lb_i<-c(init_lb_i,m_s$LowerBound)
          init_niter_i<-c(init_niter_i,m_s$niter)
          init_bm_i<-append(init_bm_i,list(m_s$BlockModel))
          m_s_list<-append(m_s_list,list(m_s))
          #    if (m_s$LowerBound > best_lower_bound) {
          #      best_lower_bound <- m_s$LowerBound
          #      best_model <- m_s
          #    }
        }
        # find the best initialization based on the permuted matrix's distance to bm in year 1

        find_best_init<-function(init_out){
          block_models <- init_out
          target_ind <- which.max(sapply(soc_mats, ncol))
          # perms_temp <- .findPerm(block_models, target_mat = block_models[[target_ind]], use_perms = ctrl$permute)
          bm1<-block_models #original bm
          bm1<-lapply(bm1,plogis)
          #  print(bm1)
          #   print(class(bm1))
          #  print(bm1[[1]])
          # print(unlist(bm1))
          # print(length(bm1))

          bm<-bm1[1:length(seeds)]
          permute_matrix <- function(mat) {
            index <- 1
            perms_temp_store<-list()
            permuted_matrices<-list()
            #   print(mat)
            m<-nrow(mat)
            n<-ncol(mat)
            #   print(m)
            #  print(n)
            all_row_perms <- gtools::permutations(m, m, v=1:m)
            all_col_perms <- gtools::permutations(n, n, v=1:n)
            for (i in 1:nrow(all_row_perms)) {
              for (j in 1:nrow(all_col_perms)) {
                # Create permutation matrices for the i-th row permutation and the j-th column permutation
                row_perm_matrix <- as.matrix(as(all_row_perms[i,], "pMatrix"))
                col_perm_matrix <- as.matrix(as(all_col_perms[j,], "pMatrix"))

                # Store the pair of permutation matrices in perms_temp_store
                perms_temp_store[[index]] <- list(row_perm_matrix, col_perm_matrix)
                index <- index + 1
              }
            }
            for (i in 1:length(perms_temp_store)){
              permuted_matrices[[i]]<-perms_temp_store[[i]][[1]]%*%mat%*%perms_temp_store[[i]][[2]]
            }

            return(permuted_matrices)
          }
          # bm2<-lapply(bm1, permute_matrix)
          bm2<-list()
          #   print(length(seeds))
          #  print("starting problematic chunck")
          #  print(length(seeds))
          #  print(length(bm1))
          for (i in 1:length(seeds)){
            bm2[[i]]<-permute_matrix(bm1[[i]])
          }
          #   print("end problematic chunck")
          #   print("finished bm2 step")
          #   print(bm2)

          bm_base<-plogis(bm_year1)
          #bm_base<-matrix((c(0.9, 0.2, 0.05, 0.35)), ncol = 2) #if want to use the truth

          # Define a function to find the closest matrix to bm_base in a list
          calculate_norm <- function(matrix1, matrix2) {
            return(base::norm(matrix1 - matrix2, type = "f"))
          }

          m<-nrow(bm_base)
          n<-ncol(bm_base)
          all_row_perms <- gtools::permutations(m, m, v=1:m)
          all_col_perms <- gtools::permutations(n, n, v=1:n)

          find_closest_matrix <- function(m, t_mat) {
            permuted_matrix_list <- m
            smallest_norm <- Inf
            perms_temp <- NULL
            #   cat("length(permuted_matrix_list)",length(permuted_matrix_list),"\n")
            # Loop through each matrix in the list
            for (i in 1:length(permuted_matrix_list)) {
              current_matrix <- permuted_matrix_list[[i]]
              current_norm <- calculate_norm(current_matrix, t_mat)

              # Check if the current norm is smaller than the smallest found so far
              if (current_norm < smallest_norm) {
                smallest_norm <- current_norm
                perms_temp_id <-  i 
              }
            }
            # Initialize perms_temp_store to hold pairs of permutation matrices
            perms_temp_store <- list()

            # Populate perms_temp_store with all possible combinations of row and column permutation matrices
            index <- 1
            for (i in 1:nrow(all_row_perms)) {
              for (j in 1:nrow(all_col_perms)) {
                # Create permutation matrices for the i-th row permutation and the j-th column permutation
                row_perm_matrix <- as.matrix(as(all_row_perms[i,], "pMatrix"))
                col_perm_matrix <- as.matrix(as(all_col_perms[j,], "pMatrix"))

                # Store the pair of permutation matrices in perms_temp_store
                perms_temp_store[[index]] <- list(row_perm_matrix, col_perm_matrix)
                index <- index + 1
              }
            }
            perms_temp<-permuted_matrix_list[[perms_temp_id]]

            #  cat("Permutation id",perms_temp_id,"\n")

            # cat("BM original",i, m$BlockModel, "\n")
            return(perms_temp)
          }

          #  print("start line 683")
          perms_temp <- list()
          for (j in 1:length(seeds)){
            perms_temp[[j]]<-find_closest_matrix(bm2[[j]],t_mat=bm_base)
          }
          #  print("finished line 683")
          # print(perms_temp)
          dist<-lapply(perms_temp,calculate_norm,matrix2=bm_base)
          best_init_id<-which.min(dist)
          best_init_bm<-bm1[[best_init_id]]
          #cat("all distance",i, dist, "\n")
          #  cat("best distance",i, perms_temp[[best_init_id]], "\n")

          return(list(best_init_bm,dist,best_init_id))
        }
        if(moretimes==TRUE){
          #        m<-best_model
          out3<-find_best_init(init_bm_i)
          #  print("out3:")
          # print(out3)
          bestid<-out3[[3]]}
        else{
          bestid<-1
        }

        m<-m_s_list[bestid] #best init model (before realignment)
        #  print(bestid)
        m<-m[[1]] #now a list object
        # print(m$niter)

        #    cat("BM original (best)",i, m$BlockModel, "\n")
        #  cat("Best LB",i, m$LowerBound, "\n")
        init_lb_best<- append(init_lb,list(best_lower_bound))
        init_niter_best<-append(init_niter,list(m$niter))
        init_bm_best<-append(init_bm_best,list(m$BlockModel))
        init_seed_best<-append(init_bm_best,list(m$seed))

        init_lb<- append(init_lb,list(init_lb_i))
        init_niter<-append(init_niter,list(init_niter_i))
        init_bm<-append(init_bm,list(init_bm_i))
        init_seed<-append(init_seed,list(init_seed_i))
        if(moretimes==TRUE){
          init_distance<-append(init_distance,list(out3[[2]]))
          init_distance_best<-append(init_distance_best,list(out3[[3]]))}

        #PredS = matrix(c(t(m$MixedMembership1)),nrow=2,byrow=T)
        #PredB = matrix(c(t(m$MixedMembership2)),nrow=2,byrow = T)
        PredS =m$MixedMembership1
        PredB=m$MixedMembership2

        out[[i]][[1]]<-PredS
        out[[i]][[2]]<-PredB

        out2[[i]]<-m$BlockModel
      }
      #res[[1]] <- do.call(cbind, lapply(out, `[[`, 1))#1st matrix of each element of big list
      #res[[2]] <- do.call(cbind, lapply(out, `[[`, 2))#2nd matrix of each element of big list

      #Node1
      block_models <- out2
      target_ind <- which.max(sapply(soc_mats, ncol))
      # perms_temp <- .findPerm(block_models, target_mat = block_models[[target_ind]], use_perms = ctrl$permute)
      bm1<-block_models #original bm
      bm1<-lapply(bm1,plogis)

      permute_matrix <- function(mat) {
        index <- 1
        perms_temp_store<-list()
        permuted_matrices<-list()
        m<-nrow(mat)
        n<-ncol(mat)
        all_row_perms <- gtools::permutations(m, m, v=1:m)
        all_col_perms <- gtools::permutations(n, n, v=1:n)
        for (i in 1:nrow(all_row_perms)) {
          for (j in 1:nrow(all_col_perms)) {
            # Create permutation matrices for the i-th row permutation and the j-th column permutation
            row_perm_matrix <- as.matrix(as(all_row_perms[i,], "pMatrix"))
            col_perm_matrix <- as.matrix(as(all_col_perms[j,], "pMatrix"))

            # Store the pair of permutation matrices in perms_temp_store
            perms_temp_store[[index]] <- list(row_perm_matrix, col_perm_matrix)
            index <- index + 1
          }
        }
        for (i in 1:length(perms_temp_store)){
          permuted_matrices[[i]]<-perms_temp_store[[i]][[1]]%*%mat%*%perms_temp_store[[i]][[2]]
        }

        return(permuted_matrices)
      }
      #   bm2<-lapply(bm1, permute_matrix)
      bm2<-list()
      #   print(length(seeds))
      #print("starting second problematic chunck")
      # print(length(seeds))
      # print(length(bm1))
      for (i in 1:periods){
        bm2[[i]]<-permute_matrix(bm1[[i]])
      }  
      #   print("end second problematic chunck")

      bm_base<-plogis(block_models[[1]])
      # Create lists to store outputs
perms_temp <- vector("list", length = periods)
best_matrix_list <- list()

# Identity permutation for year 1
perms_temp[[1]] <- list(diag(n.blocks[1]), diag(n.blocks[2]))
best_matrix_list[[1]] <- bm1[[1]]  # Year 1 unpermuted reference

      #bm_base<-matrix((c(0.9, 0.2, 0.05, 0.35)), ncol = 2) #if want to use the truth

      # Define a function to find the closest matrix to bm_base in a list
      calculate_norm <- function(matrix1, matrix2) {
        return(base::norm(matrix1 - matrix2, type = "f"))
      }

      m<-nrow(bm_base)
      n<-ncol(bm_base)
      all_row_perms <- gtools::permutations(m, m, v=1:m)
      all_col_perms <- gtools::permutations(n, n, v=1:n)

      find_closest_matrix <- function(m, t_mat) {
        permuted_matrix_list <- permute_matrix(m)
        smallest_norm <- Inf
        best_matrix<-NULL
        perms_temp <- NULL
        cat("length(permuted_matrix_list)",length(permuted_matrix_list),"\n")
        # Loop through each matrix in the list
        for (i in 1:length(permuted_matrix_list)) {
          current_matrix <- permuted_matrix_list[[i]]
          current_norm <- calculate_norm(current_matrix, t_mat)

          # Check if the current norm is smaller than the smallest found so far
          if (current_norm < smallest_norm) {
            smallest_norm <- current_norm
            best_matrix <- current_matrix
            perms_temp_id <-  i 
          }
        }
        # Initialize perms_temp_store to hold pairs of permutation matrices
        perms_temp_store <- list()

        # Populate perms_temp_store with all possible combinations of row and column permutation matrices
        index <- 1
        for (i in 1:nrow(all_row_perms)) {
          for (j in 1:nrow(all_col_perms)) {
            # Create permutation matrices for the i-th row permutation and the j-th column permutation
            row_perm_matrix <- as.matrix(as(all_row_perms[i,], "pMatrix"))
            col_perm_matrix <- as.matrix(as(all_col_perms[j,], "pMatrix"))

            # Store the pair of permutation matrices in perms_temp_store
            perms_temp_store[[index]] <- list(row_perm_matrix, col_perm_matrix)
            index <- index + 1
          }
        }
        perms_temp<-perms_temp_store[[perms_temp_id]]

        cat("Permutation id",perms_temp_id,"\n")

        # cat("BM original",i, m$BlockModel, "\n")
        return(list(perms_temp = perms_temp, best_matrix = best_matrix))
      }


  # Cumulative alignment loop
for (i in 2:periods) {
  bm_base <- Reduce("+", best_matrix_list[1:(i - 1)]) / (i - 1)

  result <- find_closest_matrix(bm1[[i]], t_mat = bm_base)
  perms_temp[[i]] <- result$perms_temp
  best_matrix_list[[i]] <- result$best_matrix
}

# ---- Realign Year 1 to the Average of Years 2 to T ----
bm_base_new <- Reduce("+", best_matrix_list[2:periods]) / (periods - 1)

result_first <- find_closest_matrix(bm1[[1]], t_mat = bm_base_new)

# Update the permutation and aligned matrix for Year 1
perms_temp[[1]] <- result_first$perms_temp
best_matrix_list[[1]] <- result_first$best_matrix
   
    #  if(realign){
    #    perms_temp <- lapply(bm1, find_closest_matrix, t_mat = bm_base)
   #   }
   #   else{
    #    perms_temp<-lapply(1:periods, function(x) list(as.matrix(as(all_perms[1,], "pMatrix")), as.matrix(as(all_perms[1,], "pMatrix"))))
    #  }

      phis_temp <- lapply(out, `[[`, 1) 
      perms_temp1<-lapply(perms_temp, `[[`, 1) 

      phi.ord <- as.numeric(lapply(phis_temp, function(x)strsplit(colnames(x), "@")[[1]][2])) # to get correct temporal order
      perms_temp1[[1]] <- diag(n.blocks[1])
      mm_init_t1 <- do.call(cbind,mapply(function(phi,perm){perm %*% phi },
                                         phis_temp, perms_temp1, SIMPLIFY = FALSE))

      #cat("dimension of mm_init_t1",dim(mm_init_t1),"\n")
      #cat("n.blocks[1]",n.blocks[1],"\n")
      rownames(mm_init_t1) <- 1:n.blocks[1]
      res[[1]] <- mm_init_t1

      #Node2
      #block_models <- out2
      #target_ind <- which.max(sapply(soc_mats, ncol))
      #perms_temp <- .findPerm(block_models, target_mat = block_models[[target_ind]], use_perms = ctrl$permute)
      phis_temp <- lapply(out, `[[`, 2) 
      perms_temp2<-lapply(perms_temp, `[[`, 2) 
      # p1<-lapply(perms_temp, `[[`, 1) 
      #  p2<-lapply(perms_temp, `[[`, 2) 
      phi.ord <- as.numeric(lapply(phis_temp, function(x)strsplit(colnames(x), "@")[[1]][2])) # to get correct temporal order
      perms_temp2[[1]] <- diag(n.blocks[2])
      #mm_init_t2 <- do.call(cbind,mapply(function(phi,perm){perm %*% phi },
      #                                  phis_temp, perms_temp2, SIMPLIFY = FALSE))
      mm_init_t2 <- do.call(cbind,mapply(function(phi,perm){t(t(phi) %*% perm ) },
                                         phis_temp, (perms_temp2), SIMPLIFY = FALSE))


      # cat("dimension of mm_init_t2",dim(mm_init_t2),"\n")
      #  cat("n.blocks[2]",n.blocks[2],"\n")
      rownames(mm_init_t2) <- 1:n.blocks[2]
      res[[2]] <- mm_init_t2
    }


  } else {
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
  }
  return(list(res,init_lb,init_niter,init_lb_best,init_niter_best,init_bm,init_bm_best,init_seed,init_seed_best,init_distance,init_distance_best))
  # return(res)
} 




#' @rdname auxfuns
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
                  print("starting line 1002")
                  adj_mat <-  matrix(NA,
                                     nnode1,
                                     nnode2,
                                     dimnames = list(nodes1,
                                                     nodes2))
                  print("end line 1002")
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
                  print("starting line 1021")
                  dimnames(adj_mat) <- list(node_names1,
                                            node_names2)
                  print("end line 1021")
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

#' @rdname auxfuns
.collapsedGibbs <- function(model.obj, gibbs_iter, burn_in, n_cores) {
  dyads <- model.obj$dyadic.data
  MM1 <- model.obj$MixedMembership1
  MM2 <- model.obj$MixedMembership2
  
  block_mat <- model.obj$BlockModel
  gamma <- model.obj$DyadCoef
  
  K1 <- nrow(MM1)
  K2 <- nrow(MM2)
  
  dyads$tid <- as.integer(dyads[["(tid)"]])
  dyads$sid <- dyads[["(sid)"]]
  dyads$rid <- dyads[["(rid)"]]
  
  all_t <- sort(unique(dyads$tid))
  
  logit_inv <- function(x) 1 / (1 + exp(-x))
  
  # Parallelized Gibbs for each time period
  time_results <- mclapply(all_t, function(t) {
    dyads_t <- subset(dyads, tid == t)
    N_dyads <- nrow(dyads_t)
    
    senders <- sort(unique(dyads_t$sid))
    receivers <- sort(unique(dyads_t$rid))
    
    sid_map <- setNames(seq_along(senders) - 1, senders)
    rid_map <- setNames(seq_along(receivers) - 1, receivers)
    
    dyads_t$sid_idx <- sid_map[dyads_t$sid]
    dyads_t$rid_idx <- rid_map[dyads_t$rid]
    
    alpha1_t <- sapply(senders, function(s) MM1[, paste0(s, "@", t)])
    alpha2_t <- sapply(receivers, function(r) MM2[, paste0(r, "@", t)])
    
    alpha1_t <- t(alpha1_t)
    alpha2_t <- t(alpha2_t)
    
    y <- dyads_t$Y
    dyad_covars <- setdiff(names(dyads_t), c("Y", "sid", "rid", "tid", "(sid)", "(rid)", "(tid)", "sid_idx", "rid_idx"))
    D <- as.matrix(dyads_t[, dyad_covars, drop = FALSE])
    d_id <- cbind(dyads_t$sid_idx, dyads_t$rid_idx)

    if (ncol(D) == 0 || length(gamma) == 0) {
    D <- matrix(0, nrow = nrow(D), ncol = 1)
    gamma <- 0
    }
    
    z <- sample(1:K1, N_dyads, replace = TRUE)
    u <- sample(1:K2, N_dyads, replace = TRUE)
    z_counts <- matrix(0, N_dyads, K1)
    u_counts <- matrix(0, N_dyads, K2)
    
    for (iter in 1:gibbs_iter) {
      for (i in 1:N_dyads) {
        p <- d_id[i, 1] + 1
        q <- d_id[i, 2] + 1
        
        # z sampling
        log_pz <- sapply(1:K1, function(g) {
          eta <- block_mat[g, u[i]] + D[i, ] %*% gamma
          theta <- logit_inv(eta)
          log(alpha1_t[p, g]) + y[i] * log(theta) + (1 - y[i]) * log(1 - theta)
        })
        prob_z <- exp(log_pz - max(log_pz)); prob_z <- prob_z / sum(prob_z)
        z[i] <- sample(1:K1, 1, prob = prob_z)
        
        # u sampling
        log_pu <- sapply(1:K2, function(h) {
          eta <- block_mat[z[i], h] + D[i, ] %*% gamma
          theta <- logit_inv(eta)
          log(alpha2_t[q, h]) + y[i] * log(theta) + (1 - y[i]) * log(1 - theta)
        })
        prob_u <- exp(log_pu - max(log_pu)); prob_u <- prob_u / sum(prob_u)
        u[i] <- sample(1:K2, 1, prob = prob_u)
        
        if (iter > burn_in) {
          z_counts[i, z[i]] <- z_counts[i, z[i]] + 1
          u_counts[i, u[i]] <- u_counts[i, u[i]] + 1
        }
      }
    }
    
    z_map <- apply(z_counts, 1, which.max)
    u_map <- apply(u_counts, 1, which.max)
    
    C_mat1 <- matrix(0, nrow(alpha1_t), K1)
    C_mat2 <- matrix(0, nrow(alpha2_t), K2)
    
    for (i in 1:N_dyads) {
      p <- d_id[i, 1] + 1
      q <- d_id[i, 2] + 1
      C_mat1[p, z_map[i]] <- C_mat1[p, z_map[i]] + 1
      C_mat2[q, u_map[i]] <- C_mat2[q, u_map[i]] + 1
    }
    
    return(list(C_mat1 = C_mat1, C_mat2 = C_mat2, z_map = z_map, u_map = u_map))
  }, mc.cores = n_cores)
  
  return(list(
    C_mat1_list = lapply(time_results, `[[`, "C_mat1"),
    C_mat2_list = lapply(time_results, `[[`, "C_mat2"),
    z_map_list  = lapply(time_results, `[[`, "z_map"),
    u_map_list  = lapply(time_results, `[[`, "u_map")
  ))
}
