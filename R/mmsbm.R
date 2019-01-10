#' Dynamic mixed-membership stochastic blockmodel with covariates
#'
#' The function estimates a dynamic mixed-membership stochastic
#' blockmodel that incorporates covariates. 
#'
#' @param formula.dyad A \code{formula} object. The variable in \code{data.dyad} that contains 
#'     binary edges should be used as a LHS, and any dyadic predictors 
#'     can be included on the RHS (when no dyadic covariates are available, use \code{y ~ 1}).
#'     Same syntax as a \code{glm} formula. 
#' @param formula.monad An optional \code{formula} object. LHS is ignored. RHS contains 
#'     names of nodal atrributes found in \code{data.monad}.  
#' @param senderID Character string. Quoted name of the variable in \code{data.dyad} identifying 
#'     the sender node. For undirected networks, the variable simply contains name of first node 
#'     in dyad.
#' @param receiverID Character string. Quoted name of the variable in \code{data.dyad} identifying 
#'     the receiver node. For undirected networks, the variable simply contains name of second node 
#'     in dyad.
#' @param nodeID Character string. Quoted name of the variable in \code{data.monad} identifying 
#'     a node in either \code{data.dyad[,senderID]} or \code{data.dyad[,senderID]}. If not \code{NULL},
#'     every node \code{data.dyad[,senderID]} or \code{data.dyad[,senderID]} must be present in 
#'     \code{data.monad[,nodeID]}.
#' @param timeID Character string. Quoted name of the variable in both \code{data.dyad} and
#'     \code{data.monad} indicating the time in which network (and correspding nodal atrributes)
#'     were observed. The variable itself must be composed of integers.
#' @param data.dyad Data frame. Sociomatrix in ``long'' (i.e. dyadic) format. Must contain at
#'    least three variables: the sender identifier (or identifier of the first node in an undirected networks dyad),
#'    the receiver identifier (or identifier of the second node in an undirected network dyad), and the value
#'    of the edge between them. Currently, only edges between zero and one (inclusive) are supported.
#' @param data.monad Data frame. Nodal atributes. Must contain a node identifier matching the names of nodes
#'    used in the \code{data.dyad} data frame. 
#' @param n.groups Integer value. How many latent groups should be used to estimate the model?
#' @param n.hmmstates Integer value. How many hidden Markov state should be used in the HMM? Defaults 
#'    to 1 (i.e. no HMM).  
#' @param directed Boolean. Is the network directed? Defaults to \code{TRUE}.
#' @param missing Means of handling missing data. One of "indicator method" (default) or "listwise deletion".
#' @param mmsbm.control A named list of optional algorithm control parameters.
#'     \describe{
#'        \item{init}{Type of initialization algorithm for mixed-membership vectors. One of
#'                    \code{spectral} (default), \code{random}, or \code{lda} (see
#'                     \code{\link[lda:lda.collapsed.gibbs.sampler]{mmsb.collapsed.gibbs.sampler}} for details about this function.)}
#'        \item{lda_iter}{If \code{init="lda"}, number of MCMC iterations to obtain initial values. Defaults to 250}
#'        \item{lda_alpha}{If \code{init="lda"}, value of \code{alpha} hyperparameter. Defaults to 1}
#'        \item{max_em_iter}{Number of maximum iterations in variational EM. Defaults to 5e3}
#'        \item{max_opt_iter}{Number of maximum iterations of BFGS in M-step. Defaults to 10e3}
#'        \item{mu_b}{Numeric vector with two elements: prior mean of blockmodel's main diagonal elements, and
#'                    and prior mean of blockmodel's offdiagonal elements. Defaults to \code{c(5.0, -5.0)}}
#'        \item{var_b}{Numeric vector with two positive elements: prior variance of blockmodel's main diagonal elements, and
#'                    and prior variance of blockmodel's offdiagonal elements. Defaults to \code{c(1.0, 1.0)}}
#'        \item{var_beta}{Numeric positive value. (Gaussian) Prior variance of monadic coefficients. Defaults to 5.0.}
#'        \item{var_gamma}{Numeric positive value. (Gaussian) Prior variance of dyadic coefficients. Defaults to 5.0.}
#'        \item{var_xi}{Numeric positive value. (Gaussian) Prior variance of log-concentration paramemeter for mixed-
#'                      membership vector. Defaults to 1}
#'        \item{eta}{Numeric positive value. Concentration hyper-parameter for HMM. Defaults to 10.3}
#'        \item{threads}{Numeric integer. Number of available cores for paralellization. Defaults to 4}
#'        \item{conv_tol}{Numeric value. Absolute tolerance for VI convergence. Defaults to 1e-4}
#'        \item{verbose}{Boolean. Should extra information be printed as model iterates? Defaults to FALSE}
#'        }
#'       
#'     
#' @return Object of class \code{mmsbm}. List with named components:
#'     \describe{
#'       \item{MixedMembership}{Matrix of variational posterior of mean of mixed-membership vectors. \code{nodes} by \
#'                              \code{n.groups}}
#'       \item{PhiSend,PhiRec}{Matrices of estimated variational parameters for the group indicators. \code{n.groups}
#'                             by total number of observed dyads.}
#'       \item{MMConcentration}{Estimated concentration parameter of mixed-membership vectors}
#'       \item{BlockModel}{\code{n.groups} by \code{n.groups} matrix of estimated tie log-odds between members
#'                         of corresponding latent groups. The blockmodel.}
#'       \item{MonadCoef}{Array of estimated coefficient values for monadic covariates. Has \code{n.groups} columns,
#'                        and \code{n.hmmstates} slices.}
#'       \item{DyadCoef}{Vector estimated coefficient values for dyadic covariates}
#'       \item{TransitionKernel}{Matrix of estimated HMM transition probabilities}
#'       \item{Kappa}{Matrix of marginal probabilities of being in an HMM state at any given point in time. 
#'                    \code{n.hmmstates} by years (or whatever time interval networks are observed at)}
#'       \item{LowerBound}{Value of the lower bound at final iteration}
#'       \item{niter}{Final number of VI iterations}
#'       \item{monadic.data,dyadic.data,n_states,n_blocks, directed}{Original values of parameters used during estimation}
#'     }
#' @author Kosuke Imai (imai@@harvard.edu), Tyler Pratt (tyler.pratt@@yale.edu), Santiago Olivella (olivella@@unc.edu)
#' 
#' @example tests/Examples/MIDColdWar.R
#' 

mmsbm <- function(formula.dyad, formula.monad=~1, senderID, receiverID,
                  nodeID = NULL, timeID = NULL, data.dyad, data.monad = NULL,
                  n.groups, n.hmmstates = 1, directed = TRUE, missing = "indicator method",
                  mmsbm.control = list()){
  
  stopifnot(class(formula.dyad) == "formula",
            class(formula.monad) == "formula",
            is.data.frame(data.dyad),
            length(unique(data.monad[,timeID])) > n.hmmstates)
  if(!is.null(data.monad)){
    stopifnot(is.data.frame(data.monad),
              !is.null(nodeID))
  }
  
  
  
  ## Preprocess data
  
  ## Add time variable if only one period
  if(is.null(timeID)){
    timeID <- "tid"
    data.dyad[timeID] <- 1
    if(!is.null(data.monad)) {
      data.monad[timeID] <- 1  
    }
  }
  
  ## Address missing data 
  if(missing=="indicator method"){
    # dyadic dataset
    if(length(all.vars(formula.dyad[[3]]))){
      miss.d <- apply(data.dyad[,all.vars(formula.dyad[[3]])], 2, function(x){length(na.omit(x))}) < nrow(data.dyad)
      md <- names(miss.d[miss.d])
      if(length(md)>0){
        m.ind <- apply(data.dyad[,md], 2, function(x){
          ifelse(is.na(x), 1, 0)
        })
        colnames(m.ind) <- paste(colnames(m.ind), "_missing", sep="")
        data.dyad[,md] <- apply(data.dyad[,md], 2, function(x){
          x[is.na(x)] <- 0
          return(x)
        })
        data.dyad <- cbind(data.dyad, m.ind)
        fc <- paste(as.character(formula.dyad[[2]]), as.character(formula.dyad[[1]]),
                    paste(c(all.vars(formula.dyad)[-1], colnames(m.ind)), collapse=" + "))
        formula.dyad <- eval(parse(text=fc))
      }
    }
    # monadic dataset
    if(length(all.vars(formula.monad[[2]]))){
    miss.m <- apply(data.monad[,all.vars(formula.monad[[2]])], 2, function(x){length(na.omit(x))}) < nrow(data.monad)
    mm <- names(miss.m[miss.m])
    if(length(mm)>0){
      m.ind <- apply(as.data.frame(data.monad[,mm]), 2, function(x){
        ifelse(is.na(x), 1, 0)
      })
      colnames(m.ind) <- paste(mm, "_missing", sep="")
      data.monad[,mm] <- apply(as.data.frame(data.monad[,mm]), 2, function(x){
        x[is.na(x)] <- 0
        return(x)
      })
      data.monad <- cbind(data.monad, m.ind)
      fc <- paste("~", paste(c(all.vars(formula.monad), colnames(m.ind)),  collapse=" + "))
      formula.monad <- eval(parse(text=fc))
    }
    }
  }
  if(missing=="listwise deletion"){
    if(length(all.vars(formula.dyad[[3]]))){
      mdyad <- apply(data.dyad[,all.vars(formula.dyad[[3]])], 1, function(x){!any(is.na(x))})
    } else {
      mdyad <- TRUE
    }
    if(length(all.vars(data.monad[[2]]))){
      mmonad <- apply(data.monad[,all.vars(formula.monad[[2]])], 1, function(x){!any(is.na(x))})
    } else {
      mmonad <- TRUE
    }
    data.dyad <- data.dyad[mdyad,]
    data.monad <- data.monad[mmonad,]
    d.keep <- lapply(unique(data.dyad[,timeID]), function(x){
      nts <- data.monad[data.monad[,timeID]==x,nodeID]
      dd <- data.dyad[data.dyad[,timeID]==x,]
      dd <- dd[dd[,senderID] %in% nts & dd[,receiverID] %in% nts,]
      return(dd)
    })
    data.dyad <- do.call("rbind", d.keep)
  }
              
  
  ## Form dyadic model frame and model matrix
  dyadic <- do.call(model.frame, list(formula = formula.dyad,
                                      data = data.dyad,
                                      drop.unused.levels = TRUE,
                                      tid = as.name(timeID),
                                      sid = as.name(senderID),
                                      rid = as.name(receiverID)))
  dyadic_order <- with(dyadic, order(`(tid)`, `(sid)`, `(rid)`))
  time_order <- with(dyadic, order(unique(`(tid)`)))
  dyadic <- dyadic[dyadic_order, ]
  dyadic$dyad_id <- do.call(paste, c(dyadic[c("(sid)","(rid)")], sep = "->"))
  dyadic$send_time_id <- do.call(paste, c(dyadic[c("(sid)","(tid)")], sep = "@"))
  dyadic$rec_time_id <- do.call(paste, c(dyadic[c("(rid)","(tid)")], sep = "@"))
  dyadic$dyad_time_id <- do.call(paste, c(dyadic[c("send_time_id",
                                                   "rec_time_id")],
                                          sep = "->"))
  Z <- scale(model.matrix(formula.dyad, dyadic))
  Z_mean <- attr(Z, "scaled:center")
  Z_sd <- attr(Z, "scaled:scale")
  if(any(Z_sd==0)){
    constz <- which(Z_sd==0)
    Z <- as.matrix(Z[,-constz])
  }
  
  ##Form monadic model frame and model matrix
  if(!is.null(data.monad)){
    monadic <- do.call(model.frame, list(formula = formula.monad,
                                         data = data.monad,
                                         drop.unused.levels = TRUE,
                                         tid = as.name(timeID),
                                         nid = as.name(nodeID)))
  } else {
    monadic_s <- unique(dyadic[,c("(sid)","(tid)")])
    monadic_r <- unique(dyadic[,c("(rid)","(tid)")])
    names(monadic_s)[1] <- names(monadic_r)[1] <- "(nid)"
    monadic <- unique(rbind(monadic_s, monadic_r))
  }
  monadic_order <- with(monadic, order(`(tid)`, `(nid)`)) 
  monadic <- monadic[monadic_order, ]
  monadic$node_time_id <- do.call(paste, c(monadic[c("(nid)","(tid)")], sep="@"))
  monadic <- monadic[monadic$node_time_id %in% dyadic$send_time_id |
                       monadic$node_time_id %in% dyadic$rec_time_id,]
  
  ## Check that all nodes in dyadic data 
  ## are also in monadic data
  stopifnot(all(dyadic$send_time_id %in% monadic$node_time_id),
            all(dyadic$rec_time_id %in% monadic$node_time_id))
  
  X <- scale(model.matrix(formula.monad, monadic))
  X_mean <- attr(X, "scaled:center")
  X_sd <- attr(X, "scaled:scale")
  if(any(X_sd==0)){
    constx <- which(X_sd==0)
    if(length(constx)>1)
      stop("Multiple constants in matrix of monadic predictors.")
    X[,constx] <- 1
  }
  
  
  
  ## Form default control list
  ctrl <- list(blocks = n.groups,
               states = n.hmmstates,
               times = length(unique(dyadic[["(tid)"]])),
               directed = directed,
               phi_init_t = NULL,
               kappa_init_t = NULL,
               b_init_t = NULL,
               beta_init = NULL,
               gamma_init = NULL,
               init = "spectral", # should this be "kmeans"?
               lda_iter = 250,
               lda_alpha = 1,
               max_em_iter = 5000,
               max_opt_iter = 10e3,
               mu_b = c(5.0, -5.0),
               var_b = c(1.0, 1.0),
               var_beta = 1.0,
               var_gamma = 1.0,
               eta = 1,
               threads = 4,
               conv_tol = 1e-2,
               verbose = FALSE)
  ctrl[names(mmsbm.control)] <- mmsbm.control
  if(!is.null(ctrl$phi_init_t)&!all(grepl("@", colnames(ctrl$phi_init_t)))){
    stop("Initial values of phi_init_t must provide node names in columns, 
         in format 'node@time'.")
  }
  mu_b <- var_b <- array(NA, c(n.groups, n.groups))
  diag(mu_b) <- ctrl[["mu_b"]][1]
  mu_b[upper.tri(mu_b)|lower.tri(mu_b)] <- ctrl[["mu_b"]][2]
  diag(var_b) <- ctrl[["var_b"]][1]
  var_b[upper.tri(var_b)|lower.tri(var_b)] <- ctrl[["var_b"]][2]
  
  
  
  ## Create initial values
  if(ctrl$verbose){
    cat("Finding optimal starting values...\n")
  }
  
  ## Initial values for hidden Markov states
  if(is.null(ctrl$kappa_init_t)){
    if(n.hmmstates > 1){
      indeces <- as.matrix(dyadic[,c("(tid)","dyad_id")])  
      dyad_time <- matrix(NA,
                          length(unique(indeces[,1])),
                          length(unique(indeces[,2])),
                          dimnames = list(unique(indeces[,1]),
                                          unique(indeces[,2]))
      )
      dyad_time[indeces] <- model.response(dyadic)
      dyad_time[is.na(dyad_time)] <- sample(0:1, sum(is.na(dyad_time)), replace = TRUE)
      state_internal <- fitted(kmeans(dyad_time,
                                      n.hmmstates,
                                      nstart = 15),"classes")
      kappa_internal <- model.matrix(~ as.factor(state_internal) - 1)
      kappa_internal <- prop.table(kappa_internal + 1e-3, 1)
      ctrl$kappa_init_t <- t(kappa_internal)
    } else {
      ctrl$kappa_init_t <- matrix(1, 
                                  ncol = ctrl$times,
                                  dimnames = list(1, rep(1, ctrl$times)))
    }
  } else {
    ctrl$kappa_init_t <- matrix(ctrl$kappa_init_t[,time_order], ncol = ctrl$times)
  } 
  state_init <- apply(ctrl$kappa_init_t, 2, which.max)
  
  ## Initial values for mixed-membership vectors
  all.nodes <- unique(unlist(dyadic[,c("(sid)","(rid)")]))
  node.cols <- which(names(dyadic)%in%c("(sid)","(rid)", "(tid)"))
  dyads <- split.data.frame(dyadic[,c(node.cols, 1)], dyadic[, "(tid)"])
  soc_mats <- lapply(dyads,
                     function(dyad_df, 
                              nnode = length(all.nodes),
                              nodes = all.nodes)
                     {
                       indeces <- as.matrix(dyad_df[,c("(sid)","(rid)")])
                       adj_mat <-  matrix(NA, 
                                          nnode,
                                          nnode,
                                          dimnames = list(nodes,
                                                          nodes))
                       diag(adj_mat) <- 0
                       adj_mat[indeces] <- dyad_df[,4] # out of bounds
                       if(!directed){
                         adj_mat[indeces[,c(2,1)]] <- dyad_df[,4]
                       }
                       adj_mat[is.na(adj_mat)] <- sample(0:1, sum(is.na(adj_mat)), replace = TRUE)
                       if(!directed){
                         mat_ind <- which(upper.tri(adj_mat), arr.ind = TRUE)
                         adj_mat[mat_ind[,c(2,1)]] <- adj_mat[upper.tri(adj_mat)]
                       }
                       return(adj_mat)
                     })
  
  if(is.null(ctrl$phi_init_t)) {
    if(ctrl$init=="spectral"){
      phi_init  <- lapply(soc_mats,
                          function(W){
                              G <- W %*% t(W) + t(W) %*% W ##bibliometric symmetrization
                              degs  <- rowSums(G)
                            D <- diag(1/sqrt(degs+0.1*max(degs)))
                            L <- D %*% G %*% D 
                            res <- eigen(L, symmetric = TRUE)
                            ord_eigen <- order(abs(res$values), decreasing = TRUE)
                            eta_spectral <- res$vectors %*% diag(res$values)
                            X_eigen <- eta_spectral[,ord_eigen[2:n.groups]]/eta_spectral[,ord_eigen[1]] # divide by zero
                            d <- eta_spectral[,ord_eigen[1]] # temp
                            if(any(d==0)){d[d==0] <- -0.001} # temp
                            X_eigen <- eta_spectral[,ord_eigen[2:n.groups]]/d # temp
                            clust_internal <- fitted(kmeans(X_eigen,
                                                            n.groups,
                                                            nstart = 15),"classes")
                            phi_internal <- model.matrix(~ as.factor(clust_internal) - 1)
                            phi_internal <- prop.table(phi_internal + 1e-3, 1)
                            rownames(phi_internal) <- colnames(W)
                            return(t(phi_internal))
                          })
    } else if(ctrl$init=="random") {
      phi_init <- lapply(soc_mats,
                         function(mat){
                           rand_phis <- prop.table(matrix(runif(ncol(mat)*n.groups), c(n.groups, ncol(mat))), 2)
                           colnames(rand_phis) <- colnames(mat)
                           return(rand_phis)}
      )
    } else {
      phi_init <- lapply(soc_mats,
                         function(mat){
                           a_mat <- b_mat <- mu_b
                           a_mat[a_mat<0] <- b_mat[b_mat>0] <- 1
                           b_mat[b_mat<0] <- abs(b_mat[b_mat<0])
                           mat <- replace(mat, is.na(mat), 0)
                           phi_internal <- lda::mmsb.collapsed.gibbs.sampler(mat > 0,
                                                                             num.iterations = as.integer(ctrl$lda_iter),
                                                                             burnin = as.integer(ctrl$lda_iter/2),
                                                                             alpha = ctrl$lda_alpha,
                                                                             K = n.groups,
                                                                             beta.prior=list(a_mat,
                                                                                             b_mat)
                           )$document_expects
                           phi_internal <- prop.table(phi_internal, 2)
                           colnames(phi_internal) <- rownames(mat)
                           return(phi_internal)
                         })
    }
    blockmodel_temp <- mapply(approxB,
                              soc_mats,
                              phi_init,
                              MoreArgs = list(directed = directed),
                              SIMPLIFY = FALSE)
    if(is.null(ctrl$b_init_t)){
      right_perm <- .findPerm(blockmodel_temp, 100)
    } else {
      right_perm <- .findPerm(blockmodel_temp, 100, plogis(t(ctrl$b_init_t)))
    }
    phi_init_temp <- mapply(function(dyad,phi){
      ind <- match(unique(unlist(dyad[,c("(sid)","(rid)")])), 
                   colnames(phi))
      
      return(matrix(phi[,ind],
                    ncol=length(ind),
                    dimnames = list(rep(NA, n.groups),
                                    paste(colnames(phi)[ind], 
                                          dyad[1,"(tid)"],sep="@"))))
    }, dyad = dyads, phi = phi_init, SIMPLIFY = FALSE)
    
    ctrl$phi_init_t <- do.call(cbind, 
                               mapply(function(ind, mat){(mat[ind,])},
                                      right_perm,
                                      phi_init_temp,
                                      SIMPLIFY = FALSE))
  } 
  
  ## Define component node index in each dyad
  nodes_in_dyads <- cbind(match(dyadic[,"send_time_id"],
                                colnames(ctrl$phi_init_t)) - 1,
                          match(dyadic[,"rec_time_id"],
                                colnames(ctrl$phi_init_t)) - 1)
  
  ## Initial value of blockmodel
  if(is.null(ctrl$b_init_t)){
    ctrl$b_init_t <- qlogis(t(approxBdyad(model.response(dyadic),
                                          nodes_in_dyads,
                                          ctrl$phi_init_t,
                                          directed)))
    if(any(!is.finite(ctrl$b_init_t))){
      inf_ind <- which(!is.finite(ctrl$b_init_t))
      ctrl$b_init_t[inf_ind] <- sign(ctrl$b_init_t[inf_ind]) * 100
    }
  }
  
  print(ctrl$b_init_t)
  ## Initial value of monadic coefficients
  if(is.null(ctrl$beta_init)){
    ## Attach mm matrix to monadic df
    monadic$time_id_int <- match(monadic[,"(tid)"], unique(monadic[,"(tid)"])) 
    phi_df <- t(ctrl$phi_init_t)[match(colnames(ctrl$phi_init_t),
                                       monadic$node_time_id),]
    colnames(phi_df) <- paste("G",1:n.groups,sep="_")
    monadic <- cbind(monadic, phi_df)
    monadic_m <- split(monadic, state_init[monadic$time_id_int])
    alpha_par_init <- lapply(monadic_m,
                             function(dat){
                               n_row <- nrow(dat)
                               X_internal <-  model.matrix(formula.monad, data=dat)
                               Y_internal <- as.matrix(dat[,tail(names(dat), n.groups)])
                               #Y_transf <- as.matrix(log(Y_internal))
                               beta <- lm.fit(X_internal, Y_internal)$coefficients
                               return(list(beta=beta))
                             }) 
    
    ctrl$beta_init <- sapply(alpha_par_init,
                             function(x)cbind(x$beta),
                             simplify = "array")
  } else {
  ctrl$beta_init <- vapply(lapply(seq(dim(ctrl$beta_init)[3]), function(x) ctrl$beta_init[ , , x]), 
                           function(mat, sd_vec, mean_vec, X_mat){
                             mat <- matrix(mat, ncol=n.groups, nrow=ncol(X_mat))
                             constx <- which(sd_vec==0)
                             if(length(constx)!=0){
                               mat[constx, ] <- mat[constx, ] + mean_vec[-constx] %*% mat[-constx, ]
                             }
                             mat[-constx, ] <- mat[-constx, ] * sd_vec[-constx]
                             return(mat)
                           },
                           array(0.0, c(ncol(X), n.groups)),
                           sd_vec = X_sd,
                           mean_vec = X_mean,
                           X_mat = X)
  }
  ## Initial value of dyadic coefficients
  if(is.null(ctrl$gamma_init)){
    if(ncol(Z) > 0){
      ctrl$gamma_init <- glm.fit(Z, model.response(dyadic),family=binomial())$coefficients  
    } else {
      ctrl$gamma_init <- 0
    }
  } else {
    ctrl$gamma_init <- ctrl$gamma_init * Z_sd[-which(Z_sd==0)]
  }
  
  if(ncol(Z) == 0)
    Z <- matrix(0, nrow = nrow(Z), ncol = 1)
  
  ## Compute additional quantities needed by main est. routine
  t_id_d <- as.integer(match(dyadic[,c("(tid)")], unique(dyadic[,c("(tid)")])) - 1) 
  t_id_n <- as.integer(match(monadic[,c("(tid)")], unique(monadic[,c("(tid)")])) - 1)
  nodes_pp <- by(monadic$node_time_id, monadic[,"(tid)"], length)
  Y <- model.response(dyadic)
  ## Estimate model
  fit <- mmsbm_fit(t(Z),
                   t(X),
                   Y,
                   t_id_d,
                   t_id_n,
                   nodes_pp,
                   nodes_in_dyads,
                   mu_b,
                   var_b,
                   ctrl$phi_init_t,
                   ctrl$kappa_init_t,
                   ctrl$b_init_t,
                   c(ctrl$beta_init),
                   ctrl$gamma_init,
                   ctrl)
  
  
  
  ##Reorder and rename to match original inputs
  colnames(fit[["MixedMembership"]]) <- colnames(ctrl$phi_init_t)
  fit[["MixedMembership"]] <- fit[["MixedMembership"]][,order(monadic_order)]
  colnames(fit[["Kappa"]]) <- unique(monadic[,"(tid)"])
  fit[["Kappa"]] <- fit[["Kappa"]][,order(time_order), drop = FALSE]
  fit[["BlockModel"]] <- t(fit[["BlockModel"]])
  fit[["TransitionKernel"]] <- t(fit[["TransitionKernel"]])
  
  
  
  
  ## Rescale and name coefficients
  fit[["DyadCoef"]] <- fit[["DyadCoef"]] / Z_sd[-which(Z_sd==0)]
  if(length(fit[["DyadCoef"]])){
    fit[["BlockModel"]] <- fit[["BlockModel"]] - c(Z_mean[-constz] %*% fit[["DyadCoef"]])
  }
  dimnames(fit[["BlockModel"]]) <- replicate(2,paste("Group",1:n.groups), simplify = FALSE)
  dimnames(fit[["TransitionKernel"]]) <- replicate(2,paste("State",1:n.hmmstates), simplify = FALSE)
  
  if(ncol(Z)>1){
    names(fit[["DyadCoef"]]) <- colnames(Z) 
  }
  fit[["MonadCoef"]] <- vapply(fit[["MonadCoef"]],
                               function(mat, sd_vec, mean_vec){
                                 constx <- which(sd_vec==0)
                                 mat[-constx, ] <- mat[-constx, ] / sd_vec[-constx]
                                 if(length(constx)!=0)
                                   mat[constx, ] <- mat[constx, ] - mean_vec[-constx] %*% mat[-constx, ]
                                 return(mat)
                               },
                               array(0.0, c(ncol(X), n.groups)),
                               sd_vec = X_sd,
                               mean_vec = X_mean)
  rownames(fit[["MonadCoef"]]) <- colnames(X)
  colnames(fit[["MonadCoef"]]) <- paste("Group",1:n.groups)
  ## Include used data in original order
  fit$monadic.data <- monadic[order(monadic_order),]
  fit$dyadic.data <- dyadic[order(dyadic_order),]
  fit$Y <- Y[order(dyadic_order)]
  
  ## Include whether network is directed 
  fit$directed <- directed
  
  ## Include internal node id's
  fit$InternalNodeIndex <- nodes_in_dyads[order(dyadic_order),] + 1
  
  ## Include original call
  fit$call <- match.call()
  
  ##Assign class for methods
  class(fit) <- "mmsbm"
  
  return(fit)
}







