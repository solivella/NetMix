################################################################
## Wrapper function to process data and initialize parameters ##
################################################################

mmsbm <- function(formula.dyad, formula.monad=~1, senderID, receiverID,
                  nodeID = NULL, timeID = NULL, data.dyad, data.monad = NULL,
                  n.blocks, n.hmmstates, directed = TRUE,
                  mmsbm.control = list()){

  stopifnot(class(formula.dyad) == "formula",
            class(formula.monad) == "formula",
            is.data.frame(data.dyad))
  if(!is.null(data.monad)){
    stopifnot(is.data.frame(data.monad),
              !is.null(nodeID))
  }
  
  
  ## Preprocess data to pass to mmsbm
  if(is.null(timeID)){
    timeID <- "tid"
    data.dyad[timeID] <- 1
  }
  mfd <- do.call(model.frame, list(formula = formula.dyad,
                                   data = data.dyad,
                                   drop.unused.levels = TRUE,
                                   tid = as.name(timeID),
                                   sid = as.name(senderID),
                                   rid = as.name(receiverID)))
  dyadic_order <- with(mfd, order(`(tid)`, `(sid)`, `(rid)`))
  time_order <- with(mfd, order(unique(`(tid)`)))
  mfd <- mfd[dyadic_order, ]
  
  
  ut <- unique(mfd[["(tid)"]])
  periods <- length(ut)
  dntid <- cbind(do.call(paste, c(mfd[c("(sid)","(tid)")], sep = "@")),
                 do.call(paste, c(mfd[c("(rid)","(tid)")], sep = "@")))
  udntid <- unique(c(dntid))
  udnid <- unique(unlist(mfd[c("(sid)","(rid)")]))
  
  if(is.null(data.monad)){
    data.monad <- data.frame(nid = rep(udnid, periods))
    nodeID <- "nid"
    data.monad[timeID] <- rep(ut, each = length(udnid))
  }
  mfm <- do.call(model.frame, list(formula = formula.monad,
                                   data = data.monad,
                                   drop.unused.levels = TRUE,
                                   tid = as.name(timeID),
                                   nid = as.name(nodeID)))
  monadic_order <- with(mfm, order(`(tid)`, `(nid)`)) 
  mfm <- mfm[monadic_order, ]
  ntid <- do.call(paste, c(mfm[c("(nid)","(tid)")], sep="@"))
  mfm <- mfm[ntid %in% unique(c(dntid)), ]
  if(!all(udntid %in% ntid))
    stop("Nodes in dyadic dataset missing from monadic dataset. Are node and time identifiers identical in data.dyad and data.monad?")
  
  
  Y <- model.response(mfd)
  X <- scale(model.matrix(terms(mfm), mfm))
  X_mean <- attr(X, "scaled:center")
  X_sd <- attr(X, "scaled:scale")
  if(any(X_sd==0)){
    constx <- which(X_sd==0)
    if(length(constx)>1)
      stop("Singularities in matrix of monadic predictors.")
    X[,constx] <- 1
  }
  n_monad_pred <- ncol(X)
  Z <- scale(model.matrix(terms(mfd), mfd))
  Z_mean <- attr(Z, "scaled:center")
  Z_sd <- attr(Z, "scaled:scale")
  if(any(Z_sd==0)){
    constz <- which(Z_sd==0)
    Z <- Z[,-constz]
  }
  n_dyad_pred <- ncol(Z)
  
  nt_id <- cbind(.mapID(udntid, dntid[, 1]) - 1, .mapID(udntid, dntid[, 2]) - 1)
  t_id_d <- .mapID(ut, mfd[["(tid)"]]) - 1
  t_id_n <- .mapID(ut, mfm[["(tid)"]]) - 1
  nodes_pp <- c(by(mfm, mfm[["(tid)"]], nrow))
  
  ## Form default control list
  ctrl <- list(blocks = n.blocks,
               states = n.hmmstates,
               times = periods,
               directed = directed,
               phi_init_t = NULL,
               kappa_init_t = NULL,
               b_init_t = NULL,
               beta_init = NULL,
               gamma_init = NULL,
               spectral = TRUE,
               em_iter = 500,
               opt_iter = 10e3,
               mu_b = c(1.0, 0.0),
               var_b = c(1.0, 1.0),
               var_beta = 1.0,
               var_gamma = 1.0,
               eta = 1.3,
               threads = 2,
               conv_tol = 1e-4,
               verbose = FALSE)
  ctrl[names(mmsbm.control)] <- mmsbm.control
  mu_b <- var_b <- array(NA, c(n.blocks, n.blocks))
  diag(mu_b) <- ctrl[["mu_b"]][1]
  mu_b[upper.tri(mu_b)|lower.tri(mu_b)] <- ctrl[["mu_b"]][2]
  diag(var_b) <- ctrl[["var_b"]][1]
  var_b[upper.tri(var_b)|lower.tri(var_b)] <- ctrl[["var_b"]][2]
  
  
  
  ## Create initial values
  dyads <- split.data.frame(dntid, mfd[, "(tid)"])
  edges <- split(Y, mfd[, "(tid)"])
  soc_mats <- mapply(function(mat, y){
    nnode <- length(unique(c(mat)))
    ifelse(directed,
           adj_mat <- matrix(NA, nnode, nnode,
                             dimnames = list(unique(mat[,1]),
                                            unique(mat[,2]))),
           adj_mat <- matrix(NA, nnode, nnode,
                             dimnames = list(unique(c(mat)),
                                             unique(c(mat)))))
    adj_mat[mat] <- y
    if(!directed)
      adj_mat[mat[,c(2,1)]] <- y
    return(adj_mat)
  },
  dyads, edges,
  SIMPLIFY = FALSE)
  
  td_id <- cbind(mfd[,"(tid)"],paste(mfd[,"(sid)"],mfd[,"(rid)"], sep = "->"))
  dyad_time <- matrix(NA, periods, length(unique(td_id[,2])),
                      dimnames = list(ut,
                                      unique(td_id[,2])))
  dyad_time[td_id] <- Y
  if(is.null(ctrl$kappa_init_t)){
    if(n.hmmstates > 1){
      state_internal <- cluster::clara(dyad_time,
                                       n.hmmstates,
                                       samples = 5,
                                       correct.d = TRUE)$clustering
      kappa_internal <- model.matrix(~ as.factor(state_internal) - 1)
      kappa_internal <- prop.table(kappa_internal+runif(length(kappa_internal)), 1)
      ctrl$kappa_init_t <- t(kappa_internal)
    } else {
      ctrl$kappa_init_t <- t(matrix(1, nrow = length(ut)))
    }
  } else {
    ctrl$kappa_init_t <- matrix(ctrl$kappa_init_t[,time_order], ncol = periods)
  } 
  if(n.hmmstates==1){
    state_init <- matrix(1, ncol = periods)
    names(state_init) <- rep(1, periods)
  } else {
    state_init <- apply(ctrl$kappa_init_t, 2, which.max)
  }
  
  if(is.null(ctrl$phi_init)) {
    phi_init_temp <- lapply(soc_mats,
                            function(mat){
                              if(ctrl$spectral) {
                                mat[is.na(mat)] <- sample(1, sum(is.na(mat)), replace = TRUE)
                                sv <- svd(mat)
                                d <- matrix(0, nrow(mat), ncol(mat))
                                diag(d) <- sv$d
                                target <- d %*% t(sv$v)
                              } else {
                                target <- mat
                              }
                              clust_internal <- cluster::clara(target,
                                                               n.blocks,
                                                               samples = 5,
                                                               correct.d = TRUE)$clustering
                              phi_internal <- model.matrix(~ as.factor(clust_internal) - 1)
                              phi_internal <- prop.table(phi_internal + runif(length(phi_internal)), 1)
                              rownames(phi_internal) <- rownames(mat)
                              t(phi_internal)
                            })
    dyad_internal <- mapply(function(mat1, mat2)
    {
      cbind(match(mat1[,1], colnames(mat2)) - 1,
            match(mat1[,2], colnames(mat2)) - 1) 
    }, dyads, phi_init_temp, SIMPLIFY = FALSE)
    blockmodel_temp <- mapply(approxB,
                              edges,
                              dyad_internal,
                              phi_init_temp,
                              SIMPLIFY = FALSE)
    if(is.null(ctrl$b_init_t)){
      right_perm <- .findPerm(blockmodel_temp, 100)
    } else {
      right_perm <- .findPerm(blockmodel_temp, 100, plogis(t(ctrl$b_init_t)))
    }
    ctrl$phi_init_t <- do.call(cbind,
                               mapply(function(perm_mat, phi_mat){perm_mat%*%phi_mat},
                              right_perm,
                              phi_init_temp,
                              SIMPLIFY = FALSE))
  } else {
    ctrl$phi_init_t <- ctrl$phi_init_t[, monadic_order]
  }
  
  if(is.null(ctrl$b_init_t)){
    ctrl$b_init_t <- qlogis(approxB(Y, nt_id, ctrl$phi_init_t))
  }
  
  
  if(is.null(ctrl$beta_init)){
    ctrl$beta_init <- sapply(1:n.hmmstates,
                             function(m, dm, df, phi_i, states){
                               dinm <- with(df, `(tid)` %in% unique(`(tid)`)[states == m])
                               phi_temp <- t(phi_i[,dinm])
                               X_sub <- dm[dinm, , drop = FALSE]
                               ytemp <- suppressWarnings(DirichletReg::DR_data(phi_temp))
                               do.call(cbind,coef(suppressWarnings(DirichletReg::DirichReg(ytemp~.-1, as.data.frame(X_sub)))))
                             },
                             dm = X, df = mfm, phi_i = ctrl$phi_init_t, states = state_init,
                             simplify = "array")
  }
  
  if(is.null(ctrl$gamma_init)){
    ctrl$gamma_init <- if(ncol(Z) > 0){
      coef(glm(Y ~ Z - 1, family = binomial()))
    } else {
      0
    }
  }
  
  if(ncol(Z) == 0)
    Z <- matrix(0, nrow = nrow(Z), ncol = 1)
  
  ## Estimate model
  fit <- mmsbm_fit(t(Z),
                   t(X),
                   Y,
                   t_id_d,
                   t_id_n,
                   nodes_pp,
                   nt_id,
                   mu_b,
                   var_b,
                   ctrl$phi_init_t,
                   ctrl$kappa_init_t,
                   ctrl$b_init_t,
                   c(ctrl$beta_init),
                   ctrl$gamma_init,
                   ctrl
  )
  
  
  ##Reorder to match original order
  fit[["MixedMembership"]] <- fit[["MixedMembership"]][,order(monadic_order)] 
  fit[["Kappa"]] <- fit[["Kappa"]][,order(time_order)]
  fit[["BlockModel"]] <- t(fit[["BlockModel"]])
  
  ## Rescale coefficients
  fit[["DyadCoef"]] <- fit[["DyadCoef"]] / Z_sd[-which(Z_sd==0)]
  fit[["MonadCoef"]] <- vapply(fit[["MonadCoef"]],
                               function(mat, sd_vec, mean_vec){
                                 constx <- which(sd_vec==0)
                                 mat[-constx, ] <- mat[-constx, ] / sd_vec[-constx]
                                 if(length(constx)!=0)
                                   mat[constx, ] <- mat[constx, ] - mean_vec[-constx] %*% mat[-constx, ]
                                 return(mat)
                               },
                               array(0.0, c(ncol(X), n.blocks)),
                               sd_vec = X_sd,
                               mean_vec = X_mean)
  ## Include used data in original order
  fit$monadic.data <- mfm[order(monadic_order),]
  fit$dyadic.data <- mfd[order(dyadic_order),]
  fit$Y <- Y[order(dyadic_order)]
  
  ## Include internal node id's
  fit$InternalNodeIndex <- nt_id[order(dyadic_order),]
  
  ## Include original call
  fit$call <- match.call()
  
  ##Assign class for methods
  class(fit) <- "mmsbm"
  
  return(fit)
}







