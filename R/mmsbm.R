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
#' @param n.blocks Integer value. How many latent groups should be used to estimate the model?
#' @param n.hmmstates Integer value. How many hidden Markov state should be used in the HMM? Defaults 
#'    to 1 (i.e. no HMM).  
#' @param directed Boolean. Is the network directed? Defaults to \code{TRUE}.
#' @param missing Means of handling missing data. One of "indicator method" (default) or "listwise deletion".
#' @param mmsbm.control A named list of optional algorithm control parameters.
#'     \describe{
#'        \item{spectral}{Boolean. Type of initialization algorithm for mixed-membership vectors. If \code{TRUE},
#'                    use spectral clustering with degree correction; otherwise, use kmeans algorithm.}
#'        \item{seed}{RNG seed. Defaults to \code{NULL}, which does not seed the RNG}            
#'        \item{em_iter}{Number of maximum iterations in variational EM. Defaults to 5e3}
#'        \item{opt_iter}{Number of maximum iterations of BFGS in M-step. Defaults to 10e3}
#'        \item{mu_b}{Numeric vector with two elements: prior mean of blockmodel's main diagonal elements, and
#'                    and prior mean of blockmodel's offdiagonal elements. Defaults to \code{c(5.0, -5.0)}}
#'        \item{var_b}{Numeric vector with two positive elements: prior variance of blockmodel's main diagonal elements, and
#'                    and prior variance of blockmodel's offdiagonal elements. Defaults to \code{c(1.0, 1.0)}}
#'        \item{var_beta}{Numeric positive value. (Gaussian) Prior variance of monadic coefficients. Defaults to 5.0.}
#'        \item{var_gamma}{Numeric positive value. (Gaussian) Prior variance of dyadic coefficients. Defaults to 5.0.}
#'        \item{eta}{Numeric positive value. Concentration hyper-parameter for HMM. Defaults to 10.3}
#'        \item{phi_init_t}{Matrix, \code{n.blocks} by total number of nodes across years. Optional initial values for variational
#'                       parameters for mixed-membership vectors. Column names must be of the form \code{nodeid\@year}
#'                       }
#'        \item{kappa_init_t}{Matrix, \code{n.hmmstates} by number of years. Optional initial values for variational 
#'                       parameters for state probabilities.}
#'        \item{b_init_t}{Matrix, \code{n.blocks} by \code{n.blocks}. Optional initial values for blockmodel.}
#'        \item{beta_init}{Array, predictors by \code{n.blocks} by \code{n.hmmstates}. Optional initial values for monadic coefficients.}
#'        \item{gamma_init}{Vector. Optional initial values for dyadic coefficients.}
#'        \item{permute}{Boolean. Should all permutations be tested to realign initial block models in dynamic case? If \code{FALSE}, realignment is 
#'                      done via faster graph matching algorithm, but may not be exact. Defaults to \code{TRUE}.}
#'        \item{threads}{Numeric integer. Number of available cores for paralellization. Defaults to 4}
#'        \item{conv_tol}{Numeric value. Absolute tolerance for VI convergence. Defaults to 1e-4}
#'        \item{verbose}{Boolean. Should extra information be printed as model iterates? Defaults to FALSE}
#'        }
#'       
#' @return Object of class \code{mmsbm}. List with named components:
#'     \describe{
#'       \item{MixedMembership}{Matrix of variational posterior of mean of mixed-membership vectors. \code{nodes} by \
#'                              \code{n.groups}}
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
#'       \item{NodeIndex}{Order in which nodes are stored in all return objects.}
#'       \item{monadic.data, dyadic.data, n_states, n_blocks, directed}{Original values of parameters used during estimation}
#'     }
#' @author Kosuke Imai (imai@@harvard.edu), Tyler Pratt (tyler.pratt@@yale.edu), Santiago Olivella (olivella@@unc.edu)
#' 
#' @example tests/Examples/MIDColdWar.R
#' 

mmsbm <- function(formula.dyad,
                  formula.monad=~1,
                  senderID, 
                  receiverID,
                  nodeID = NULL,
                  timeID = NULL,
                  data.dyad,
                  data.monad = NULL,
                  n.blocks,
                  n.hmmstates,
                  directed = TRUE,
                  missing="indicator method",
                  mmsbm.control = list()){
  stopifnot(class(formula.dyad) == "formula",
            class(formula.monad) == "formula",
            is.data.frame(data.dyad))
  if(!is.null(data.monad)){
    stopifnot(is.data.frame(data.monad),
              !is.null(nodeID))
  }
  if(((length(formula.monad)>2) || (formula.monad[2] != "1()")) & is.null(data.monad)){
    stop("Monadic dataset not defined.")
  }
  
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
    # monadic dataset
    if(length(all.vars(formula.monad[[2]]))){
      miss.m <- apply(data.monad[,all.vars(formula.monad[[2]]), drop = FALSE], 2, function(x){length(na.omit(x))}) < nrow(data.monad)
      mm <- names(miss.m[miss.m])
      if(length(mm)>0){
        m.ind <- apply(as.data.frame(data.monad[,mm]), 2, function(x){
          ifelse(is.na(x), 1, 0)
        })
        colnames(m.ind) <- paste(mm, "_missing", sep="")
        data.monad[,mm] <- as.vector(apply(as.data.frame(data.monad[,mm]), 2, function(x){
          x[is.na(x)] <- 0
          return(x)
        }))
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
    if(length(all.vars(formula.monad[[2]]))){
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
  
  mfd <- do.call(model.frame, list(formula = formula.dyad,
                                   data = data.dyad,
                                   drop.unused.levels = TRUE,
                                   tid = as.name(timeID),
                                   sid = as.name(senderID),
                                   rid = as.name(receiverID)))
  #dyadic_order <- with(mfd, order(`(tid)`, `(sid)`, `(rid)`))
  #mfd <- mfd[dyadic_order, ]
  
  
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
  #monadic_order <- with(mfm, order(`(tid)`, `(nid)`)) 
  #mfm <- mfm[monadic_order, ]
  ntid <- do.call(paste, c(mfm[c("(nid)","(tid)")], sep="@"))
  mfm <- mfm[ntid %in% unique(c(dntid)), ]
  data.monad <- data.monad[ntid %in% unique(c(dntid)), ]
  #monadic_order <- with(mfm, order(`(tid)`, `(nid)`)) 
  
  if(!all(udntid %in% ntid))
    stop("Nodes in dyadic dataset missing from monadic dataset. Are node and time identifiers identical in data.dyad and data.monad?")
  
  
  Y <- model.response(mfd)
  X <- scale(model.matrix(terms(mfm), mfm))
  X_mean <- attr(X, "scaled:center")
  X_sd <- attr(X, "scaled:scale")
  if(any(X_sd==0)){
    constx <- which(X_sd==0)
    X[,constx] <- 1
  }
  n_monad_pred <- ncol(X)
  Z <- scale(model.matrix(terms(mfd), mfd))
  Z_mean <- attr(Z, "scaled:center")
  Z_sd <- attr(Z, "scaled:scale")
  if(any(Z_sd==0)){
    constz <- which(Z_sd==0)
    Z <- Z[,-constz, drop = FALSE]
  }
  n_dyad_pred <- ncol(Z)
  
  #nt_id <- cbind(.mapID(udntid, dntid[, 1]) - 1, .mapID(udntid, dntid[, 2]) - 1) # this was creating problems
  nt_id <- cbind(.mapID(ntid, dntid[, 1]) - 1, .mapID(ntid, dntid[, 2]) - 1) # this works
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
               seed = NULL,
               em_iter = 50,
               opt_iter = 10e3,
               mu_b = c(1.0, 0.0),
               var_b = c(1.0, 1.0),
               var_beta = 1.0,
               var_gamma = 1.0,
               eta = 1.3,
               permute = TRUE,
               threads = 1,
               conv_tol = 1e-2,
               verbose = FALSE)
  ctrl[names(mmsbm.control)] <- mmsbm.control
  if(!is.null(ctrl$seed)) set.seed(ctrl$seed)
  mu_b <- var_b <- array(NA, c(n.blocks, n.blocks))
  diag(mu_b) <- ctrl[["mu_b"]][1]
  mu_b[upper.tri(mu_b)|lower.tri(mu_b)] <- ctrl[["mu_b"]][2]
  diag(var_b) <- ctrl[["var_b"]][1]
  var_b[upper.tri(var_b)|lower.tri(var_b)] <- ctrl[["var_b"]][2]
  
  
  
  ## Create initial values
  if(ctrl$verbose)
    cat("Obtaining initial values...\n")
  all.nodes <- unique(unlist(mfd[,c("(sid)","(rid)")]))
  node.cols <- which(names(mfd)%in%c("(sid)","(rid)", "(tid)"))
  dyads <- split.data.frame(mfd[,c(node.cols, 1)], mfd[, "(tid)"])
  edges <- split(Y, mfd[, "(tid)"])
  
  soc_mats <- lapply(dyads,
                     function(dyad_df, 
                              nnode = length(all.nodes),
                              nodes = all.nodes)
                     {
                       indeces <- as.matrix(dyad_df[,c("(sid)","(rid)")])
                       time_ind <- unique(dyad_df[,"(tid)"]) 
                       adj_mat <-  matrix(NA, 
                                          nnode,
                                          nnode,
                                          dimnames = list(nodes,
                                                          nodes))
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
                       node_names <- paste(nodes, "@", time_ind, sep="")
                       dimnames(adj_mat) <- list(node_names,
                                                 node_names)
                       return(adj_mat)
                     })
  
  td_id <- cbind(mfd[,"(tid)"],paste(mfd[,"(sid)"],mfd[,"(rid)"], sep = "->"))
  dyad_time <- matrix(NA, periods, length(unique(td_id[,2])),
                      dimnames = list(ut,
                                      unique(td_id[,2])))
  dyad_time[td_id] <- Y
  if(any(is.na(dyad_time))){
    dyad_time <- apply(dyad_time, 2, function(x){
      x[is.na(x)] <- rbinom(sum(is.na(x)), 1, mean(x, na.rm=TRUE))
      return(x)
    })
  }
  
  if(is.null(ctrl$kappa_init_t)){
    if((periods > 1) & (n.hmmstates > 1)){
      state_internal <- fitted(kmeans(dyad_time,
                                      n.hmmstates,
                                      nstart = 15), "classes")
      kappa_internal <- model.matrix(~ as.factor(state_internal) - 1)
      kappa_internal <- .transf(kappa_internal)
      ctrl$kappa_init_t <- t(kappa_internal)
    } else {
      ctrl$kappa_init_t <- t(matrix(1, nrow = periods))
    }
  } 
  state_init <- apply(ctrl$kappa_init_t, 2, which.max)
  if(n.hmmstates==1){
    names(state_init) = 1
  }
  
  ## Adjusted this for when edges are never formed for a given node/year
  if(is.null(ctrl$phi_init_t)) {
    phi_init_temp <- lapply(soc_mats,
                            function(mat){
                              if(directed){
                                D_o <- diag(1/sqrt(rowSums(mat)))
                                D_o[D_o==Inf] <- 1 # check
                                D_i <- diag(1/sqrt(colSums(mat)))
                                D_i[D_i==Inf] <- 1 # check
                                U <- D_o %*% mat %*% D_i %*% t(mat) %*% D_o +
                                  D_i %*% t(mat) %*% D_o %*% mat %*% D_i
                              } else {
                                degs <- rowSums(mat)
                                max_d <- max(degs)
                                D <- 1 / sqrt(degs)
                                U <- t(D * mat) * D
                              }
                              if(ctrl$spectral) {
                                n_elem <- n.blocks + 1
                                res <- RSpectra::eigs_sym(U, n_elem)
                                # div <- res$vectors[,1] # check
                                # div[div==0] <- -1e-6 # check
                                # target <- res$vectors[,2:n_elem]/div # check
                                eta <- res$vectors[,1:n_elem] %*% diag(res$values[1:n_elem])
                                target <- eta[,2:n_elem] / eta[,1]
                                sig <- 1 - res$values[n_elem] / res$values[n.blocks]
                                if ( sig > 0.1 ){
                                  target <- target[,1:(n_elem - 2)]
                                }
                              } else {
                                target <- U
                              }
                              clust_internal <- fitted(kmeans(target,
                                                              n.blocks,
                                                              nstart = 50), "classes")
                              phi_internal <- model.matrix(~ as.factor(clust_internal) - 1)
                              phi_internal <- .transf(phi_internal)
                              rownames(phi_internal) <- rownames(mat)
                              colnames(phi_internal) <- 1:n.blocks
                              t(phi_internal)
                            })
    dyad_internal <- mapply(function(mat1, mat2)
    {
      cbind(match(paste(mat1[,2],"@",mat1[,1], sep=""), colnames(mat2)) - 1,
            match(paste(mat1[,3],"@",mat1[,1], sep=""), colnames(mat2)) - 1) 
    }, dyads, phi_init_temp, SIMPLIFY = FALSE)
    
    ctrl$phi_init_t <- do.call(cbind, phi_init_temp)
    ctrl$phi_list <- phi_init_temp
  } #else {
  #ctrl$phi_init_t <- ctrl$phi_init_t[, monadic_order]
  #}
  
  if(is.null(ctrl$gamma_init)){
    ctrl$gamma_init <- if(ncol(Z) > 0){
      coef(glm(Y ~ Z - 1, family = binomial()))
    } else {
      0
    }
  }
  if(any(is.na(ctrl$gamma_init))){
    stop("Nearly singular design matrix; check dyadic predictors.")
  }
  
  if(ncol(Z) == 0)
    Z <- matrix(0, nrow = nrow(Z), ncol = 1)
  
  ## Estimate models for multi-year inits
  if(periods > 1 & !is.null(ctrl$phi_list)){
    mfd_list <- split(mfd, mfd[,c("(tid)")])
    mfm_list <- split(mfm, mfm[,c("(tid)")])
    temp_res <- lapply(1:periods, function(x){
      ret <- mmsbm(update(formula.dyad, .~1),
                   formula.monad = ~ 1,
                   senderID = "(sid)",
                   receiverID = "(rid)",
                   nodeID = "(nid)",
                   timeID = "(tid)",
                   data.dyad = mfd_list[[x]],
                   data.monad = mfm_list[[x]][,c("(nid)", "(tid)")],
                   n.blocks = n.blocks,
                   n.hmmstates = 1,
                   directed = directed,
                   missing = missing,
                   mmsbm.control = list(em_iter = 5,
                                        mu_b = ctrl$mu_b,
                                        var_b = ctrl$var_b,
                                        phi_init_t = ctrl$phi_list[[x]],
                                        verbose = FALSE))
      #print(Cstack_info())
      # print(gc())
      return(ret)
    })
    block_models <- Map(function(x)x$BlockModel, temp_res)
    perms_temp <- .findPerm(block_models, perms = ctrl$permute)
    phis_temp <- Map(function(x)x$MixedMembership, temp_res)
    ctrl$phi_init_t <- do.call(cbind,mapply(function(phi,perm){perm %*% phi},
                                            phis_temp, perms_temp, SIMPLIFY = FALSE))
  }
  
  if(is.null(ctrl$b_init_t)){
    ctrl$b_init_t <- qlogis(approxB(Y, nt_id, ctrl$phi_init_t))
    if(any(is.infinite(ctrl$b_init_t))){
      which.inf <- which(is.infinite(ctrl$b_init_t))
      ctrl$b_init_t[which.inf] <- -25
    }
  }
  if(is.null(ctrl$beta_init)){
    ctrl$beta_init <- sapply(1:n.hmmstates,
                             function(m, dm, df, phi_i, states){
                               obsinm <- with(df, `(tid)` %in% unique(`(tid)`)[states == m])
                               phi_inm <- c(sapply(length(all.nodes) * (which(states==m)-1) + 1, 
                                                   function(x){seq(x, length.out=length(all.nodes))}))
                               
                               phi_temp <- .transf(t(phi_i[, phi_inm]))
                               phi_temp <- phi_temp[rownames(phi_temp) %in% paste(df[,"(nid)"], df[,"(tid)"], sep="@"),]
                               X_sub <- dm[obsinm, , drop = FALSE]
                               lm.fit(X_sub, log(phi_temp))$coefficients
                             },
                             dm = X, df = mfm, phi_i = ctrl$phi_init_t, states = state_init,
                             simplify = "array")
  }
  if(any(is.na(ctrl$beta_init))){
    stop("Nearly singular design matrix; check monadic predictors.")
  }
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
  
  ## Add names
  colnames(fit[["Kappa"]]) <- unique(mfm[,"(tid)"])
  dimnames(fit[["BlockModel"]]) <- replicate(2,paste("Group",1:n.blocks), simplify = FALSE)
  dimnames(fit[["TransitionKernel"]]) <- replicate(2,paste("State",1:n.hmmstates), simplify = FALSE)
  
  ##Reorder mixmem to match original order
  colnames(fit[["MixedMembership"]]) <- ntid
  #fit[["MixedMembership"]] <- fit[["MixedMembership"]][,order(monadic_order)] 
  
  ## Rescale and name coefficients
  fit[["DyadCoef"]] <- fit[["DyadCoef"]] / Z_sd[-which(Z_sd==0)]
  if(length(fit[["DyadCoef"]])){
    fit[["BlockModel"]] <- fit[["BlockModel"]] - c(Z_mean[-constz] %*% fit[["DyadCoef"]])
  }
  
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
                               array(0.0, c(ncol(X), n.blocks)),
                               sd_vec = X_sd,
                               mean_vec = X_mean)
  rownames(fit[["MonadCoef"]]) <- colnames(X)
  colnames(fit[["MonadCoef"]]) <- paste("Group",1:n.blocks)
  
  ## Include used data in original order
  #fit$monadic.data <- mfm[order(monadic_order),]
  #fit$dyadic.data <- mfd[order(dyadic_order),]
  #fit$Y <- Y[order(dyadic_order)]
  fit$monadic.data <- mfm
  fit$dyadic.data <- mfd
  fit$Y <- Y
  
  ## Include node id's
  fit$NodeIndex <- nt_id
  
  ## Include original call
  fit$call <- match.call()
  
  ##Assign class for methods
  class(fit) <- "mmsbm"
  
  return(fit)
}
