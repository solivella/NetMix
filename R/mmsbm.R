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
#'     in dyad. Cannot contain special charecter "`@`". 
#' @param receiverID Character string. Quoted name of the variable in \code{data.dyad} identifying 
#'     the receiver node. For undirected networks, the variable simply contains name of second node 
#'     in dyad. Cannot contain special charecter "`@`".
#' @param nodeID Character string. Quoted name of the variable in \code{data.monad} identifying 
#'     a node in either \code{data.dyad[,senderID]} or \code{data.dyad[,senderID]}. If not \code{NULL},
#'     every node \code{data.dyad[,senderID]} or \code{data.dyad[,senderID]} must be present in 
#'     \code{data.monad[,nodeID]}. Cannot contain special charecter "`@`".
#' @param timeID Character string. Quoted name of the variable in both \code{data.dyad} and
#'     \code{data.monad} indicating the time in which network (and correspding nodal atrributes)
#'     were observed. The variable itself must be composed of integers. Cannot contain special charecter "`@`".
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
#' @param mmsbm.control A named list of optional algorithm control parameters.
#'     \describe{
#'        \item{seed}{Integer value. Seed the RNG. By default, a random seed is generated and returned for reproducibility purposes.}
#'        \item{spectral}{Boolean. Type of initialization algorithm for mixed-membership vectors in static case. If \code{TRUE} (default),
#'                    use spectral clustering with degree correction; otherwise, use kmeans algorithm.}
#'        \item{init_gibbs}{Boolean. Should a collapsed Gibbs sampler of non-regression mmsbm be used to initialize
#'                    mixed-membership vectors, instead of a spectral or simple kmeans initialization?
#'                    Setting to \code{TRUE} will result in slower initialization and faster model estimation. When \code{TRUE}, results are typically very sensitive to
#'                    choice of alpha (see below).}            
#'        \item{alpha}{Numeric positive value. Concentration parameter for collapsed Gibbs sampler to find initial
#'                     mixed-membership values when \code{init_gibbs=TRUE}. Defaults to 1.0.}            
#'        \item{missing}{Means of handling missing data. One of "indicator method" (default) or "listwise deletion".}       
#'        \item{em_iter}{Number of maximum iterations in variational EM. Defaults to 5e3.}
#'        \item{opt_iter}{Number of maximum iterations of BFGS in M-step. Defaults to 10e3.}
#'        \item{hessian}{Boolean indicating whether the Hessian matrix of regression coefficients should e returned. Defaults to \code{TRUE}.}
#'        \item{mu_b}{Numeric vector with two elements: prior mean of blockmodel's main diagonal elements, and
#'                    and prior mean of blockmodel's offdiagonal elements. Defaults to \code{c(5.0, -5.0)}.}
#'        \item{var_b}{Numeric vector with two positive elements: prior variance of blockmodel's main diagonal elements, and
#'                    and prior variance of blockmodel's offdiagonal elements. Defaults to \code{c(1.0, 1.0)}.}
#'        \item{var_beta}{Numeric positive value. (Gaussian) Prior variance of monadic coefficients. Defaults to 5.0.}
#'        \item{var_gamma}{Numeric positive value. (Gaussian) Prior variance of dyadic coefficients. Defaults to 5.0.}
#'        \item{eta}{Numeric positive value. Concentration hyper-parameter for HMM. Defaults to 10.3.}
#'        \item{phi_init_t}{Matrix, \code{n.blocks} by total number of nodes across years. Optional initial values for variational
#'                       parameters for mixed-membership vectors. Column names must be of the form \code{nodeid\@year }.}
#'        \item{kappa_init_t}{Matrix, \code{n.hmmstates} by number of years. Optional initial values for variational 
#'                       parameters for state probabilities.}
#'        \item{b_init_t}{Matrix, \code{n.blocks} by \code{n.blocks}. Optional initial values for blockmodel.}
#'        \item{beta_init}{Array, predictors by \code{n.blocks} by \code{n.hmmstates}. Optional initial values for monadic coefficients.}
#'        \item{gamma_init}{Vector. Optional initial values for dyadic coefficients.}
#'        \item{permute}{Boolean. Should all permutations be tested to realign initial block models in dynamic case? If \code{FALSE}, realignment is 
#'                      done via faster graph matching algorithm, but may not be exact. Defaults to \code{TRUE}.}
#'        \item{conv_tol}{Numeric value. Absolute tolerance for VI convergence. Defaults to 1e-3.}
#'        \item{verbose}{Boolean. Should extra information be printed as model iterates? Defaults to FALSE.}
#'        }
#'       
#' @return Object of class \code{mmsbm}. List with named components:
#'     \describe{
#'       \item{MixedMembership}{Matrix of variational posterior of mean of mixed-membership vectors. \code{nodes} by \
#'                              \code{n.groups}.}
#'       \item{BlockModel}{\code{n.groups} by \code{n.groups} matrix of estimated tie log-odds between members
#'                         of corresponding latent groups. The blockmodel.}
#'       \item{vcov_blockmodel}{If \code{hessian=TRUE}, variance-covariance matrix of parameters in blockmodel, ordered in column-major order.}
#'       \item{MonadCoef}{Array of estimated coefficient values for monadic covariates. Has \code{n.groups} columns,
#'                        and \code{n.hmmstates} slices.}
#'       \item{vcov_monad}{If \code{hessian=TRUE}, variance-covariance matrix of monadic coefficients.}                  
#'       \item{DyadCoef}{Vector estimated coefficient values for dyadic covariates.}
#'       \item{vcov_dyad}{If \code{hessian=TRUE}, variance-covariance matrix of dyadic coefficients.}
#'       \item{TransitionKernel}{Matrix of estimated HMM transition probabilities.}
#'       \item{Kappa}{Matrix of marginal probabilities of being in an HMM state at any given point in time. 
#'                    \code{n.hmmstates} by years (or whatever time interval networks are observed at).}
#'       \item{niter}{Final number of VI iterations.}
#'       \item{converged}{Convergence indicator; zero indicates failure to converge.}
#'       \item{NodeIndex}{Order in which nodes are stored in all return objects.}
#'       \item{monadic.data, dyadic.data}{Model frames used during estimation (stripped of attributes).}
#'       \item{forms}{Values of selected formal arguments used by other methods.}
#'       \item{seed}{The value of RNG seed used during estimation.}
#'       \item{call}{Original (unevaluated) function call.}
#'     }
#' @author Santiago Olivella (olivella@@unc.edu), Adeline Lo (adelinel@@princeton.edu), Tyler Pratt (tyler.pratt@@yale.edu), Kosuke Imai (imai@@harvard.edu)
#' 
#' @examples 
#' library(NetMix)
#' ## Load datasets
#' data("lazega_dyadic")
#' data("lazega_monadic")
#' ## Estimate model with 2 groups
#' ## Setting to `hessian=TRUE` increases computation time
#' ## but is needed if standard errors are to be computed. 
#' lazega_mmsbm <- mmsbm(SocializeWith ~ Coworkers,
#'                       ~  School + Practice + Status,
#'                       senderID = "Lawyer1",
#'                       receiverID = "Lawyer2",
#'                       nodeID = "Lawyer",
#'                       data.dyad = lazega_dyadic,
#'                       data.monad = lazega_monadic,
#'                       n.blocks = 2,
#'                       mmsbm.control = list(seed = 123,
#'                                            hessian = FALSE))
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
                  n.hmmstates = 1,
                  directed = TRUE,
                  mmsbm.control = list()){
  
  cl <- match.call(expand.dots = FALSE)
  formulas <- cl[match(c("formula.dyad","formula.monad"), names(cl))]
  
  ## Form default control list
  ctrl <- list(blocks = n.blocks,
               states = n.hmmstates,
               times = 1,
               directed = directed,
               seed = sample(500, 1),
               phi_init_t = NULL,
               kappa_init_t = NULL,
               b_init_t = NULL,
               beta_init = NULL,
               gamma_init = NULL,
               spectral = TRUE,
               init_gibbs = ifelse(n.hmmstates > 1, TRUE, FALSE),
               alpha = 1.0,
               missing="indicator method",
               em_iter = 50,
               hessian = TRUE,
               se_sim = 10,
               opt_iter = 10e3,
               mu_b = c(1.0, 0.0),
               var_b = c(1.0, 1.0),
               var_beta = 1.0,
               var_gamma = 1.0,
               eta = 1.0,
               permute = TRUE,
               conv_tol = 1e-3,
               verbose = FALSE)
  ctrl[names(mmsbm.control)] <- mmsbm.control
  set.seed(ctrl$seed)
  
  mu_b <- var_b <- array(NA, c(n.blocks, n.blocks))
  diag(mu_b) <- ctrl[["mu_b"]][1]
  mu_b[upper.tri(mu_b)|lower.tri(mu_b)] <- ctrl[["mu_b"]][2]
  diag(var_b) <- ctrl[["var_b"]][1]
  var_b[upper.tri(var_b)|lower.tri(var_b)] <- ctrl[["var_b"]][2]
  
  if(ctrl$verbose){
    cat("Pre-processing data...\n")
  }
  
  stopifnot(identical(class(formula.dyad), "formula"),
            identical(class(formula.monad), "formula"),
            is.data.frame(data.dyad))
  if(!is.null(data.monad)){
    stopifnot(is.data.frame(data.monad),
              !is.null(nodeID))
  }
  if(((length(formula.monad)>2) || (formula.monad[2] != "1()")) & is.null(data.monad)){
    stop("Monadic dataset not defined.")
  }
  
  ## Add time variable if null or single period
  if(is.null(timeID) || (length(unique(data.dyad[[timeID]])) == 1)){
    timeID <- "tid"
    data.dyad[timeID] <- 1
    if(!is.null(data.monad)) {
      data.monad[timeID] <- 1  
    }
  }
  
  ## Address missing data 
  if(identical(ctrl$missing, "indicator method")){
    # dyadic dataset
    if(length(all.vars(formula.dyad[[3]]))){
      miss.d <- apply(as.matrix(data.dyad[,all.vars(formula.dyad[[3]]), drop = FALSE]), 2, function(x){length(na.omit(x))}) < nrow(data.dyad)
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
  if(identical(ctrl$missing, "listwise deletion")){
    if(length(all.vars(formula.dyad[[3]]))){
      mdyad <- apply(as.matrix(data.dyad[,all.vars(formula.dyad[[3]])]), 1, function(x){!any(is.na(x))})
    } else {
      mdyad <- TRUE
    }
    if(length(all.vars(formula.monad[[2]]))){
      mmonad <- apply(as.matrix(data.monad[,all.vars(formula.monad[[2]])]), 1, function(x){!any(is.na(x))})
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

  
  ut <- unique(mfd[["(tid)"]])
  periods <- length(ut)
  if(periods > 1){
    ctrl$times <- periods
    if((n.hmmstates > 1) & is.null(mmsbm.control$eta)){
      ctrl$eta <- periods/n.hmmstates 
    }
  }
  dntid <- cbind(do.call(paste, c(mfd[c("(sid)","(tid)")], sep = "@")),
                 do.call(paste, c(mfd[c("(rid)","(tid)")], sep = "@")))
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
  ntid <- do.call(paste, c(mfm[c("(nid)","(tid)")], sep="@"))
  if(!all(dntid %in% ntid))
    stop("Nodes in dyadic dataset missing from monadic dataset. Are node and time identifiers identical in data.dyad and data.monad?")
  match_ids <- ntid %in% dntid
  if(any(!match_ids)){
    if(ctrl$verbose){
      cat("\tSome nodes in data.monad not present in data.dyad; dropping them.\n")
    }
    mfm <- mfm[match_ids, ]
    ntid <- do.call(paste, c(mfm[c("(nid)","(tid)")], sep="@"))
  }
  
  Y <- stats::model.response(mfd)
  X <- base::scale(model.matrix(terms(mfm), mfm))
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
  
  nt_id <- cbind(match(dntid[,1], ntid) - 1, match(dntid[,2], ntid) - 1)
  t_id_d <- match(mfd[["(tid)"]], ut) - 1
  t_id_n <- match(mfm[["(tid)"]], ut) - 1
  nodes_pp <- c(by(mfm, mfm[["(tid)"]], nrow))
  dyads_pp <- c(by(mfd, mfd[["(tid)"]], nrow))
  
  
  
  ## Create initial values
  if(ctrl$verbose){
    cat("Obtaining initial values...\n")
  }
  
  dyads <- split.data.frame(dntid, mfd[, "(tid)"])
  edges <- split(Y, mfd[, "(tid)"])
  soc_mats <- Map(function(dyad_mat, edge_vec){
    nodes <- unique(c(dyad_mat))
    nnode <- length(nodes)
    adj_mat <- matrix(NA,
                      nnode,
                      nnode,
                      dimnames = list(nodes,
                                      nodes))
    adj_mat[dyad_mat] <- edge_vec
    if(!directed){
      adj_mat[dyad_mat[,c(2,1)]] <- edge_vec
    }
    obs_prop <- mean(adj_mat, na.rm = TRUE)
    if(anyNA(adj_mat)){
      if(is.nan(obs_prop)){
        obs_prop <- 0.01
      }
      adj_mat[is.na(adj_mat)] <- rbinom(sum(is.na(adj_mat)), 1, obs_prop)
    }
    diag(adj_mat) <- 0
    if(!directed){
      mat_ind <- which(upper.tri(adj_mat), arr.ind = TRUE)
      adj_mat[mat_ind[,c(2,1)]] <- adj_mat[upper.tri(adj_mat)]
    }
    return(adj_mat)
  }, dyads, edges)
  
  
  if(is.null(ctrl$kappa_init_t)){
    if((periods > 1) & (n.hmmstates > 1)){
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
      state_init <- fitted(kmeans(x = dyad_time,
                                  centers = n.hmmstates,
                                  iter.max = 15,
                                  nstart = 15), "classes")
      kappa_internal <- model.matrix(~ factor(state_init, 1:n.hmmstates) - 1)
      kappa_internal <- .transf(kappa_internal)
      ctrl$kappa_init_t <- t(kappa_internal)
    } else {
      ctrl$kappa_init_t <- t(matrix(1, nrow = periods))
      state_init <- apply(ctrl$kappa_init_t, 2, which.max)
    }
  } else {
    state_init <- apply(ctrl$kappa_init_t, 2, which.max)
  } 
  if(identical(n.hmmstates, 1)){
    names(state_init) = 1
  }
  
  
  if(is.null(ctrl$phi_init_t)){
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
                                          nstart = 10), "classes")
          
        } else {
          init_c <- sample(1:nrow(target), n.blocks, replace = FALSE)
          cents <- jitter(target[init_c, ])
          clust_internal <- fitted(suppressWarnings(kmeans(x = target,
                                                           centers = cents,
                                                           iter.max = 15,
                                                           algorithm = "Lloyd",
                                                           nstart = 1)), "classes")
        }
        
        phi_internal <- model.matrix(~ factor(clust_internal, 1:n.blocks) - 1)
        phi_internal <- .transf(phi_internal)
        rownames(phi_internal) <- rownames(soc_mats[[i]])
        colnames(phi_internal) <- 1:n.blocks
        MixedMembership <- t(phi_internal)
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
        temp_res[[i]] <- list(BlockModel = BlockModel,
                              MixedMembership = MixedMembership)
        
        
      }
    }
    block_models <- lapply(temp_res, function(x)x$BlockModel)
    target_ind <- which.max(sapply(soc_mats, ncol))
    perms_temp <- .findPerm(block_models, target_mat = block_models[[target_ind]], use_perms = ctrl$permute)
    phis_temp <- lapply(temp_res, function(x)x$MixedMembership) 
    phi.ord <- as.numeric(lapply(phis_temp, function(x)strsplit(colnames(x), "@")[[1]][2])) # to get correct temporal order
    ctrl$phi_init_t <- do.call(cbind,mapply(function(phi,perm){perm %*% phi},
                                            phis_temp[order(phi.ord)], perms_temp, SIMPLIFY = FALSE)) 
    rownames(ctrl$phi_init_t) <- 1:n.blocks
  } 
  
  if(is.null(ctrl$gamma_init)){
    if(ncol(Z) > 0){
      ctrl$gamma_init <- coef(lm.fit(as.matrix(Z),log((Y+1e-5)/(1-Y+1e-5))))
    } else {
      ctrl$gamma_init <- 0
    }
  }
  if(ncol(Z) == 0)
    Z <- matrix(0, nrow = nrow(Z), ncol = 1)
  if(anyNA(ctrl$gamma_init)){
    stop("Singular design matrix; check dyadic predictors.")
  }
  
  if(is.null(ctrl$b_init_t)){
    ctrl$b_init_t <- qlogis(approxB(Y, nt_id, ctrl$phi_init_t))
    if(any(is.infinite(ctrl$b_init_t))){
      which.inf <- which(is.infinite(ctrl$b_init_t))
      ctrl$b_init_t[which.inf] <- ifelse(ctrl$b_init_t[which.inf] > 0, 25, -25) 
    }
  }
  if(is.null(ctrl$beta_init)){
    X_state <- split.data.frame(X, state_init[t_id_n + 1])
    phi_state <- split.data.frame(t(ctrl$phi_init_t), state_init[t_id_n + 1])
    ctrl$beta_init <- mapply(function(X_sub, phi_sub){
      phi_temp <- .transf(phi_sub)
      lm.fit(X_sub, log(phi_temp))$coefficients
    },
    X_state, phi_state,
    SIMPLIFY = "array")
  }
  if(anyNA(ctrl$beta_init)){
    stop("Nearly singular design matrix; check monadic predictors.")
  }
  ## Estimate model
  if(ctrl$verbose){
    cat("Estimating model...\n");
  }
  X_t <- t(X)
  Z_t <- t(Z)
  fit <- mmsbm_fit(Z_t,
                    X_t,
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
                    ctrl$beta_init,
                    ctrl$gamma_init,
                    ctrl
  )
  if((!fit[["converged"]]) & ctrl$verbose)
    warning(paste("Model did not converge after", fit[["niter"]], "iterations.\n"))
  else if (ctrl$verbose){
    cat("done after", fit[["niter"]], "iterations.\n")
  }
  
  
  ##Return transposes 
  fit[["TransitionKernel"]] <- t(fit[["TransitionKernel"]])
  fit[["BlockModel"]] <- t(fit[["BlockModel"]])
  
  
  ## Rescale and name coefficients
  fit[["DyadCoef"]] <- fit[["DyadCoef"]] / Z_sd[-which(Z_sd==0)]
  if(length(fit[["DyadCoef"]])){
    Z <- t(t(Z) * Z_sd[-1] + Z_mean[-1])
    fit[["BlockModel"]] <- fit[["BlockModel"]] - c(Z_mean[-constz] %*% fit[["DyadCoef"]])
    names(fit[["DyadCoef"]]) <- colnames(Z) 
 }

  fit[["MonadCoef"]] <- vapply(1:n.hmmstates,
                               function(ind, coefs, sd_vec, mean_vec){
                                 mat <- coefs[,,ind, drop=FALSE]
                                 constx <- which(sd_vec==0)
                                 mat[-constx, , 1] <- mat[-constx, , 1] / sd_vec[-constx]
                                 if(length(constx)!=0){
                                   mat[constx, ,1] <- mat[constx, ,1] - mean_vec[-constx] %*% mat[-constx, , 1]
                                 }
                                 return(mat)
                               },
                               array(0.0, c(ncol(X), n.blocks)),
                               coefs = fit[["MonadCoef"]],
                               sd_vec = X_sd,
                               mean_vec = X_mean)
  rownames(fit[["MonadCoef"]]) <- colnames(X)
  colnames(fit[["MonadCoef"]]) <- paste("Group",1:n.blocks)
  X <- t(t(X) * X_sd + X_mean) #unscale
  
  ## Add other names
  colnames(fit[["Kappa"]]) <- unique(mfm[,"(tid)"])
  dimnames(fit[["BlockModel"]]) <- replicate(2,paste("Group",1:n.blocks), simplify = FALSE)
  dimnames(fit[["TransitionKernel"]]) <- replicate(2,paste("State",1:n.hmmstates), simplify = FALSE)
  colnames(fit[["MixedMembership"]]) <- ntid
  
  
  if(ctrl$hessian){
    if(ctrl$verbose){
      cat("Computing approximate vcov. matrices...\n")
    }
    
    ## Compute approximate standard errors
    ## for monadic coefficients
    all_phi <- split.data.frame(rbind(t(fit[["SenderPhi"]]),
                                      t(fit[["ReceiverPhi"]])),
                                c(nt_id))
    sampleC_perm <- lapply(all_phi,
                           function(mat){
                             apply(mat, 2, function(vec)poisbinom::rpoisbinom(ctrl$se_sim, vec))
                           })
    sampleC_perm <- cbind(do.call(rbind, sampleC_perm), # samples
                          rep(1:length(all_phi), each = ctrl$se_sim), #node id
                          rep(1:ctrl$se_sim, times = length(all_phi))) #sample id
    sampleC_perm <- sampleC_perm[order(sampleC_perm[,n.blocks + 2], sampleC_perm[,n.blocks + 1]),]
    C_samples <- split.data.frame(sampleC_perm[,1:n.blocks], sampleC_perm[,n.blocks + 2])
    S_samples <- replicate(ctrl$se_sim, apply(fit[["Kappa"]], 2, function(x)sample(1:n.hmmstates, 1, prob = x)), simplify = FALSE)
    hessBeta_list <- mapply(
      function(C_samp, S_samp, tidn, X_i, Nvec, beta_vec, vbeta, periods)
      {
        if(n.hmmstates > 1) {
          s_matrix <- t(model.matrix(~factor(S_samp, 1:n.hmmstates) - 1))
        } else {
          s_matrix <- matrix(1, ncol=periods)
        }
        tot_in_state <- rowSums(s_matrix)
        if(any(tot_in_state == 0.0)){
          stop("Some HMM states are empty; consider reducing n.hmmstates, or increasing eta.")
        }  
        return(optimHess(c(beta_vec),alphaLB,
                         tot_nodes = Nvec,
                         c_t = t(C_samp),
                         x_t = X_i,
                         s_mat = s_matrix,
                         t_id = tidn,
                         var_beta = vbeta))
      },
      C_samples, S_samples,
      MoreArgs = list(tidn = t_id_n,
                      X_i = t(X),
                      Nvec = fit[["TotNodes"]],
                      beta_vec = fit[["MonadCoef"]], 
                      vbeta = ctrl$var_beta, 
                      periods = periods),
      SIMPLIFY=FALSE)
    hessBeta <- Reduce("+", hessBeta_list)/ctrl$se_sim
    
    fit$vcov_monad <- solve(hessBeta)
    colnames(fit$vcov_monad) <- rownames(fit$vcov_monad) <- paste(rep(paste("State",1:n.hmmstates), each = prod(dim(fit[["MonadCoef"]])[1:2])),
                                                                  rep(colnames(fit[["MonadCoef"]]), each = nrow(fit[["MonadCoef"]]), times = n.hmmstates),
                                                                  rep(rownames(fit[["MonadCoef"]]), times = n.blocks*n.hmmstates),
                                                                  sep=":") 
    
    ## and for dyadic coefficients
    z_samples <- replicate(ctrl$se_sim, getZ(fit[["SenderPhi"]]), simplify = FALSE)
    w_samples <- replicate(ctrl$se_sim, getZ(fit[["ReceiverPhi"]]), simplify = FALSE)
    all_theta_par <- c(
      if(directed){
        c(fit[["BlockModel"]]) }
      else{
        fit[["BlockModel"]][lower.tri(fit[["BlockModel"]], diag = TRUE)]
      }, 
      fit[["DyadCoef"]])
    hessTheta_list <- mapply(
      function(send_samp, rec_samp, y_vec, Z_d, par_theta, mu_b_mat, var_b_mat, var_g, dir_net)
      {
        optimHess(par_theta,
                  thetaLB, 
                  y = y_vec,
                  z_t = Z_d,
                  send_phi = send_samp,
                  rec_phi = rec_samp,
                  mu_b_t = mu_b_mat,
                  var_b_t = var_b_mat,
                  var_gamma = var_g,
                  directed = dir_net)
      },
      z_samples, w_samples,
      MoreArgs = list(par_theta = all_theta_par, 
                      y_vec = Y,
                      Z_d = t(Z),
                      mu_b_mat = mu_b,
                      var_b_mat = var_b,
                      var_g = ctrl$var_gamma, 
                      dir_net = directed),
      SIMPLIFY = FALSE)
    
    hessTheta <- Reduce("+", hessTheta_list)/ctrl$se_sim
    vcovTheta <- solve(hessTheta)
    
    
    N_B_PAR <- ifelse(directed, n.blocks^2 , n.blocks * (1 + n.blocks) / 2)
    fit$vcov_blockmodel <- vcovTheta[1:N_B_PAR, 1:N_B_PAR, drop = FALSE]
    bm_names <- outer(rownames(fit[["BlockModel"]]), colnames(fit[["BlockModel"]]), paste, sep=":")
    colnames(fit$vcov_blockmodel) <- rownames(fit$vcov_blockmodel) <- if(directed){c(bm_names)}else{c(bm_names[lower.tri(bm_names, TRUE)])}
    
    
    if(any(Z_sd > 0)){
      fit$vcov_dyad <- vcovTheta[(N_B_PAR + 1):nrow(vcovTheta),
                                 (N_B_PAR + 1):ncol(vcovTheta),
                                 drop = FALSE]
      colnames(fit$vcov_dyad) <- rownames(fit$vcov_dyad) <- names(fit[["DyadCoef"]])
    }
    
    
    if(ctrl$verbose){
      cat("done.\n")
    }
    
  }
  
  
  ## Include used data 
  attr(mfm, "terms") <- NULL
  attr(mfd, "terms") <- NULL
  fit$monadic.data <- mfm
  fit$dyadic.data <- mfd
  fit$Y <- Y
  
  ## Include node id's
  fit$NodeIndex <- nt_id
  
  ## Include a few formals needed by other methods
  fit$forms <- list(directed = directed,
                    senderID = senderID,
                    receiverID = receiverID,
                    timeID = timeID,
                    nodeID = nodeID,
                    t_id_d = t_id_d,
                    hessian = ctrl$hessian,
                    formula.dyad = formulas[[1]],
                    formula.monad = formulas[[2]])
  
  ## Include used seed
  fit$seed <- ctrl$seed
  
  fit$call <- cl
  
  ##Assign class for methods
  class(fit) <- "mmsbm"
  
  return(fit)
}
