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
#'        \item{seed}{Integer. Seed the RNG. By default, a random seed is generated and returned for reproducibility purposes.}
#'        \item{nstart}{Integer. Number of random initialization trials. Defaults to 5.}
#'        \item{spectral}{Boolean. Type of initialization algorithm for mixed-membership vectors in static case. If \code{TRUE} (default),
#'                    use spectral clustering with degree correction; otherwise, use kmeans algorithm.}
#'        \item{init_gibbs}{Boolean. Should a collapsed Gibbs sampler of non-regression mmsbm be used to initialize
#'                    mixed-membership vectors, instead of a spectral or simple kmeans initialization?
#'                    Setting to \code{TRUE} will result in slower initialization and faster model estimation. When \code{TRUE}, results are typically very sensitive to
#'                    choice of alpha (see below).}            
#'        \item{alpha}{Numeric positive value. Concentration parameter for collapsed Gibbs sampler to find initial
#'                     mixed-membership values when \code{init_gibbs=TRUE}. Defaults to 1.0.}            
#'        \item{missing}{Means of handling missing data. One of "indicator method" (default) or "listwise deletion".}  
#'        \item{svi}{Boolean; should stochastic variational inference be used? Defaults to \code{TRUE}.}     
#'        \item{vi_iter}{Number of maximum iterations in stochastic variational updates. Defaults to 5e2.}
#'        \item{batch_size}{When \code{svi=TRUE}, proportion of nodes sampled in each local. Defaults to 0.05 when \code{svi=TRUE}, and to 1.0 otherwise.}                                 
#'        \item{forget_rate}{When \code{svi=TRUE}, value between (0.5,1], controlling speed of decay of weight of prior
#'                            parameter values in global steps. Defaults to 0.75 when \code{svi=TRUE}, and to 0.0 otherwise.}
#'        \item{delay}{When \code{svi=TRUE}, non-negative value controlling weight of past iterations in global steps. Defaults to 1.0 when \code{svi=TRUE},
#'                     and ignored otherwise.}                    
#'        \item{opt_iter}{Number of maximum iterations of BFGS in global step. Defaults to 10e3.}
#'        \item{hessian}{Boolean indicating whether the Hessian matrix of regression coefficients should e returned. Defaults to \code{TRUE}.}
#'        \item{assortative}{Boolean indicating whether blockmodel should be assortative (i.e. stronger connections within groups) or disassortative
#'                           (i.e. stronger connections between groups). Defaults to \code{TRUE}.}
#'        \item{mu_block}{Numeric vector with two elements: prior mean of blockmodel's main diagonal elements, and
#'                    and prior mean of blockmodel's offdiagonal elements. Defaults to \code{c(5.0, -5.0)} if \code{assortative=TRUE} (default)
#'                    and to \code{c(-5.0, 5.0)} otherwise.}
#'        \item{var_block}{Numeric vector with two positive elements: prior variance of blockmodel's main diagonal elements, and
#'                    and prior variance of blockmodel's offdiagonal elements. Defaults to \code{c(5.0, 5.0)}.}
#'        \item{mu_beta}{Either single numeric value, in which case the same prior mean is applied to all monadic coefficients, or
#'                       an array that is \code{npredictors} by \code{n.blocks} by \code{n.hmmstates}, where \code{npredictors}
#'                       is the number of monadic predictors for which a prior mean is being set (prior means need not be set for all)
#'                       predictors). The rows in the array should be named to identify which variables a prior mean is being set for.
#'                       Defaults to a common prior mean of 0.0 for all monadic coefficients.}            
#'        \item{var_beta}{See \code{mu_beta}. Defaults to a single common prior variance of 5.0 for all (standardized) monadic coefficients.}
#'        \item{mu_gamma}{Either a single numeric value, in which case the same prior mean is applied to all dyadic coefficients, or
#'                        a named vector of numeric values (with names corresponding to the name of the variable 
#'                       for which a prior mean is being set). Defaults to a common prior mean of 0.0 for all dyadic coefficients.}
#'        \item{var_gamma}{See \code{mu_gamma}. Defaults to a single common prior variance of 5.0 for all (standardized) dyadic coefficients.}
#'        \item{eta}{Numeric positive value. Concentration hyper-parameter for HMM. Defaults to 1.0.}
#'        \item{se_sim}{Number of samples from variational posterior of latent variables on which approximation to variance-covariance
#'                      matrices are based. Defaults to 10.}
#'        \item{dyad_vcov_samp}{Maximum number of dyads to sample in computation of variance-covariance of dyadic and blockmodel parameters, when compared to 
#'                              ten percent of the observed dyads. Defaults to 1000.}
#'        \item{fixed_mm}{Optional character vector, with \code{"nodeID@timeID"} as elements, indicating which mixed-membership vectors
#'                        should remain constant at their initial values throughout estimation. When only one year is observed, elements should be 
#'                         \code{"nodeID@1"}. Typically used with \code{mm_init_t}.}                      
#'        \item{mm_init_t}{Matrix, \code{n.blocks} by nodes across years. Optional initial values for mixed-membership vectors.
#'                           Although initial values need not be provided for all nodes, column names must have a \code{nodeID@timeID} format to 
#'                           avoid ambiguity. When only one year is observed, names should be \code{"nodeID@1"}.}
#'        \item{kappa_init_t}{Matrix, \code{n.hmmstates} by number of years. Optional initial values for variational 
#'                       parameters for state probabilities. Columns must be named according to unique year values.}
#'        \item{b_init_t}{Matrix, \code{n.blocks} by \code{n.blocks}. Optional initial values for blockmodel.}
#'        \item{beta_init}{Array, \code{predictors} by \code{n.blocks} by \code{n.hmmstates}. Optional initial values for monadic coefficients. If }
#'        \item{gamma_init}{Vector. Optional initial values for dyadic coefficients.}
#'        \item{permute}{Boolean. Should all permutations be tested to realign initial block models in dynamic case? If \code{FALSE}, realignment is 
#'                      done via faster graph matching algorithm, but may not be exact. Defaults to \code{TRUE}.}
#'        \item{conv_tol}{Numeric value. Absolute tolerance for VI convergence. Defaults to 1e-3.}
#'        \item{verbose}{Boolean. Should extra information be printed as model iterates? Defaults to FALSE.}
#'        }
#'       
#' @return Object of class \code{mmsbm}. List with named components:
#'     \describe{
#'       \item{MixedMembership}{Matrix of variational posterior of mean of mixed-membership vectors. \code{nodes} by
#'                              \code{n.blocks}.}
#'       \item{BlockModel}{\code{n.blocks} by \code{n.blocks} matrix of estimated tie log-odds between members
#'                         of corresponding latent groups. The blockmodel.}
#'       \item{vcov_blockmodel}{If \code{hessian=TRUE}, variance-covariance matrix of parameters in blockmodel, ordered in column-major order.}
#'       \item{MonadCoef}{Array of estimated coefficient values for monadic covariates. Has \code{n.blocks} columns,
#'                        and \code{n.hmmstates} slices.}
#'       \item{vcov_monad}{If \code{hessian=TRUE}, variance-covariance matrix of monadic coefficients.}                  
#'       \item{DyadCoef}{Vector estimated coefficient values for dyadic covariates.}
#'       \item{vcov_dyad}{If \code{hessian=TRUE}, variance-covariance matrix of dyadic coefficients.}
#'       \item{TransitionKernel}{Matrix of estimated HMM transition probabilities.}
#'       \item{Kappa}{Matrix of marginal probabilities of being in an HMM state at any given point in time. 
#'                    \code{n.hmmstates} by years (or whatever time interval networks are observed at).}
#'       \item{LowerBound}{Final LB value}
#'       \item{lb}{Vector of all LB across iterations, useful to check early convergence issues.}              
#'       \item{niter}{Final number of VI iterations.}
#'       \item{converged}{Convergence indicator; zero indicates failure to converge.}
#'       \item{NodeIndex}{Order in which nodes are stored in all return objects.}
#'       \item{monadic.data, dyadic.data}{Model frames used during estimation (stripped of attributes).}
#'       \item{forms}{Values of selected formal arguments used by other methods.}
#'       \item{seed}{The value of RNG seed used during estimation.}
#'       \item{call}{Original (unevaluated) function call.}
#'     }
#' @author Santiago Olivella (olivella@@unc.edu), Adeline Lo (aylo@@wisc.edu), Tyler Pratt (tyler.pratt@@yale.edu), Kosuke Imai (imai@@harvard.edu)
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
#'                                            conv_tol = 1e-2,
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
               svi = TRUE,
               nstarts = 5,
               spectral = TRUE,
               init_gibbs = ifelse(n.hmmstates > 1, TRUE, FALSE),
               threads = 1,
               alpha = 1.0,
               forget_rate = 0.75,
               delay = 1.0,
               batch_size = 0.05,
               missing="indicator method",
               vi_iter = 500,
               hessian = TRUE,
               se_sim = 10,
               dyad_vcov_samp = 100,
               opt_iter = 10e3,
               assortative = TRUE,
               mu_block = c(5.0, -5.0),
               var_block = c(5.0, 5.0),
               mu_beta = 0.0,
               mu_gamma = 0.0,
               var_beta = 5.0,
               var_gamma = 5.0,
               eta = 1.0,
               permute = TRUE,
               conv_tol = 1e-3,
               verbose = FALSE)
  ctrl[names(mmsbm.control)] <- mmsbm.control
  ctrl$conv_window <- floor(4 + 1/(ctrl$batch_size))
  set.seed(ctrl$seed)
  if(((ctrl$assortative == FALSE) & (diff(ctrl$mu_block) < 0.0)) | ((ctrl$assortative == TRUE) & (diff(ctrl$mu_block) > 0.0))){
    if(!is.null(mmsbm.control$mu_block)){
      warning("Value of mu_block is not consistent with assortative argument. Switching signs.")
    }
    ctrl$mu_block <- ctrl$mu_block * -1
  } 
  
  
  ## Perform control checks
  if(ctrl$svi){
    if((ctrl$forget_rate <= 0.5) | (ctrl$forget_rate > 1.0)){
      stop("For stochastic VI, forget_rate must be in (0.5,1].")
    }
    if(ctrl$delay < 0.0){
      stop("For stochastic VI, delay must be non-negative.")
    }
  } else {
    ctrl$forget_rate <- 0.0
    ctrl$batch_size <- 1.0
  }
  
  
  
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
  if(anyDuplicated(mfd[,c("(tid)","(sid)","(rid)")])){
    stop("timeID, senderID, and receiverID do not uniquely identify observations in data.dyad.")
  }
  
  
  
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
  colnames(dntid) <- c("(sid)","(rid)")
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
  if(anyDuplicated(mfm[,c("(tid)","(nid)")])){
    stop("timeID and nodeID do not uniquely identify observations in data.monad.")
  }
  ntid <- do.call(paste, c(mfm[c("(nid)","(tid)")], sep="@"))
  if(!all(dntid %in% ntid))
    stop("Nodes in data.dyad missing from data.monad. Are node and time identifiers identical in data.dyad and data.monad?")
  match_ids <- ntid %in% dntid
  if(any(!match_ids)){
    if(ctrl$verbose){
      cat("\tSome nodes in data.monad not present in data.dyad; dropping them.\n")
    }
    mfm <- mfm[match_ids, ]
    ntid <- do.call(paste, c(mfm[c("(nid)","(tid)")], sep="@"))
  }
  if(!is.null(ctrl$fixed_mm)){
    ctrl$node_est <- !(ntid %in% ctrl$fixed_mm) 
  } else{
    ctrl$node_est <- rep(1, length(ntid))
  }
  Y <- stats::model.response(mfd)
  X <- .scaleVars(mfm)
  X_mean <- attr(X, "scaled:center")
  X_sd <- attr(X, "scaled:scale")
  n_monad_pred <- ncol(X)
  
  Z <- .scaleVars(mfd, FALSE)
  Z_mean <- attr(Z, "scaled:center")
  Z_sd <- attr(Z, "scaled:scale")
  n_dyad_pred <- ncol(Z)
 
  
  
  ## Modify prior means and variances to match transformed model matrix
  
  ctrl$mu_gamma <- .transf_muvar(ctrl$mu_gamma, FALSE, FALSE, Z)
  ctrl$var_gamma <- .transf_muvar(ctrl$var_gamma, TRUE, FALSE, Z)
  ctrl$mu_beta <- .transf_muvar(ctrl$mu_beta, FALSE, TRUE, X, n.blocks, n.hmmstates)
  ctrl$var_beta <- .transf_muvar(ctrl$var_beta, TRUE, TRUE, X, n.blocks, n.hmmstates)
  mu_block <- var_block <- array(NA, c(n.blocks, n.blocks))
  diag(mu_block) <- ctrl[["mu_block"]][1]
  mu_block[upper.tri(mu_block)|lower.tri(mu_block)] <- ctrl[["mu_block"]][2]
  diag(var_block) <- ctrl[["var_block"]][1]
  var_block[upper.tri(var_block)|lower.tri(var_block)] <- ctrl[["var_block"]][2]
  
  
  nt_id <- cbind(match(dntid[,1], ntid) - 1, match(dntid[,2], ntid) - 1)
  t_id_d <- match(mfd[["(tid)"]], ut) - 1
  t_id_n <- match(mfm[["(tid)"]], ut) - 1
  nodes_pp <- c(by(mfm, mfm[["(tid)"]], nrow))
  dyads_pp <- c(by(mfd, mfd[["(tid)"]], nrow))
  node_id_period <- split(1:nrow(X), t_id_n)
  
  ## Translate batch size to number of nodes
  if(periods == 1){
    ctrl$batch_size = max(1, floor(ctrl$batch_size * sum(nodes_pp)))
  } else {
    ctrl$batch_size = sapply(nodes_pp, function(x)max(1, floor(ctrl$batch_size * x)))
  }
  
  ## Create initial values
  if(ctrl$verbose){
    cat("Obtaining initial values...\n")
  }
  
  ##Initial HMM states
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
    if(isFALSE(all.equal(colSums(ctrl$kappa_init_t), rep(1.0, ncol(ctrl$kappa_init_t)), check.names = FALSE)) || any(ctrl$kappa_init_t < 0.0)){
      stop("Elements in kappa_init_t must be positive, and its columns must sum to one.")
    }
    ctrl$kappa_init_t <- t(.transf(t(ctrl$kappa_init_t)))
    state_init <- apply(ctrl$kappa_init_t, 2, which.max)
  } 
  if(identical(n.hmmstates, 1)){
    names(state_init) = 1
  }
  
  ##Initial mm 
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
  
  ## Initialize mm
  if(is.null(ctrl$mm_init_t) | !(all(dntid %in% colnames(ctrl$mm_init_t)))){
    mm_init_t <- .initPi(soc_mats,
                        dyads,
                        edges,
                        nodes_pp,
                        dyads_pp,
                        n.blocks, periods, directed, ctrl)[[1]]
    if(!(is.null(ctrl$mm_init_t)) & !(all(dntid %in% colnames(ctrl$mm_init_t)))) {
      sum_mm <- mm_init_t[,colnames(ctrl$mm_init_t)]
      loss.mat <- sum_mm %*% t(ctrl$mm_init_t)
      right.perm <- clue::solve_LSAP(t(loss.mat), TRUE)
      mm_init_t <- mm_init_t[right.perm,]
      mm_init_t[,colnames(ctrl$mm_init_t)] <- ctrl$mm_init_t
    } 
    ctrl$mm_init_t <- mm_init_t
  }
  
  
  


  if(is.null(ctrl$gamma_init)){
    if(ncol(Z) > 0){
      ctrl$gamma_init <- rnorm(length(ctrl$mu_gamma), ctrl$mu_gamma, sqrt(ctrl$var_gamma))
    } else {
      ctrl$gamma_init <- 0
    }
    names(ctrl$gamma_init) <- names(ctrl$mu_gamma)
  }
  if(ncol(Z) == 0)
    Z <- matrix(0, nrow = nrow(Z), ncol = 1)

  if(is.null(ctrl$b_init_t)){
    # ctrl$b_init_t <- qlogis(approxB(Y, nt_id, ctrl$mm_init_t, directed))
    # if(any(is.infinite(ctrl$b_init_t))){
    #   which.inf <- which(is.infinite(ctrl$b_init_t))
    #   ctrl$b_init_t[which.inf] <- ifelse(ctrl$b_init_t[which.inf] > 0, 25, -25)
    # }
    ctrl$b_init_t <- array(rnorm(mu_block, mu_block, sqrt(var_block)), c(n.blocks, n.blocks))
  }
  if(is.null(ctrl$beta_init)){
    prot <- array(.1, dim(ctrl$mu_beta)[-3], dimnames=dimnames(ctrl$mu_beta)[-3])
    ctrl$beta_init <- vapply(seq.int(n.hmmstates),
           function(m){
             array(rnorm(ctrl$mu_beta[,,m], ctrl$mu_beta[,,m],sqrt(ctrl$var_beta[,,m])), dim(ctrl$mu_beta)[-3])
           }, prot)
  }
  # ## Estimate model
  if(ctrl$verbose){
    cat("Estimating model...\n");
  }
  X_t <- t(X)
  Z_t <- t(Z)
  ## For HO sample, 5% (or 10 min) of each type of tie
  ##ho_ind_lo <- sample(seq_len(ncol(Z_t))[Y < 0.5], max(sum(Y < 0.5)*0.05, 10))
  ##ho_ind_hi <- sample(seq_len(ncol(Z_t))[Y >= 0.5], max(sum(Y >= 0.5)*0.05, 10))
  ##ho_ind <- c(ho_ind_lo, ho_ind_hi)
  ##sparsity <- mean(Y >= 0.5)
  fit <- mmsbm_fit(Z_t,
                   ##Z_t[, ho_ind, drop=FALSE],
                   X_t,
                   Y,
                   ##Y[ho_ind, drop=FALSE],
                   t_id_d,##[-ho_ind, drop=FALSE],
                   t_id_n,
                   nodes_pp,
                   nt_id,
                   ##nt_id[ho_ind,, drop=FALSE],
                   node_id_period,
                   mu_block,
                   var_block,
                   ctrl$mu_beta,
                   ctrl$var_beta,
                   ctrl$mu_gamma,
                   ctrl$var_gamma,
                   ctrl$mm_init_t,
                   ctrl$kappa_init_t,
                   ctrl$b_init_t,
                   ctrl$beta_init,
                   ctrl$gamma_init,
                   ##sparsity,
                   ctrl)
  if(!fit[["converged"]])
    warning(paste("Model did not converge after", fit[["niter"]] - 1, "iterations.\n"))
  else if (ctrl$verbose){
    cat("done after", fit[["niter"]] - 1, "iterations.\n")
  }
  
  
  ##Return transposes 
  fit[["TransitionKernel"]] <- t(fit[["TransitionKernel"]])
  fit[["BlockModel"]] <- t(fit[["BlockModel"]])
  
  
  ## Rescale and name coefficients
  fit[["DyadCoef"]] <- c(fit[["DyadCoef"]]) / Z_sd
  if(length(fit[["DyadCoef"]])){
    Z <- t(t(Z) * Z_sd + Z_mean)
    fit[["BlockModel"]] <- fit[["BlockModel"]] - c(Z_mean %*% fit[["DyadCoef"]])
    names(fit[["DyadCoef"]]) <- colnames(Z) 
  }
  
  fit[["MonadCoef"]] <- vapply(1:n.hmmstates,
                               function(ind, coefs, sd_vec, mean_vec){
                                 mat <- coefs[,,ind, drop=FALSE]
                                 constx <- 1
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
    
    fit$vcov_monad <- .vcovBeta(all_phi, fit[["MonadCoef"]], ctrl$se_sim, n.blocks,
                                 n.hmmstates, fit[["TotNodes"]], periods,
                                 ctrl$mu_beta, ctrl$var_beta, fit[["Kappa"]], t_id_n, X_t) 
    
    
    ## and for dyadic coefficients
    z_samples <- replicate(ctrl$se_sim, getZ(fit[["SenderPhi"]]), simplify = FALSE)
    w_samples <- replicate(ctrl$se_sim, getZ(fit[["ReceiverPhi"]]), simplify = FALSE)
    
    if(ctrl$directed){
      all_theta_par<-c(fit[["BlockModel"]], fit[["DyadCoef"]])
    }else{
      all_theta_par<-c(fit[["BlockModel"]][lower.tri(fit[["BlockModel"]], diag = TRUE)], fit[["DyadCoef"]])
    }
    group_mat <- matrix(1:(n.blocks*n.blocks), n.blocks, n.blocks)
    lambda_vec <- c(c(var_block), ctrl$var_gamma) #var_b changed to var_block
    if(!directed){
      group_mat[upper.tri(group_mat)] <- group_mat[lower.tri(group_mat)]
      lambda_vec <- c(c(var_block[lower.tri(var_block, TRUE)]), ctrl$var_gamma)
    } 
    #hessTheta
    hessTheta_list <- mapply(
      function(send_samp, rec_samp, y_vec, Z_d, par_theta, mu_b_mat, var_b_mat, var_g, mu_g, dir_net, group_mat, lambda_vec)
      {
        n_samp <- max(ctrl$dyad_vcov_samp, floor(ncol(Z_d)*0.10))
        samp_ind <- sample(1:ncol(Z_d), n_samp)
        tries <- 0
        if(any(Z_d!=0)){  
          while(any(apply(Z_d[,samp_ind,drop=FALSE], 1, stats::sd) == 0.0) & (tries < 100)){
            samp_ind <- sample(1:ncol(Z_d), n_samp)
            tries <- tries + 1
          }
          
        }
        if(tries >= 100){
          stop("Bad sample for dyadic vcov computation; too little variation in dyadic covariates.")
        }
        group_vec <- model.matrix(~factor(diag(t(send_samp[, samp_ind]) %*% group_mat %*% rec_samp[,samp_ind]), levels = unique(c(group_mat)))-1)
        mod_Z <- group_vec
        if(any(Z_d!=0)){
          mod_Z <- cbind(mod_Z, t(Z_d[,samp_ind, drop=FALSE]))
        }
        if(directed){
          mod_gamma <- c(c(fit$BlockModel),fit$DyadCoef)
        } else {
          mod_gamma <- c(c(fit$BlockModel[lower.tri(fit$BlockModel, TRUE)]),fit$DyadCoef)
        }
        s_eta <- plogis(mod_Z %*% mod_gamma)
        D_mat <- diag(c(s_eta*(1-s_eta)))
        hess_tmp <- ((t(mod_Z) %*% D_mat %*% mod_Z) - diag(1/lambda_vec))*(ncol(Z_d)/n_samp)
        vc_tmp <- Matrix::forceSymmetric(solve(hess_tmp)) 
        ev <- eigen(vc_tmp)$value
        if(any(ev<0)){
          vc_tmp <- vc_tmp - diag(min(ev) - 1e-4, ncol(vc_tmp))
        }
        ch_vc <- chol(vc_tmp)
        return(t(ch_vc) %*% ch_vc)
      },
      z_samples, w_samples,
      MoreArgs = list(par_theta = all_theta_par, 
                      y_vec = Y,
                      Z_d = t(Z),
                      mu_b_mat = mu_block,
                      var_b_mat = var_block,
                      var_g = ctrl$var_gamma, 
                      mu_g = ctrl$mu_gamma,
                      dir_net = ctrl$directed,
                      group_mat = group_mat,
                      lambda_vec = lambda_vec),
      SIMPLIFY = FALSE)
    
    vcovTheta <- Reduce("+", hessTheta_list)/ctrl$se_sim
    N_B_PAR <- ifelse(directed, n.blocks*n.blocks , n.blocks * (1 + n.blocks) / 2)
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
    
  }#end Hessian portion
  
  
  
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
                    n.blocks = n.blocks,
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
