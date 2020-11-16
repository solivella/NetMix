#' Dynamic mixed-membership stochastic blockmodel with covariates
#'
#' The function estimates a dynamic mixed-membership stochastic
#' blockmodel that incorporates covariates. 
#'
#' @param formula.dyad A \code{formula} object. The variable in \code{data.dyad} that contains 
#'     binary edges should be used as a LHS, and any dyadic predictors 
#'     can be included on the RHS (when no dyadic covariates are available, use \code{y ~ 1}).
#'     Same syntax as a \code{glm} formula. 
#' @param formula.monad1 An optional \code{formula} object. LHS is ignored. RHS contains 
#'     names of nodal atrributes found in \code{data.monad1}.  
#' @param formula.monad2 An optional \code{formula} object. LHS is ignored. RHS contains 
#'     names of nodal atrributes found in \code{data.monad2}.  
#' @param senderID Character string. Quoted name of the variable in \code{data.dyad} identifying 
#'     the sender node. For undirected networks, the variable simply contains name of first node 
#'     in dyad. Cannot contain special charecter "`@`". 
#' @param receiverID Character string. Quoted name of the variable in \code{data.dyad} identifying 
#'     the receiver node. For undirected networks, the variable simply contains name of second node 
#'     in dyad. Cannot contain special charecter "`@`".
#' @param nodeID1 Character string. Quoted name of the variable in \code{data.monad1} identifying 
#'     a node in either \code{data.dyad[,senderID]} or \code{data.dyad[,receiverIDD]}. If not \code{NULL},
#'     every node \code{data.dyad[,senderID]} or \code{data.dyad[,senderID]} must be present in 
#'     \code{data.monad1[,nodeID1]}. Cannot contain special charecter "`@`".
#' @param nodeID2 Character string. Quoted name of the variable in \code{data.monad2} identifying 
#'     a node in either \code{data.dyad[,senderID]} or \code{data.dyad[,receiverID]}. If not \code{NULL},
#'     every node \code{data.dyad[,senderID]} or \code{data.dyad[,senderID]} must be present in 
#'     \code{data.monad2[,nodeID2]}. Cannot contain special charecter "`@`".
#' @param timeID Character string. Quoted name of the variable in both \code{data.dyad} and
#'     \code{data.monad} indicating the time in which network (and correspding nodal atrributes)
#'     were observed. The variable itself must be composed of integers. Cannot contain special charecter "`@`".
#' @param data.dyad Data frame. Sociomatrix in ``long'' (i.e. dyadic) format. Must contain at
#'    least three variables: the sender identifier (or identifier of the first node in an undirected networks dyad),
#'    the receiver identifier (or identifier of the second node in an undirected network dyad), and the value
#'    of the edge between them. Currently, only edges between zero and one (inclusive) are supported.
#' @param data.monad1 Data frame for Family 1. Nodal atributes. Must contain a node identifier matching the names of nodes
#'    used in the \code{data.dyad} data frame for nodes of Family 1. 
#' @param data.monad2 Data frame for Family 2. Nodal atributes. Must contain a node identifier matching the names of nodes
#'    used in the \code{data.dyad} data frame for nodes of Family 2. 
#' @param n.blocks1 Integer value. How many latent groups in Family 1 should be used to estimate the model?
#' @param n.blocks2 Integer value. How many latent groups in Family 2 should be used to estimate the model?
#' @param n.hmmstates Integer value. How many hidden Markov state should be used in the HMM? Defaults 
#'    to 1 (i.e. no HMM).  
#' @param directed Boolean. Is the network directed? Defaults to \code{TRUE}.
#' @param missing Means of handling missing data. One of "indicator method" (default) or "listwise deletion".
#' @param mmsbm.control A named list of optional algorithm control parameters.
#'     \describe{
#'        \item{seed}{RNG seed. Defaults to \code{NULL}, which does not seed the RNG}
#'        \item{bipartite}{Boolean; is the network bipartite.}
#'        \item{nstart}{Integer. Number of random initialization trials. Defaults to 5.}            
#'        \item{spectral}{Boolean. Type of initialization algorithm for mixed-membership vectors in static case. If \code{TRUE} (default),
#'                    use spectral clustering with degree correction; otherwise, use kmeans algorithm}
#'        \item{init_gibbs}{Boolean. Should a collapsed Gibbs sampler of non-regression mmsbmB be used to initialize
#'                    each time period? Setting to \code{TRUE} will result in slower initialization and faster model estimation. Setting to \code{TRUE} will be result in faster estimation that is very sensitive to
#'                    choice of alpha (see below)}            
#'        \item{alpha}{Numeric positive value. Concentration parameter for collapsed Gibbs sampler to find initial
#'                     mixed-membership values in dynamic case when \code{init_gibbs=TRUE}. Defaults to 1.0}
#'        \item{missing}{Means of handling missing data. One of "indicator method" (default) or "listwise deletion".}           
#'        \item{svi}{Boolean; should stochastic variational inference be used? Defaults to \code{TRUE}.}   
#'        \item{vi_iter}{Number of maximum iterations in stochastic variational updates. Defaults to 5e2.}
#'        \item{batch_size1}{When \code{svi=TRUE}, proportion of family 1 nodes sampled in each local. Defaults to 0.05 when \code{svi=TRUE}, and to 1.0 otherwise.}
#'        \item{batch_size2}{When \code{svi=TRUE}, proportion of family 2 nodes sampled in each local. Defaults to 0.05 when \code{svi=TRUE}, and to 1.0 otherwise.}
#'        \item{forget_rate}{When \code{svi=TRUE}, value between (0.5,1], controlling speed of decay of weight of prior
#'                            parameter values in global steps. Defaults to 0.75 when \code{svi=TRUE}, and to 0.0 otherwise.}
#'        \item{delay}{When \code{svi=TRUE}, non-negative value controlling weight of past iterations in global steps. Defaults to 1.0 when \code{svi=TRUE},
#'                     and ignored otherwise.}
#'        \item{opt_iter}{Number of maximum iterations of BFGS in global step. Defaults to 10e3.}
#'        \item{hessian}{Boolean indicating whether the Hessian matrix of regression coefficients should e returned. Defaults to \code{TRUE}.}
#'        \item{mu_b}{Numeric vector with two elements: prior mean of blockmodel's main diagonal elements, and
#'                    and prior mean of blockmodel's offdiagonal elements. Defaults to \code{c(5.0, -5.0)}}
#'        \item{var_b}{Numeric vector with two positive elements: prior variance of blockmodel's main diagonal elements, and
#'                    and prior variance of blockmodel's offdiagonal elements. Defaults to \code{c(1.0, 1.0)}}
#'        \item{mu_beta1}{Either single numeric value, in which case the same prior mean is applied to all family 1 monadic coefficients, or
#'                       an array with that is \code{npredictors1} by \code{n.blocks1} by \code{n.hmmstates}, where \code{npredictors1}
#'                       is the number of family 1 monadic predictors for which a prior mean is being set (prior means need not be set for all)
#'                       predictors). The rows in the array should be named to identify which variables a prior mean is being set for.
#'                       Defaults to a common prior mean of 0.0 for all monadic coefficients.}           
#'        \item{var_beta1}{See \code{mu_beta1}. Defaults to a single common prior variance of 1.0 for all family 1 monadic coefficients.}
#'        \item{mu_beta2}{Either single numeric value, in which case the same prior mean is applied to all family 2 monadic coefficients, or
#'                       an array with that is \code{npredictors2} by \code{n.blocks2} by \code{n.hmmstates}, where \code{npredictors2}
#'                       is the number of family 2 monadic predictors for which a prior mean is being set (prior means need not be set for all)
#'                       predictors). The rows in the array should be named to identify which variables a prior mean is being set for.
#'                       Defaults to a common prior mean of 0.0 for all monadic coefficients.}   
#'        \item{var_beta2}{See \code{mu_beta2}. Defaults to a single common prior variance of 1.0 for all family 2 monadic coefficients.}  
#'        \item{mu_gamma}{Either a single numeric value, in which case the same prior mean is applied to all dyadic coefficients, or
#'                        a named vector of numeric values (with names corresponding to the name of the variable 
#'                       for which a prior mean is being set). Defaults to a common prior mean of 0.0 for all dyadic coefficients.}
#'        \item{var_gamma}{See \code{mu_gamma}. Defaults to a single common prior variance of 1.0 for all dyadic coefficients.}
#'        \item{eta}{Numeric positive value. Concentration hyper-parameter for HMM. Defaults to 10.3}
#'        \item{se_sim}{Number of samples from variational posterior of latent variables on which approximation to variance-covariance
#'                      matrices are based. Defaults to 10.}
#'        \item{dyad_vcov_samp}{Number of dyads to sample in computation of variance-covariance of dyadic and blockmodel parameters. 
#'                              Defaults to 1000.}
#'        \item{fixed_mm1}{Optional character vector, with \code{"nodeID1@timeID"} as elements, indicating which mixed-membership vectors
#'                        should remain constant at their initial values throughout estimation. When only one year is observed, elements should be 
#'                         \code{"nodeID1@1"}. Typically used with \code{mm_init_t1}.}  
#'        \item{fixed_mm2}{Optional character vector, with \code{"nodeID2@timeID"} as elements, indicating which mixed-membership vectors
#'                        should remain constant at their initial values throughout estimation. When only one year is observed, elements should be 
#'                         \code{"nodeID2@1"}. Typically used with \code{mm_init_t2}.}
#'        \item{mm_init_t1}{Matrix, \code{n.blocks1} by nodes across years. Optional initial values for mixed-membership vectors.
#'                           Although initial values need not be provided for all nodes, column names must have a \code{nodeID1@timeID} format to 
#'                           avoid ambiguity. When only one year is observed, names should be \code{"nodeID1@1"}.}
#'        \item{mm_init_t2}{Matrix, \code{n.blocks2} by nodes across years. Optional initial values for mixed-membership vectors.
#'                           Although initial values need not be provided for all nodes, column names must have a \code{nodeID2@timeID} format to 
#'                           avoid ambiguity. When only one year is observed, names should be \code{"nodeID2@1"}.}
#'        \item{kappa_init_t}{Matrix, \code{n.hmmstates} by number of years. Optional initial values for variational 
#'                       parameters for state probabilities.}
#'        \item{b_init_t}{Matrix, \code{n.blocks1} by \code{n.blocks2}. Optional initial values for blockmodel.}
#'        \item{beta1_init}{Array, predictors by \code{n.blocks1} by \code{n.hmmstates}. Optional initial values for family 1 monadic coefficients.}
#'        \item{beta2_init}{Array, predictors by \code{n.blocks2} by \code{n.hmmstates}. Optional initial values for family 2 monadic coefficients.}
#'        \item{gamma_init}{Vector. Optional initial values for dyadic coefficients.}
#'        \item{permute}{Boolean. Should all permutations be tested to realign initial block models in dynamic case? If \code{FALSE}, realignment is 
#'                      done via faster graph matching algorithm, but may not be exact. Defaults to \code{TRUE}.}
#'        \item{conv_tol}{Numeric value. Absolute tolerance for VI convergence. Defaults to 1e-4}
#'        \item{verbose}{Boolean. Should extra information be printed as model iterates? Defaults to FALSE}
#'        }
#'       
#' @return Object of class \code{mmsbmB}. List with named components:
#'     \describe{
#'       \item{MixedMembership 1}{Matrix of variational posterior of mean of mixed-membership vectors for family 1. \code{nodes1} by \
#'                              \code{n.groups1}}
#'       \item{MixedMembership 2}{Matrix of variational posterior of mean of mixed-membership vectors for family 2. \code{nodes2} by \
#'                              \code{n.groups2}}
#'       \item{BlockModel}{\code{n.groups1} by \code{n.groups2} matrix of estimated tie log-odds between members
#'                         of corresponding latent groups. The blockmodel.}
#'       \item{vcov_blockmodel}{If \code{hessian=TRUE}, variance-covariance matrix of parameters in blockmodel, ordered in column-major order.}
#'       \item{MonadCoef1}{Array of estimated coefficient values for family 1 monadic covariates. Has \code{n.groups1} columns,
#'                        and \code{n.hmmstates} slices.}
#'       \item{MonadCoef2}{Array of estimated coefficient values for family 2 monadic covariates. Has \code{n.groups2} columns,
#'                        and \code{n.hmmstates} slices.}
#'       \item{DyadCoef}{Vector estimated coefficient values for dyadic covariates}
#'       \item{vcov_dyad}{If \code{hessian=TRUE}, variance-covariance matrix of dyadic coefficients.}
#'       \item{TransitionKernel}{Matrix of estimated HMM transition probabilities}
#'       \item{Kappa}{Matrix of marginal probabilities of being in an HMM state at any given point in time. 
#'                    \code{n.hmmstates} by years (or whatever time interval networks are observed at)}
#'       \item{LowerBound}{Final LB value}
#'       \item{niter}{Final number of VI iterations}
#'       \item{converged}{Convergence indicator; zero indicates failure to converge.}
#'       \item{NodeIndex1}{Order in which family 1 nodes are stored in all return objects.}
#'       \item{NodeIndex2}{Order in which family 2 nodes are stored in all return objects.}
#'       \item{monadic.data1, monadic.data2, dyadic.data, n_states, n_blocks1, n_blocks2, directed}{Original values of parameters used during estimation}
#'     }
#' @author Santiago Olivella (olivella@@unc.edu), Adeline Lo (aylo@@wisc.edu), Tyler Pratt (tyler.pratt@@yale.edu), Kosuke Imai (imai@@harvard.edu)
#' 
#' @example tests/Examples/MIDColdWar.R
#' 

mmsbmB <- function(formula.dyad,
                  formula.monad1=~1,
                  formula.monad2=~1,
                  senderID, 
                  receiverID,
                  nodeID1 = NULL,
                  nodeID2 = NULL,
                  timeID = NULL,
                  data.dyad,
                  data.monad1 = NULL,
                  data.monad2 = NULL,
                  n.blocks1,
                  n.blocks2,
                  n.hmmstates,
                  directed = TRUE,
                  mmsbm.control = list()){
  require(blockcluster)
  require(snowboot)
  #
  cl <- match.call(expand.dots = FALSE)
  formulas1 <- cl[match(c("formula.dyad","formula.monad1"), names(cl))]
  formulas2 <- cl[match(c("formula.dyad","formula.monad2"), names(cl))]
  ## Form default control list
  ctrl <- list(blocks1 = n.blocks1,
               blocks2 = n.blocks2,
               #nodes2 = NULL, ## remove at some point
               states = n.hmmstates,
               times = 1,
               directed = TRUE,
               seed = sample(500,1),
               bipartite = TRUE,
               svi = TRUE,
               nstarts = 5,
               spectral = TRUE,
               init_gibbs = if (n.hmmstates > 1) TRUE else FALSE,
               threads = 1,
               alpha = 1.0,
               forget_rate = 0.75,
               delay = 1.0,
               batch_size1 = 0.05,
               batch_size2 = 0.05,
               missing = "indicator method",
               vi_iter = 500,
               hessian = TRUE,
               se_sim = 10,
               dyad_vcov_samp = 100,
               opt_iter = 10e3,
               mu_b = c(1.0, 0.0),
               var_b = c(1.0, 1.0),
               mu_beta1 = 0.0,
               var_beta1 = 1.0,
               mu_beta2 = 0.0,
               var_beta2 = 1.0,
               mu_gamma = 0.0,
               var_gamma = 1.0,
               #phi1_init_t = NULL,
               #phi2_init_t = NULL,
               #kappa_init_t = NULL,
               #b_init_t = NULL,
               #beta1_init = NULL,
               #beta2_init = NULL,
               gamma_init = NULL,
               #theta_init = NULL,
               eta = 1.0,
               permute = TRUE,
               threads = 1,
               conv_tol = 1e-3,
               verbose = FALSE)
  
  
  directed<-ctrl$directed ## INSERTED HACK 
  
  ctrl[names(mmsbm.control)] <- mmsbm.control
  ctrl$conv_window <- floor(4 + 1/(ctrl$batch_size1)) #currently just for batch size of 1
  set.seed(ctrl$seed)
  
  #if(ctrl$bipartite == FALSE){
    #ctrl$nodes2 <- 0
    #ctrl$npred2 <- 0
  #}else{
    #ctrl$nodes2 <- nodes2
    #ctrl$npred2 <- npred2
  #}
  
  mu_b <- var_b <- array(NA, c(n.blocks2, n.blocks1))
  diag(mu_b) <- ctrl[["mu_b"]][1]
  mu_b[upper.tri(mu_b)|lower.tri(mu_b)] <- ctrl[["mu_b"]][2]
  diag(var_b) <- ctrl[["var_b"]][1]
  var_b[upper.tri(var_b)|lower.tri(var_b)] <- ctrl[["var_b"]][2]
  
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
            identical(class(formula.monad1), "formula"),
            identical(class(formula.monad2), "formula"),
            is.data.frame(data.dyad))
  
  if(!is.null(data.monad1)){
    stopifnot(is.data.frame(data.monad1),
              !is.null(nodeID1))
  }
  if(!is.null(data.monad2)){
    stopifnot(is.data.frame(data.monad2),
              !is.null(nodeID2))
  }
  if(((length(formula.monad1)>2) || (formula.monad1[2] != "1()")) & is.null(data.monad1)){
    stop("Monadic dataset 1 not defined.")
  }
  if(((length(formula.monad2)>2) || (formula.monad2[2] != "1()")) & is.null(data.monad2)){
    stop("Monadic dataset 2 not defined.")
  }
  
  ## Add time variable if null or single period
  if(is.null(timeID) || (length(unique(data.dyad[[timeID]])) == 1)){
    timeID <- "tid"
    data.dyad[timeID] <- 1
    if(!is.null(data.monad1)) {
      data.monad1[timeID] <- 1  
    }
    if(!is.null(data.monad2)) {
      data.monad2[timeID] <- 1  
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
    # monadic dataset 1
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
    # monadic dataset 2
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
  if(identical(ctrl$missing, "listwise deletion")){
    if(length(all.vars(formula.dyad[[3]]))){
      mdyad <- apply(as.matrix(data.dyad[,all.vars(formula.dyad[[3]])]), 1, function(x){!any(is.na(x))})
    } else {
      mdyad <- TRUE
    }
    #monadic 1
    if(length(all.vars(formula.monad1[[2]]))){
      mmonad1 <- apply(as.matrix(data.monad1[,all.vars(formula.monad1[[2]])]), 1, function(x){!any(is.na(x))})
    } else {
      mmonad1 <- TRUE
    }
    #monadic 2
    if(length(all.vars(formula.monad2[[2]]))){
      mmonad2 <- apply(as.matrix(data.monad2[,all.vars(formula.monad2[[2]])]), 1, function(x){!any(is.na(x))})
    } else {
      mmonad2 <- TRUE
    }
    data.dyad <- data.dyad[mdyad,]
    data.monad1 <- data.monad1[mmonad1,]
    data.monad2 <- data.monad2[mmonad2,]
    d.keep <- lapply(unique(data.dyad[,timeID]), function(x){
      nts1 <- data.monad1[data.monad1[,timeID]==x,nodeID1]
      nts2 <- data.monad2[data.monad2[,timeID]==x,nodeID2]
      dd <- data.dyad[data.dyad[,timeID]==x,]
      dd <- dd[dd[,senderID] %in% nts1 & dd[,receiverID] %in% nts2,]
      return(dd)
    })
    data.dyad <- do.call("rbind", d.keep)
  }
  
  mfd <- do.call(model.frame, list(formula = formula.dyad,
                                   data = data.dyad,
                                   drop.unused.levels = TRUE,
                                   tid = as.name(timeID),
                                   sid = as.name(senderID),
                                   rid = as.name(receiverID))) #ok for differential names
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
  dntid1 <- dntid[,1]
  dntid2 <- dntid[,2]
  udnid1 <- unique(unlist(mfd[c("(sid)")]))
  udnid2 <- unique(unlist(mfd[c("(rid)")]))
  #if no monadic 1 data entered
  if(is.null(data.monad1)){
    data.monad1 <- data.frame(nid1 = rep(udnid1, periods))
    nodeID1 <- "nid1"
    data.monad1[timeID] <- rep(ut, each = length(udnid1))
  }
  #if no monadic 2 data entered
  if(is.null(data.monad2)){
    data.monad2 <- data.frame(nid2 = rep(udnid2, periods))
    nodeID2 <- "nid2"
    data.monad2[timeID] <- rep(ut, each = length(udnid2))
  }
  
  #Monadic data 1: mfm1
  mfm1 <- do.call(model.frame, list(formula = formula.monad1,
                                    data = data.monad1,
                                    drop.unused.levels = TRUE,
                                    tid = as.name(timeID),
                                    nid1 = as.name(nodeID1)))
  mfm1[,"(nid1)"] <- as.character(mfm1[,"(nid1)"]) 
  if(anyDuplicated(mfm1[,c("(tid)","(nid1)")])){
    stop("timeID and nodeID1 do not uniquely identify observations in data.monad1.")
  }
  
  ntid1 <- do.call(paste, c(mfm1[c("(nid1)","(tid)")], sep="@"))
  if(!all(dntid1 %in% ntid1))
    stop("Nodes in dyadic dataset missing from monadic dataset 1. Are node and time identifiers identical in data.dyad and data.monad1?")
  match_ids <- ntid1 %in% dntid1
  if(any(!match_ids)){
    if(ctrl$verbose){
      cat("\tSome nodes in data.monad1 not present in data.dyad; dropping them.\n")
    }
    mfm1 <- mfm1[match_ids, ]
    ntid1 <- do.call(paste, c(mfm1[c("(nid1)","(tid)")], sep="@"))
  }
  #Monadic data 2: mfm2
  mfm2 <- do.call(model.frame, list(formula = formula.monad2,
                                    data = data.monad2,
                                    drop.unused.levels = TRUE,
                                    tid = as.name(timeID),
                                    nid2 = as.name(nodeID2)))
  mfm2[,"(nid2)"] <- as.character(mfm2[,"(nid2)"]) 
  if(anyDuplicated(mfm2[,c("(tid)","(nid2)")])){
    stop("timeID and nodeID2 do not uniquely identify observations in data.monad2.")
  }
  ntid2 <- do.call(paste, c(mfm2[c("(nid2)","(tid)")], sep="@"))
  if(!all(dntid2 %in% ntid2))
    stop("Nodes in dyadic dataset missing from monadic dataset 2. Are node and time identifiers identical in data.dyad and data.monad2?")
  match_ids <- ntid2 %in% dntid2
  if(any(!match_ids)){
    if(ctrl$verbose){
      cat("\tSome nodes in data.monad2 not present in data.dyad; dropping them.\n")
    }
    mfm2 <- mfm2[match_ids, ]
    ntid2 <- do.call(paste, c(mfm2[c("(nid2)","(tid)")], sep="@"))
  }
  if(!is.null(ctrl$fixed_mm1)){
    ctrl$node_est1 <- !(ntid1 %in% ctrl$fixed_mm1) 
  } else{
    ctrl$node_est1 <- rep(1, length(ntid1))
  }
  if(!is.null(ctrl$fixed_mm2)){
    ctrl$node_est2 <- !(ntid2 %in% ctrl$fixed_mm2) 
  } else{
    ctrl$node_est2 <- rep(1, length(ntid2))
  }
  
  
  
  Y <- stats::model.response(mfd)
  X1 <- base::scale(model.matrix(terms(mfm1), mfm1))
  X1_mean <- attr(X1, "scaled:center")
  X1_sd <- attr(X1, "scaled:scale")
  if(any(X1_sd==0)){
    constx <- which(X1_sd==0)
    X1[,constx] <- 1
    X1_sd[constx]<-1##ADDED to remove X1_sd elements that ==0 and set 1
  }
  n_monad1_pred <- ncol(X1)
  X2 <- base::scale(model.matrix(terms(mfm2), mfm2))
  X2_mean <- attr(X2, "scaled:center")
  X2_sd <- attr(X2, "scaled:scale")
  if(any(X2_sd==0)){
    constx <- which(X2_sd==0)
    X2[,constx] <- 1
    X2_sd[constx]<-1##ADDED to remove X2_sd elements that ==0 and set 1
  }
  n_monad2_pred <- ncol(X2)
  Z <- scale(model.matrix(terms(mfd), mfd))
  Z_mean <- attr(Z, "scaled:center")
  Z_sd <- attr(Z, "scaled:scale")
  if(any(Z_sd==0)){
    constz <- which(Z_sd==0)
    Z <- Z[,-constz, drop = FALSE]
  }
  n_dyad_pred <- ncol(Z)
  
  ## Modify prior means and variances to match transformed model matrix
  ctrl$mu_gamma <- .transf_muvar(ctrl$mu_gamma, FALSE, FALSE, Z)
  ctrl$var_gamma <- .transf_muvar(ctrl$var_gamma, TRUE, FALSE, Z, devs=Z_sd[-which(Z_sd==0)])
  ctrl$mu_beta1 <- .transf_muvar(ctrl$mu_beta1, FALSE, TRUE, X1, n.blocks1, n.hmmstates)
  ctrl$mu_beta2 <- .transf_muvar(ctrl$mu_beta2, FALSE, TRUE, X2, n.blocks2, n.hmmstates)
  ctrl$var_beta1 <- .transf_muvar(ctrl$var_beta1, TRUE, TRUE, X1, n.blocks1, n.hmmstates, c(1,X1_sd[-1]))
  ctrl$var_beta2 <- .transf_muvar(ctrl$var_beta2, TRUE, TRUE, X2, n.blocks2, n.hmmstates, c(1,X2_sd[-1]))
  ## (older) Modify prior means and variances to match transformed model matrix
  #ctrl$mu_gamma <- .transf_muvar(ctrl$mu_gamma, FALSE, FALSE, Z)
  #ctrl$var_gamma <- .transf_muvar(ctrl$var_gamma, TRUE, FALSE, Z, devs=Z_sd[-which(Z_sd==0)])
  #ctrl$mu_beta1 <- .transf_muvar(ctrl$mu_beta1, FALSE, TRUE, X1, n.blocks1, n.hmmstates) 
  #ctrl$mu_beta2 <- .transf_muvar(ctrl$mu_beta2, FALSE, TRUE, X2, n.blocks2, n.hmmstates) 
  #ctrl$var_beta1 <- .transf_muvar(ctrl$var_beta1, TRUE, TRUE, X1, n.blocks1, n.hmmstates, c(1,X1_sd[-1]))
  #ctrl$var_beta1[1,,]<-0.5 #### just for the intercepts, you could define a standard deviation of 0.5
  #ctrl$var_beta2 <- .transf_muvar(ctrl$var_beta2, TRUE, TRUE, X2, n.blocks2, n.hmmstates, c(1,X2_sd[-1]))
  #ctrl$var_beta2[1,,]<-0.5 #### just for the intercepts, you could define a standard deviation of 0.5
  
  nt_id <- cbind(match(dntid[,1], ntid1) - 1, match(dntid[,2], ntid2) - 1) ## here's where it matters the swap in node id
  nt_id1 <- nt_id[,1]
  nt_id2 <- nt_id[,2]
  t_id_d <- match(mfd[["(tid)"]], ut) - 1
  t_id_n1 <- match(mfm1[["(tid)"]], ut) - 1
  t_id_n2 <- match(mfm2[["(tid)"]], ut) - 1
  nodes_pp <- c(c(by(mfm1, mfm1[["(tid)"]], nrow)),c(by(mfm2, mfm2[["(tid)"]], nrow)))
  nodes_pp1<- c(by(mfm1, mfm1[["(tid)"]], nrow))
  nodes_pp2<- c(by(mfm2, mfm2[["(tid)"]], nrow))
  dyads_pp <- c(by(mfd, mfd[["(tid)"]], nrow))
  node_id_period1 <- split(1:nrow(X1), t_id_n1)
  node_id_period2 <- split(1:nrow(X2), t_id_n2)
  
  ## Translate batch size to number of nodes
  if(periods == 1){
    ctrl$batch_size1 = max(1, floor(ctrl$batch_size1 * sum(nodes_pp1)))
    ctrl$batch_size2 = max(1, floor(ctrl$batch_size2 * sum(nodes_pp2)))
  } else {
    ctrl$batch_size1 = sapply(nodes_pp1, function(x)max(1, floor(ctrl$batch_size1 * x)))
    ctrl$batch_size2 = sapply(nodes_pp2, function(x)max(1, floor(ctrl$batch_size2 * x)))
  }
  ## Create initial values
  if(ctrl$verbose){
    cat("Obtaining initial values...\n")
  }
  
  if(!ctrl$bipartite){
    all.nodes <- unique(unlist(mfd[,c("(sid)","(rid)")]))
  }else{
    all.nodes1 <- unique(unlist(mfd[,"(sid)"]))
    all.nodes2 <- unique(unlist(mfd[,"(rid)"]))
  }
  node.cols <- which(names(mfd)%in%c("(sid)","(rid)", "(tid)"))
  
  dyads <- split.data.frame(mfd[,c(node.cols, 1)], mfd[, "(tid)"])
  edges <- split(Y, mfd[, "(tid)"])
  #create sociomatrix
  #if(ctrl$bipartite){soc_mats <- .createSocioB(dyads,all.nodes1,all.nodes2, ctrl$directed)}else{
    #soc_mats <- .createSocioB(dyads,all.nodes,all.nodes, ctrl$directed)}
  soc_mats <- Map(function(dyad_mat, edge_vec){
    #nodes <- unique(c(dyad_mat))
    nodes1 <- unique(dyad_mat[,which(names(dyad_mat)%in%"(sid)")])
    nodes2 <- unique(dyad_mat[,which(names(dyad_mat)%in%"(rid)")])
    nnode1 <- length(nodes1)
    nnode2 <- length(nodes2)
    adj_mat <- matrix(NA,
                      nnode1,
                      nnode2,
                      dimnames = list(nodes1,
                                      nodes2))
    #adj_mat[dyad_mat] <- edge_vec
    indeces <- as.matrix(dyad_mat[,c("(sid)","(rid)")])
    index<-cbind(match(indeces[,1],rownames(adj_mat)),match(indeces[,2],colnames(adj_mat)))
    adj_mat[index] <- dyad_mat[,which(names(dyad_mat)%in%"Y")] # out of bounds
    if(!directed){
      #adj_mat[dyad_mat[,c(2,1)]] <- edge_vec
      adj_mat[index[,c(2,1)]] <- dyad_mat[,which(names(dyad_mat)%in%"Y")]
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
  
  
  #null Kappa
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
                                  iter.max = 10,
                                  nstart = ctrl$nstarts), "classes")
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

  #null mm_init1 and mm_init2 (previously Phi 1 and 2): needs debugging
  # if(is.null(ctrl$mm_init_t1)|is.null(ctrl$mm_init_t2)){
  #   ctrl$mm_init_t <- .initPi(soc_mats,
  #                             Y,
  #                             dyads,
  #                             edges,
  #                             t_id_d,
  #                             t_id_n1,t_id_n2,
  #                             nodes_pp,nodes_pp1,nodes_pp2,
  #                             dyads_pp,
  #                             nt_id1,nt_id2,
  #                             node_id_period1,node_id_period2,
  #                             mu_b,
  #                             var_b,
  #                             nrow(Z), n.blocks1, n.blocks2, periods, directed, ctrl) #22
  # } else{
  #   if(isFALSE(all.equal(colSums(ctrl$mm_init_t), rep(1.0, ncol(ctrl$mm_init_t)), check.names = FALSE)) || any(ctrl$mm_init_t < 0.0)){
  #     stop("Elements in mm_init_t must be positive, and its columns must sum to one.")
  #   }
  #   ctrl$mm_init_t <- t(.transf(t(ctrl$mm_init_t)))
  # } 
  
  #null mm_init1 and mm_init2 (previously Phi 1 and 2):
  ## Adjusted this for when edges are never formed for a given node/year
    if((is.null(ctrl$mm_init_t1)|is.null(ctrl$mm_init_t2)) & (periods == 1)){
      phi_init_temp <- lapply(soc_mats, function(mat){
        clust.o<-coclusterBinary(mat,nbcocluster=c(n.blocks1,n.blocks2))
        phi1_init_temp<-matrix(0,nrow=nrow(mat),ncol=n.blocks1)
        phi2_init_temp<-matrix(0,nrow=ncol(mat),ncol=n.blocks2)
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
      ctrl$mm_init_t1 <- do.call(cbind, lapply(phi_init_temp, `[[`, 1))#1st matrix of each element of big list
      ctrl$mm_init_t2 <- do.call(cbind, lapply(phi_init_temp, `[[`, 2))#2nd matrix of each element of big list
    } 

  #null gamma
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
  
  #null B: currently using only first item in phi1_init_t/phi2_init_t -- how to condense that to two matrices only?
  if(is.null(ctrl$b_init_t)){
    cat("Estimating initial blockmodel values...\n")
    ctrl$b_init_t<-qlogis(approxB(Y, nt_id, ctrl$mm_init_t1, ctrl$mm_init_t2))
    if(any(is.infinite(ctrl$b_init_t))|any(is.nan(ctrl$b_init_t))){
      cat("Conducting bootstrapped initialization of blockmodel...\n")
      #set.seed(1)
      B <- 5 #set number of bootstraps
      ## Note AL: currently using first item in sociomatrix only!
      #II. Draw bootstrap sample from mat, keep track of rows, create phi matrices - "Phi'"
      mat <- soc_mats[[1]]
      tmp.n1<-nrow(mat)
      tmp.n2<-ncol(mat)
      mat.B <- .vertboot2(mat, B)
      #set up holders
      phi1<-t(ctrl$mm_init_t1) #t(blk x node) -> node x blk
      phi2<-t(ctrl$mm_init_t2) 
      block.boot<-vector(mode="list", length=(B+1))
      block.boot[[1]]<-ctrl$b_init_t # as first item orig sociomatrix created blockmodel; g x h
      if(any(is.infinite(block.boot[[1]]))){
        cat("Inf element in B matrix; replace... \n ")
        which.inf <- which(is.infinite(block.boot[[1]]))
        block.boot[[1]][which.inf] <- ifelse(block.boot[[1]][which.inf] > 0, 25, -25)
      }
      if(any(is.nan(block.boot[[1]]))){
        cat("NaN element in B matrix; replace... \n ")
        which.nan <- which(is.nan(block.boot[[1]]))
        block.boot[[1]][which.nan] <- 0
      }
      
      for(i in 1: B){
        phi1.prime<-matrix(0,nrow=tmp.n1,ncol=n.blocks1)
        phi2.prime<-matrix(0,nrow=tmp.n2,ncol=n.blocks2)
        cat("B iter=",i,"\n")
        clust.b<-coclusterBinary(mat.B[[i]]$x,nbcocluster=c(n.blocks1,n.blocks2))
        for(j in 1:tmp.n1){#blk x n
          val<-clust.b@rowclass[j]+1
          phi1.prime[j,val]<-1
        }
        for(j in 1:tmp.n2){
          val<-clust.b@colclass[j]+1
          phi2.prime[j,val]<-1
        }
        #III. Subset "Phi" and "Phi'" to only overlapping rows --> "sub.Phi", "sub.Phi'"
        sub.phi1<-phi1[unique(mat.B[[i]]$index1)+1,] #note index returned is in cpp, starts from 0; n' x blk
        sub.phi2<-phi2[unique(mat.B[[i]]$index2)+1,] # n' x blk
        sub.phi1.prime<-phi1.prime[!duplicated(mat.B[[i]]$index1+1),] #take first time unique row appears; n' x blk
        sub.phi2.prime<-phi2.prime[!duplicated(mat.B[[i]]$index2+1),] # n' x blk
        #IV. Realign sub.Phi to sub.Phi' using Hungarian algorithm, switch labels for Phi and Phi'.
        loss <- t(sub.phi1) %*% sub.phi1.prime
        labels<-clue::solve_LSAP(t(loss), TRUE)
        phi1<-phi1[,labels]
        phi1.prime<-phi1.prime[labels,]
        loss <- t(sub.phi2) %*% sub.phi2.prime
        labels<-clue::solve_LSAP(t(loss), TRUE)
        phi2<-phi2[,labels]
        phi2.prime<-phi2.prime[,labels]
        #V. 1st bootstrap compute approxB on Phi and Phi', store B, B'. 
        # All other bootstraps, compute approxB on aligned Phi', store B'.
        # Replace NaNs in B/B' with 0s.
        block.boot[[i+1]]<-qlogis(approxB(Y, nt_id, t(phi1.prime), t(phi2.prime))) #output is K2 x K1
        if(any(is.infinite(block.boot[[i+1]]))){
          cat("Bootstrap=",i,": Inf found in blockmodel!\n")
          which.inf <- which(is.infinite(block.boot[[i+1]]))
          block.boot[[i+1]][which.inf] <- 0#ifelse(block.boot[[i+1]][which.inf] > 0, 25, -25)
        }
        if(any(is.nan(block.boot[[i+1]]))){
          cat("Bootstrap=",i,": NaN found in blockmodel!\n")
          which.nan <- which(is.nan(block.boot[[i+1]]))
          block.boot[[i+1]][which.nan] <- 0
        }
        print(block.boot[[i]])
        #VI. Repeat II-V for set number of bootstraps.
      }
      #VII. Average across B, B's. Return bootstrapped average B.
      #avg.block<-Reduce("+", block.boot) / length(block.boot)
      avg.block<-matrix(0,nrow=n.blocks2,ncol=n.blocks1)
      for(r in 1:n.blocks1){
        for(c in 1:n.blocks2){
          avg<-rep(0,B)
          for(b in 1:B){
            avg[b]<-block.boot[[b]][c,r] #block.boot is k2 x k1
          }
          avg.block[c,r]<-ifelse(length(na.omit(avg))>0,mean(na.omit(avg)),0)
        }
      }
      ctrl$b_init_t<-avg.block
      
    }
  }
  
  
  #null Beta 1
  if(is.null(ctrl$beta1_init)){
    X_state <- split.data.frame(X1, state_init[t_id_n1 + 1])
    phi_state <- split.data.frame(t(ctrl$mm_init_t1), state_init[t_id_n1 + 1])
    ctrl$beta1_init <- mapply(function(X_sub, phi_sub){
      phi_temp <- .transf(phi_sub)
      lm.fit(X_sub, log(phi_temp))$coefficients},
      X_state, phi_state,SIMPLIFY = "array")
  }
  if(anyNA(ctrl$beta1_init)){
    stop("Nearly singular design matrix in family 1; check monadic predictors.")
  }

  #null Beta 2
  if(is.null(ctrl$beta2_init)){
    X_state <- split.data.frame(X2, state_init[t_id_n2 + 1])
    phi_state <- split.data.frame(t(ctrl$mm_init_t2), state_init[t_id_n2 + 1])
    ctrl$beta2_init <- mapply(function(X_sub, phi_sub){
      phi_temp <- .transf(phi_sub)
      lm.fit(X_sub, log(phi_temp))$coefficients},
      X_state, phi_state,SIMPLIFY = "array")
  }
  if(anyNA(ctrl$beta2_init)){
    stop("Nearly singular design matrix in family 2; check monadic predictors.")
  }
  
  # null Theta
  #if(!is.null(true.theta)){
    #ctrl$theta_init<-true.theta
  #}else{
    #if(ctrl$bipartite){
      #ctrl$theta_init<-array(rep(0, n.blocks1*n.blocks2*length(Y)), c(n.blocks2,n.blocks1,length(Y))) #dim theta: h,g,dyad
    #}else{
      #ctrl$theta_init<-rep(0, n.blocks1*n.blocks1*length(Y), c(n.blocks1,n.blocks1,length(Y)))
    #}
  #}
  
  ## Create randomizer for order of updatePhis
  ctrl$phi_order<-rbinom(nrow(Z)[1],1,0.5) #ndyad
  #test print
  #print(ctrl$b_init_t)
  #print(ctrl$beta2_init)
  
  
  ## Estimate model
  if(ctrl$verbose){
    cat("Estimating model...\n");
  }
  X1_t <- t(X1)
  X2_t <- t(X2)
  Z_t <- t(Z)
  ## For HO sample, 5% (or 10 min) of each type of tie
  ##ho_ind_lo <- sample(seq_len(ncol(Z_t))[Y < 0.5], max(sum(Y < 0.5)*0.05, 10))
  ##ho_ind_hi <- sample(seq_len(ncol(Z_t))[Y >= 0.5], max(sum(Y >= 0.5)*0.05, 10))
  ##ho_ind <- c(ho_ind_lo, ho_ind_hi)
  ##sparsity <- mean(Y >= 0.5)
  ## Estimate model
  fit <- mmsbm_fit(Z_t, # 1x100
                   ##Z_t[, ho_ind, drop=FALSE],
                   X1_t, #K1 x N1
                   X2_t, #K2 x N2
                   Y, #length N_DYAD
                   ##Y[ho_ind, drop=FALSE],
                   t_id_d,#length N_DYAD
                   ##[-ho_ind, drop=FALSE],
                   t_id_n1,#length N1
                   t_id_n2,#length N2
                   nodes_pp,#length time*2
                   nodes_pp1,#length time
                   nodes_pp2,#length time
                   nt_id,#2 x N_DYAD
                   ##nt_id[ho_ind,, drop=FALSE],
                   node_id_period1,# time x N1
                   node_id_period2,# time x N2
                   mu_b, #K1 x K2
                   var_b,#K1 x K2
                   ctrl$mu_beta1,
                   ctrl$var_beta1,#predictors1+1 x K1
                   ctrl$mu_beta2,#predictors2+1 x K2
                   ctrl$var_beta2,#predictors2+1 x K2
                   as.vector(ctrl$mu_gamma),#ctrl$mu_gamma,#array
                   as.vector(ctrl$var_gamma),#ctrl$var_gamma,#array
                   ctrl$mm_init_t1, #K1 x N1
                   ctrl$mm_init_t2, #K2 x N2
                   ctrl$kappa_init_t, #time x state
                   ctrl$b_init_t,#K1 x K2
                   ctrl$beta1_init,# predictors1+1 x K1 x time
                   ctrl$beta2_init,# predictors2+1 x K2 x time
                   ctrl$gamma_init,
                   ctrl
  )
  ##Return transposes 
  fit[["TransitionKernel"]] <- t(fit[["TransitionKernel"]])
  fit[["BlockModel"]] <- t(fit[["BlockModel"]])
  
  ## Rescale and name coefficients
  fit[["DyadCoef"]] <- fit[["DyadCoef"]] / Z_sd[-which(Z_sd==0)]
  if(length(fit[["DyadCoef"]])>0){
    Z <- t(t(Z) * Z_sd[-1] + Z_mean[-1])
    fit[["BlockModel"]] <- fit[["BlockModel"]] - c(Z_mean[-constz] %*% fit[["DyadCoef"]])
    names(fit[["DyadCoef"]]) <- colnames(Z) 
  }
  fit[["MonadCoef1"]] <- vapply(1:n.hmmstates,
                               function(ind, coefs, sd_vec, mean_vec){
                                 mat <- coefs[,,ind, drop=FALSE]
                                 constx <- which(sd_vec==0)
                                 mat[-constx, , 1] <- mat[-constx, , 1] / sd_vec[-constx]
                                 if(length(constx)!=0){
                                   mat[constx, ,1] <- mat[constx, ,1] - mean_vec[-constx] %*% mat[-constx, , 1]
                                 }
                                 return(mat)
                               },
                               array(0.0, c(ncol(X1), n.blocks1)),
                               coefs = fit[["MonadCoef1"]],
                               sd_vec = X1_sd,
                               mean_vec = X1_mean)
  rownames(fit[["MonadCoef1"]]) <- colnames(X1)
  colnames(fit[["MonadCoef1"]]) <- paste("Group",1:n.blocks1)
  X1 <- t(t(X1) * X1_sd + X1_mean) #unscale
  fit[["MonadCoef2"]] <- vapply(1:n.hmmstates,
                                function(ind, coefs, sd_vec, mean_vec){
                                  mat <- coefs[,,ind, drop=FALSE]
                                  constx <- which(sd_vec==0)
                                  mat[-constx, , 1] <- mat[-constx, , 1] / sd_vec[-constx]
                                  if(length(constx)!=0){
                                    mat[constx, ,1] <- mat[constx, ,1] - mean_vec[-constx] %*% mat[-constx, , 1]
                                  }
                                  return(mat)
                                },
                                array(0.0, c(ncol(X2), n.blocks2)),
                                coefs = fit[["MonadCoef2"]],
                                sd_vec = X2_sd,
                                mean_vec = X2_mean)
  rownames(fit[["MonadCoef2"]]) <- colnames(X2)
  colnames(fit[["MonadCoef2"]]) <- paste("Group",1:n.blocks2)
  X2 <- t(t(X2) * X2_sd + X2_mean) #unscale
  
  ## Add other names
  colnames(fit[["Kappa"]]) <- unique(mfm1[,"(tid)"])
  dimnames(fit[["BlockModel"]]) <- c(replicate(1,paste("1 Group",1:n.blocks1), simplify = FALSE),replicate(1,paste("2 Group",1:n.blocks2), simplify = FALSE))
  dimnames(fit[["TransitionKernel"]]) <- replicate(2,paste("State",1:n.hmmstates), simplify = FALSE)
  colnames(fit[["MixedMembership 1"]]) <- ntid1
  colnames(fit[["MixedMembership 2"]]) <- ntid2 
  
  if(ctrl$hessian){
    if(ctrl$verbose){
      cat("Computing approximate vcov. matrices...\n")
    }
    ## Compute approximate standard errors
    ## for monadic coefficients
    
  all_phi1 <- split.data.frame((t(fit[["SenderPhi"]])),
                              c(nt_id[,1]))
  all_phi2 <- split.data.frame((t(fit[["ReceiverPhi"]])),
                               c(nt_id[,2]))

  sampleC1_perm <- lapply(all_phi1,
                         function(mat){
                           apply(mat, 2, function(vec)poisbinom::rpoisbinom(ctrl$se_sim, vec))
                         })
  sampleC2_perm <- lapply(all_phi2,
                          function(mat){
                            apply(mat, 2, function(vec)poisbinom::rpoisbinom(ctrl$se_sim, vec))
                          })
  sampleC1_perm <- cbind(do.call(rbind, sampleC1_perm), # samples
                        rep(1:length(all_phi1), each = ctrl$se_sim), #node id
                        rep(1:ctrl$se_sim, times = length(all_phi1))) #sample id
  sampleC2_perm <- cbind(do.call(rbind, sampleC2_perm), # samples
                         rep(1:length(all_phi2), each = ctrl$se_sim), #node id
                         rep(1:ctrl$se_sim, times = length(all_phi2))) #sample id
  sampleC1_perm <- sampleC1_perm[order(sampleC1_perm[,n.blocks1 + 2], sampleC1_perm[,n.blocks1 + 1]),]
  sampleC2_perm <- sampleC2_perm[order(sampleC2_perm[,n.blocks2 + 2], sampleC2_perm[,n.blocks2 + 1]),]
  C1_samples <- split.data.frame(sampleC1_perm[,1:n.blocks1], sampleC1_perm[,n.blocks1 + 2])
  C2_samples <- split.data.frame(sampleC2_perm[,1:n.blocks2], sampleC2_perm[,n.blocks2 + 2])
  S_samples <- replicate(ctrl$se_sim, apply(fit[["Kappa"]], 2, function(x)sample(1:n.hmmstates, 1, prob = x)), simplify = FALSE)

  #hessBeta1
    hessBeta1_list <- mapply( function(C_samp, S_samp, tidn, X_i, Nvec, beta_vec, vbeta, mbeta, periods)
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
        hess_tmp <- optimHess(c(beta_vec),alphaLB,
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
      C1_samples, S_samples,
      MoreArgs = list(tidn = t_id_n1,
                      X_i = t(X1),
                      Nvec = fit[["TotNodes1"]],
                      beta_vec = fit[["MonadCoef1"]],
                      vbeta = ctrl$var_beta1,
                      mbeta = ctrl$mu_beta1,
                      periods = periods),
      SIMPLIFY=FALSE)
    #hessBeta2
    hessBeta2_list <- mapply( function(C_samp, S_samp, tidn, X_i, Nvec, beta_vec, vbeta, mbeta, periods)
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
      hess_tmp <- optimHess(c(beta_vec),alphaLB,
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
    C2_samples, S_samples,
    MoreArgs = list(tidn = t_id_n2,
                    X_i = t(X2),
                    Nvec = fit[["TotNodes2"]],
                    beta_vec = fit[["MonadCoef2"]],
                    vbeta = ctrl$var_beta2,
                    mbeta = ctrl$mu_beta2,
                    periods = periods),
    SIMPLIFY=FALSE)
     
    fit$vcov_monad1 <- Reduce("+", hessBeta1_list)/ctrl$se_sim
    fit$vcov_monad2 <- Reduce("+", hessBeta2_list)/ctrl$se_sim
    
    colnames(fit$vcov_monad1) <- rownames(fit$vcov_monad1) <- paste(rep(paste("State",1:n.hmmstates), each = prod(dim(fit[["MonadCoef1"]])[1:2])),
                                                                  rep(colnames(fit[["MonadCoef1"]]), each = nrow(fit[["MonadCoef1"]]), times = n.hmmstates),
                                                                  rep(rownames(fit[["MonadCoef1"]]), times = n.blocks1*n.hmmstates),
                                                                  sep=":") 
    colnames(fit$vcov_monad2) <- rownames(fit$vcov_monad2) <- paste(rep(paste("State",1:n.hmmstates), each = prod(dim(fit[["MonadCoef2"]])[1:2])),
                                                                    rep(colnames(fit[["MonadCoef2"]]), each = nrow(fit[["MonadCoef2"]]), times = n.hmmstates),
                                                                    rep(rownames(fit[["MonadCoef2"]]), times = n.blocks2*n.hmmstates),
                                                                    sep=":") 
    ## and for dyadic coefficients
    z_samples <- replicate(ctrl$se_sim, getZ(fit[["SenderPhi"]]), simplify = FALSE)
    w_samples <- replicate(ctrl$se_sim, getZ(fit[["ReceiverPhi"]]), simplify = FALSE)
    
    if(ctrl$directed){
      all_theta_par<-c(fit[["BlockModel"]], fit[["DyadCoef"]])
    }else{
      all_theta_par<-c(fit[["BlockModel"]][lower.tri(fit[["BlockModel"]], diag = TRUE)], fit[["DyadCoef"]])
    }
    group_mat <- matrix(1:(n.blocks1*n.blocks2), n.blocks1, n.blocks2)
    lambda_vec <- c(c(var_b), ctrl$var_gamma)
    if(!directed){
      group_mat[upper.tri(group_mat)] <- group_mat[lower.tri(group_mat)]
      lambda_vec <- c(c(var_b[lower.tri(var_b, TRUE)]), ctrl$var_gamma)
    } 
    #hessTheta
    hessTheta_list <- mapply(
      function(send_samp, rec_samp, y_vec, Z_d, par_theta, mu_b_mat, var_b_mat, var_g, mu_g, dir_net, group_mat, lambda_vec)
      {
        n_samp <- min(ctrl$dyad_vcov_samp, floor(ncol(Z_d)*0.25))
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
                      mu_b_mat = mu_b,
                      var_b_mat = var_b,
                      var_g = ctrl$var_gamma, 
                      mu_g = ctrl$mu_gamma,
                      dir_net = ctrl$directed,
                      group_mat = group_mat,
                      lambda_vec = lambda_vec),
      SIMPLIFY = FALSE)
    
    vcovTheta <- Reduce("+", hessTheta_list)/ctrl$se_sim
    N_B_PAR <- ifelse(directed, n.blocks1*n.blocks2 , n.blocks1 * (1 + n.blocks2) / 2)
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
  
    
  #Include used data
  attr(mfm1, "terms") <- NULL
  attr(mfm2, "terms") <- NULL
  attr(mfd, "terms") <- NULL
  #
  fit$monadic1.data <- mfm1 #nodes are in nid1,nide2 naming conventions; time in tid convention
  fit$monadic2.data <- mfm2
  fit$dyadic.data <- mfd
  fit$Y <- Y
  
  ## Include node id's
  fit$NodeIndex <- nt_id
  
  ## Include a few formals needed by other methods
  fit$forms <- list(directed = directed,
                    senderID = senderID,
                    receiverID = receiverID,
                    timeID = timeID,
                    nodeID1 = nodeID1,
                    nodeID2 = nodeID2,
                    t_id_d = t_id_d,
                    hessian = ctrl$hessian,
                    formula.dyad = formula.dyad,
                    formula.monad1 = formula.monad1,
                    formula.monad2 = formula.monad2
                    #nodes2 = nodes2,
                    #npred2 = npred2
                    #formula.dyad = formulas1[[1]],
                    #formula.monad1 = formulas1[[2]],
                    #formula.monad2 = formulas2[[2]]
                    )
  
  ## Include used seed
  fit$seed <- ctrl$seed
  
  ## Include original call
  fit$call <- match.call()
  
  fit$bipartite<-ifelse(ctrl$bipartite,TRUE,FALSE)
  
  ##Assign class for methods
  if(fit$bipartite){
  class(fit) <-"mmsbmB"
  }else{
  class(fit) <- "mmsbm"
  }
  return(fit)
}

