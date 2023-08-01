#' Dynamic mixed-membership stochastic blockmodel with covariates
#'
#' The function estimates a dynamic mixed-membership stochastic
#' blockmodel that incorporates covariates. 
#'
#' @param formula.dyad A \code{formula} object. The variable in \code{data.dyad} that contains 
#'     binary edges should be used as a LHS, and any dyadic predictors 
#'     can be included on the RHS (when no dyadic covariates are available, use \code{y ~ 1}).
#'     Same syntax as a \code{glm} formula. 
#' @param formula.monad Optional. Either a \code{formula} object (when \code{bipartite=FALSE}) or
#'     a list of two \code{formula} objects (when \code{bipartite=FALSE}). In both instances, formula RHS
#'    containing names of nodal attributes found in \code{data.monad}, and LHS is ignored. 
#'     names of nodal attributes found in \code{data.monad1}. 
#' @param bipartite Boolean. Is the network bipartite? Defaults to \code{FALSE}. 
#' @param senderID Character string. Quoted name of the variable in \code{data.dyad} identifying 
#'     the sender node. For undirected networks, the variable simply contains name of first node 
#'     in dyad. Cannot contain special charecter "`@`". 
#' @param receiverID Character string. Quoted name of the variable in \code{data.dyad} identifying 
#'     the receiver node. For undirected networks, the variable simply contains name of second node 
#'     in dyad. Cannot contain special charecter "`@`".
#' @param nodeID Either character string (when \code{bipartite = FALSE}) or list of two character strings (\code{bipartite = TRUE}).
#'     In both cases, the strings are quoted name(s) of the variable(s) in \code{data.monad} identifying 
#'     a node in either \code{data.dyad[,senderID]} or \code{data.dyad[,receiverID]}. If not \code{NULL},
#'     every node \code{data.dyad[,senderID]} or \code{data.dyad[,senderID]} must be present in 
#'     one of the dataframes contained in \code{data.monad}. Names cannot contain special character "`@`".
#' @param timeID Either character string (when \code{bipartite = FALSE}) or list of two character strings (\code{bipartite = TRUE}). 
#'  In both cases, the strings are quoted name(s) of the variable(s) in both \code{data.dyad} and the elements in
#'     \code{data.monad} indicating the time in which network (and corresponding nodal attributes)
#'     were observed. The variable itself must be composed of integers. Names cannot contain special character "`@`".
#' @param data.dyad Data frame. Sociomatrix in ``long'' (i.e. dyadic) format. Must contain at
#'    least three variables: the sender identifier (or identifier of the first node in an undirected networks dyad),
#'    the receiver identifier (or identifier of the second node in an undirected network dyad), and the value
#'    of the edge between them. Currently, only edge values between zero and one (inclusive) are supported.
#' @param data.monad Either data.frame (when \code{bipartite = FALSE}) or list of two data.frames strings (\code{bipartite = TRUE}).
#'  Nodal attributes. Must contain a node identifier matching the names of nodes
#'    used in the \code{data.dyad} data frame for nodes. 
#' @param n.blocks Integer value (when \code{bipartite = FALSE}) or vector of two integers. How many latent groups (in either family type) should be used to estimate the model?
#' @param n.hmmstates Integer value. How many hidden Markov state should be used in the HMM? Defaults 
#'    to 1 (i.e. no HMM). Currently, dynamic models are only supported for non-bipartite models.  
#' @param directed Boolean. Is the network directed? Defaults to \code{TRUE}.
#' @param missing Means of handling missing data. One of "indicator method" (default) or "listwise deletion".
#' @param mmsbm.control A named list of optional algorithm control parameters.
#'     \describe{
#'        \item{seed}{RNG seed. Defaults to \code{NULL}, which does not seed the RNG}
#'        \item{nstart}{Integer. Number of random initialization trials. Defaults to 5.}            
#'        \item{spectral}{Boolean. Type of initialization algorithm for mixed-membership vectors in static case. If \code{TRUE} (default),
#'                    use spectral clustering with degree correction; otherwise, use kmeans algorithm}
#'        \item{init_gibbs}{Boolean. Should a collapsed Gibbs sampler of non-regression mmsbmB be used to initialize
#'                    each time period? Setting to \code{TRUE} will result in slower initialization and faster model estimation. Setting to \code{TRUE} will be result in faster estimation that is very sensitive to
#'                    choice of alpha (see below)}            
#'        \item{alpha}{Numeric positive value. Concentration parameter for collapsed Gibbs sampler to find initial
#'                     mixed-membership values in dynamic case when \code{init_gibbs=TRUE}. Defaults to 1.0}
#'        \item{missing}{Means of handling missing data. One of "indicator method" (default) or "listwise deletion".}   
#'        \item{assortative}{Boolean indicating whether blockmodel should be assortative (i.e. stronger connections within groups) or disassortative
#'                           (i.e. stronger connections between groups). Defaults to \code{TRUE}.}        
#'        \item{svi}{Boolean; should stochastic variational inference be used? Defaults to \code{TRUE}.}   
#'        \item{vi_iter}{Number of maximum iterations in stochastic variational updates. Defaults to 5e2.}
#'        \item{batch_size}{Numeric vector. When \code{svi=TRUE}, proportion of family 1 and family 2 nodes sampled in each local. Defaults to c(0.05, 0.05) when \code{svi=TRUE}, and to c(1.0, 1.0) otherwise.}
#'        \item{forget_rate}{When \code{svi=TRUE}, value between (0.5,1], controlling speed of decay of weight of prior
#'                            parameter values in global steps. Defaults to 0.75 when \code{svi=TRUE}, and to 0.0 otherwise.}
#'        \item{delay}{When \code{svi=TRUE}, non-negative value controlling weight of past iterations in global steps. Defaults to 1.0 when \code{svi=TRUE},
#'                     and ignored otherwise.}
#'        \item{opt_iter}{Number of maximum iterations of BFGS in global step. Defaults to 10e3.}
#'        \item{hessian}{Boolean indicating whether the Hessian matrix of regression coefficients should e returned. Defaults to \code{TRUE}.}
#'        \item{mu_block}{Numeric vector with two elements: prior mean of blockmodel's main diagonal elements, and
#'                    and prior mean of blockmodel's offdiagonal elements. Defaults to \code{c(5.0, -5.0)}}
#'        \item{var_block}{Numeric vector with two positive elements: prior variance of blockmodel's main diagonal elements, and
#'                    and prior variance of blockmodel's offdiagonal elements. Defaults to \code{c(1.0, 1.0)}}
#'        \item{mu_beta}{List with one or two elements, depending on \code{bipartite}. Each element must be either single numeric value, in which case the same prior mean is applied to all family 1 monadic coefficients, or
#'                       an array with that is \code{npredictors} by \code{n.blocks} by \code{n.hmmstates}, where \code{npredictors}
#'                       is the number of  monadic predictors in the corresponding family for which a prior mean is being set (prior means need not be set for all 
#'                       predictors). The rows in the array should be named to identify which variables a prior mean is being set for.
#'                       Defaults to a common prior mean of 0.0 for all monadic coefficients.}           
#'        \item{var_beta}{See \code{mu_beta}. Defaults to a single common prior variance of 2.55 for all monadic coefficients, and 10 for all intercepts.}
#'        \item{mu_gamma}{Either a single numeric value, in which case the same prior mean is applied to all dyadic coefficients, or
#'                        a named vector of numeric values (with names corresponding to the name of the variable 
#'                       for which a prior mean is being set). Defaults to a common prior mean of 0.0 for all dyadic coefficients.}
#'        \item{var_gamma}{See \code{mu_gamma}. Defaults to a single common prior variance of 5.0 for all dyadic coefficients.}
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
#'        \item{mm_init_t1}{Matrix, \code{n.blocks[1]} by nodes across years. Optional initial values for mixed-membership vectors.
#'                           Although initial values need not be provided for all nodes, column names must have a \code{nodeID1@timeID} format to 
#'                           avoid ambiguity. When only one year is observed, names should be \code{"nodeID1@1"}.}
#'        \item{mm_init_t2}{Matrix, \code{n.blocks[2]} by nodes across years. Optional initial values for mixed-membership vectors.
#'                           Although initial values need not be provided for all nodes, column names must have a \code{nodeID2@timeID} format to 
#'                           avoid ambiguity. When only one year is observed, names should be \code{"nodeID2@1"}.}
#'        \item{kappa_init_t}{Matrix, \code{n.hmmstates} by number of years. Optional initial values for variational 
#'                       parameters for state probabilities.}
#'        \item{block_init_t}{Matrix, \code{n.blocks[1]} by \code{n.blocks2}. Optional initial values for blockmodel.}
#'        \item{beta1_init}{Array, predictors by \code{n.blocks[1]} by \code{n.hmmstates}. Optional initial values for family 1 monadic coefficients.}
#'        \item{beta2_init}{Array, predictors by \code{n.blocks[2]} by \code{n.hmmstates}. Optional initial values for family 2 monadic coefficients.}
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
                  bipartite = FALSE,
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
  
  if(bipartite){
    if(identical(formula.monad, ~1)){
      formula.monad <- list(formula.monad, formula.monad)
    }else if(!is.list(formula.monad) | (length(formula.monad) != 2)){
      stop("When bipartite is TRUE, formula.monad must be a list of formulas of length two.")
    }
    if(!is.null(data.monad) & (!is.list(data.monad) | (length(data.monad)!=2))){
      stop("When bipartite is TRUE, data.monad must be a list of data.frame's of length two.")
    }
    if(is.null(nodeID)){
      nodeID <- list(NULL, NULL)
    } else if(!is.list(nodeID) | (length(nodeID)!=2)){
      stop("When bipartite is TRUE, nodeID must be a list of character vectors of length two.")
    }
    if(length(n.blocks) != 2){
      stop("n.blocks must be an integer vector of length 2.")
    }
  } else{
    formula.monad[[2]] <- 0
    data.monad <- list(data.monad)
    data.monad[[2]] <- 0
    nodeID <- list(nodeID)
    n.blocks <- c(n.blocks, n.blocks)
  }
  
  cl <- match.call(expand.dots = FALSE)
  
  ## Form default control list
  ctrl <- list(blocks1 = n.blocks[1],
               blocks2 = n.blocks[2],
               states = n.hmmstates,
               times = 1,
               seed = sample(500,1),
               svi = TRUE,
               nstarts = 5,
               spectral = TRUE,
               init_gibbs = if (n.hmmstates > 1) TRUE else FALSE,
               threads = 1,
               alpha = 1.0,
               forget_rate = 0.75,
               delay = 1.0,
               batch_size = c(0.05, 0.05),
               missing = "indicator method",
               vi_iter = 500,
               hessian = TRUE,
               se_sim = 10,
               dyad_vcov_samp = 1000,
               opt_iter = 10e3,
               mu_block = c(5.0, -5.0),
               var_block = c(5.0, 5.0),
               mu_beta = list(0.0, 0.0),
               var_beta = list(5.0, 5.0),
               mu_gamma = 0.0,
               var_gamma = 5.0,
               mm1_init_t = NULL,
               mm2_init_t = NULL,
               kappa_init_t = NULL,
               b_init_t = NULL,
               assortative = TRUE,
               beta1_init = NULL,
               beta2_init = NULL,
               gamma_init = NULL,
               eta = 1.0,
               permute = TRUE,
               threads = 1,
               conv_tol = 1e-2,
               verbose = FALSE)
  ctrl[names(mmsbm.control)] <- mmsbm.control
  ctrl$bipartite <- bipartite
  ctrl$directed <- ifelse(!bipartite,directed,TRUE) #patch currently, since not doing directed bipartite
  ctrl$conv_window <- floor(4 + 1/(ctrl$batch_size[1])) #currently just for batch size of 1
  set.seed(ctrl$seed)
  if(((ctrl$assortative == FALSE) & (diff(ctrl$mu_block) < 0.0)) | ((ctrl$assortative == TRUE) & (diff(ctrl$mu_block) > 0.0))){
    if(!is.null(mmsbm.control$mu_block)){
      warning("Value of mu_block is not consistent with assortative argument. Changing signs.")
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
    if(length(ctrl$batch_size) == 1){
      ctrl$batch_size[2] <- ctrl$batch_size
    }
  } else {
    ctrl$forget_rate <- 0.0
    ctrl$batch_size <- c(1.0, 1.0)
  }
  
  if(ctrl$verbose){
    cat("Pre-processing data...\n")
  }
  
  ## Add time variable if null or single period
  if(is.null(timeID) || (length(unique(data.dyad[[timeID]])) == 1)){
    timeID <- "(tid)"
    data.dyad[timeID] <- 1
    if(!is.null(data.monad[[1]])) {
      data.monad[[1]][timeID] <- 1  
    }
    if(bipartite){
      if(!is.null(data.monad[[2]])) {
        data.monad[[2]][timeID] <- 1  
      }
    }
  }
  
  ## Address missing data 
  if(any(is.na(data.monad[[1]]))|any(is.na(data.monad[[2]]))|any(is.na(data.dyad))){
    new_dat_dyad <- .missHandle(formula.dyad, data.dyad, ctrl$missing)
    data.dyad <- new_dat_dyad$dat
    formula.dyad <- new_dat_dyad$form
    
    new_dat_monad1 <- .missHandle(formula.monad[[1]], data.monad[[1]], ctrl$missing)
    data.monad[[1]] <- new_dat_monad1$dat
    formula.monad1 <- new_dat_monad1$form
    if(bipartite){
      new_dat_monad2 <- .missHandle(formula.monad[[2]], data.monad[[2]], ctrl$missing)
      data.monad[[2]] <- new_dat_monad2$dat
      formula.monad2 <- new_dat_monad2$form
    } 
    
    ## Drop dyads with nodes not in monadic dataset
    
    if(!is.null(data.monad[[1]])){
      if(bipartite){
        d.keep <- lapply(unique(data.dyad[,timeID]), function(x){
          nts1 <- data.monad[[1]][data.monad[[1]][,timeID]==x,nodeID[[1]]]
          dd <- data.dyad[data.dyad[,timeID]==x,]
          dd <- dd[(dd[,senderID] %in% nts1),]
          return(dd)
        })
        if(!is.null(data.monad[[2]])){
          d.keep <- lapply(unique(data.dyad[,timeID]), function(x){
            nts1 <- data.monad[[1]][data.monad[[1]][,timeID]==x,nodeID[[1]]]
            nts2 <- data.monad[[2]][data.monad[[2]][,timeID]==x,nodeID[[2]]]
            dd <- data.dyad[data.dyad[,timeID]==x,]
            dd <- dd[(dd[,senderID] %in% nts1) & (dd[,receiverID] %in% nts2),]
            return(dd)
          })
        }
      } else {
        d.keep <- lapply(unique(data.dyad[,timeID]), function(x){
          nts <- data.monad[[1]][data.monad[[1]][,timeID]==x,nodeID[[1]]]
          dd <- data.dyad[data.dyad[,timeID]==x,]
          dd <- dd[(dd[,senderID] %in% nts) & (dd[,receiverID] %in% nts),]
          return(dd)
        })
      }
      data.dyad <- do.call("rbind", d.keep)
    }
  }
  ##Create model frames
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
  if(bipartite){
    dntid1 <- dntid[,1]
    dntid2 <- dntid[,2]
    udnid1 <- unique(unlist(mfd[c("(sid)")]))
    udnid2 <- unique(unlist(mfd[c("(rid)")]))
  } else {
    dntid1 <- dntid
    udnid1 <- unique(unlist(mfd[c("(sid)","(rid)")]))
  }
  
  #if no monadic 1 data entered
  if(is.null(data.monad[[1]])){
    data.monad <- list(data.monad)
    data.monad[[1]] <- data.frame("(nid)" = rep(udnid1, periods), check.names = FALSE)
    nodeID[[1]] <- "(nid)"
    data.monad[[1]][timeID] <- rep(ut, each = length(udnid1))
  }
  if(bipartite){
    #if no monadic 2 data entered
    if((length(data.monad) == 1) || is.null(data.monad[[2]])){
      data.monad[[2]] <- data.frame("(nid)" = rep(udnid2, periods), check.names = FALSE)
      nodeID[[2]] <- "(nid)"
      data.monad[[2]][timeID] <- rep(ut, each = length(udnid2))
    }
  }
  
  #Monadic data 1: mfm1
  mfm_tmp1 <- .monadData(formula.monad[[1]], data.monad[[1]], timeID, nodeID[[1]], dntid1, ctrl$verbose) 
  mfm1 <- mfm_tmp1$mf
  ntid1 <- mfm_tmp1$id
  ntid2 <- ntid1
  
  if(bipartite){
    #Monadic data 2: mfm2
    mfm_tmp2 <- .monadData(formula.monad[[2]], data.monad[[2]], timeID, nodeID[[2]], dntid2, ctrl$verbose) 
    mfm2 <- mfm_tmp2$mf
    ntid2 <- mfm_tmp2$id
  }
  
  ## Define fixed mixed-memberships 
  if(!is.null(ctrl$fixed_mm[[1]])){
    ctrl$node_est1 <- !(ntid1 %in% ctrl$fixed_mm[[1]]) 
  } else{
    ctrl$node_est1 <- rep(1, length(ntid1))
  }
  if(bipartite){
    if(!is.null(ctrl$fixed_mm[[2]])){
      ctrl$node_est2 <- !(ntid2 %in% ctrl$fixed_mm[[2]]) 
    } else{
      ctrl$node_est2 <- rep(1, length(ntid2))
    }
  }
  
  
  Y <- stats::model.response(mfd)
  
  X1 <- .scaleVars(mfm1)
  X1_mean <-attr(X1, "scaled:center")
  X1_sd <- attr(X1, "scaled:scale")
  n_monad1_pred <- ncol(X1)
  if(bipartite){
    X2 <- .scaleVars(mfm2)
    X2_mean <- attr(X2, "scaled:center")
    X2_sd <-  attr(X2, "scaled:scale")
    n_monad2_pred <- ncol(X2)
  }
  
  Z <- .scaleVars(mfd, FALSE)
  Z_mean <- attr(Z, "scaled:center")
  Z_sd <- attr(Z, "scaled:scale")
  n_dyad_pred <- ncol(Z)
  #if(n_dyad_pred == 0){
  #Z <- matrix(0, nrow = nrow(Z), ncol = 1)
  #}
  
  ctrl$mu_gamma <- .transf_muvar(ctrl$mu_gamma, FALSE, FALSE, Z)
  ctrl$var_gamma <- .transf_muvar(ctrl$var_gamma, TRUE, FALSE, Z)
  ctrl$mu_beta1 <- .transf_muvar(ctrl$mu_beta[[1]], FALSE, TRUE, X1, n.blocks[1], n.hmmstates)
  ctrl$var_beta1 <- .transf_muvar(ctrl$var_beta[[1]], TRUE, TRUE, X1, n.blocks[1], n.hmmstates)
  if(bipartite){
    ctrl$mu_beta2 <- .transf_muvar(ctrl$mu_beta[[2]], FALSE, TRUE, X2, n.blocks[2], n.hmmstates)
    ctrl$var_beta2 <- .transf_muvar(ctrl$var_beta[[2]], TRUE, TRUE, X2, n.blocks[2], n.hmmstates)
  }
  
  ## Create full blockmodel mean and variance priors
  mu_block <- var_block <- array(NA, c(n.blocks[2], n.blocks[1]))
  diag(mu_block) <- ctrl[["mu_block"]][1]
  mu_block[upper.tri(mu_block)|lower.tri(mu_block)] <- ctrl[["mu_block"]][2]
  diag(var_block) <- ctrl[["var_block"]][1]
  var_block[upper.tri(var_block)|lower.tri(var_block)] <- ctrl[["var_block"]][2]
  
  
  ## Define node id's
  nt_id <- cbind(match(dntid[,1], ntid1) - 1, match(dntid[,2], ntid2) - 1) ## here's where it matters the swap in node id
  t_id_n1 <- match(mfm1[["(tid)"]], ut) - 1
  node_id_period1 <- split(1:nrow(X1), t_id_n1)
  nodes_pp <- nodes_pp1 <- c(by(mfm1, mfm1[["(tid)"]], nrow))
  if(bipartite){
    nt_id1 <- nt_id[,1]
    nt_id2 <- nt_id[,2]
    t_id_n2 <- match(mfm2[["(tid)"]], ut) - 1
    nodes_pp <- c(c(by(mfm1, mfm1[["(tid)"]], nrow)),c(by(mfm2, mfm2[["(tid)"]], nrow)))
    nodes_pp2<- c(by(mfm2, mfm2[["(tid)"]], nrow))
    node_id_period2 <- split(1:nrow(X2), t_id_n2)
  } 
  t_id_d <- match(mfd[["(tid)"]], ut) - 1
  dyads_pp <- c(by(mfd, mfd[["(tid)"]], nrow))
  
  ## Translate batch size to number of nodes
  if(periods == 1){
    ctrl$batch_size1 <- max(1, floor(ctrl$batch_size[1] * sum(nodes_pp1)))
    if(bipartite){
      ctrl$batch_size2 <- max(1, floor(ctrl$batch_size[2] * sum(nodes_pp2)))
    }
  } else {
    ctrl$batch_size1 <- sapply(nodes_pp1, function(x)max(1, floor(ctrl$batch_size[1] * x)))
    if(bipartite){
      ctrl$batch_size2 <- sapply(nodes_pp2, function(x)max(1, floor(ctrl$batch_size[2] * x)))
    }
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
  
  if(!bipartite){
    all.nodes <- unique(unlist(mfd[,c("(sid)","(rid)")]))
  }else{
    all.nodes1 <- unique(unlist(mfd[,"(sid)"]))
    all.nodes2 <- unique(unlist(mfd[,"(rid)"]))
  }
  
  node.cols <- which(names(mfd) %in% c("(sid)","(rid)", "(tid)"))
  
  dyads <- split.data.frame(mfd[,c(node.cols, 1)], mfd[, "(tid)"])
  edges <- split(Y, mfd[, "(tid)"])
  
  #create sociomatrices
  soc_mats <- Map(function(dyad_mat, edge_vec, bipartite, y_var = all.vars(formula.dyad)[1]){
    #nodes <- unique(c(dyad_mat))
    if(bipartite){
      nodes1 <- unique(as.vector(dyad_mat[,which(names(dyad_mat)=="(sid)")]))
      nodes2 <- unique(as.vector(dyad_mat[,which(names(dyad_mat)=="(rid)")]))
      nnode1 <- length(nodes1)
      nnode2 <- length(nodes2)
    } else {
      nodes1 <- nodes2 <- unique(unlist(dyad_mat[,c("(sid)","(rid)")]))
      nnode1 <- nnode2 <- length(nodes1)
    }
    adj_mat <- matrix(NA,
                      nnode1,
                      nnode2,
                      dimnames = list(nodes1,
                                      nodes2))
    indeces <- as.matrix(dyad_mat[,c("(sid)","(rid)")])
    index <- cbind(match(indeces[,1],rownames(adj_mat)),match(indeces[,2],colnames(adj_mat)))
    adj_mat[index] <- dyad_mat[,y_var] 
    if(!directed){
      adj_mat[index[,c(2,1)]] <- dyad_mat[,y_var]
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
  }, dyads, edges, MoreArgs = list(bipartite=bipartite))
  
  ##Initial mm 
  mm_init <- .initPi(soc_mats,
                     bipartite,
                     dyads,
                     edges,
                     nodes_pp,
                     dyads_pp,
                     n.blocks, periods, directed, ctrl)
  ctrl$mm_init_t[[1]] <- mm_init[[1]]
  if(bipartite){
    ctrl$mm_init_t[[2]] <- mm_init[[2]]
  }
  
  
  
  ##Initial gamma
  if(is.null(ctrl$gamma_init)){
    if(n_dyad_pred > 0){
      ctrl$gamma_init <- rnorm(length(ctrl$mu_gamma), ctrl$mu_gamma, sqrt(ctrl$var_gamma))
    } else {
      ctrl$gamma_init <- 0
    }
    names(ctrl$gamma_init) <- names(ctrl$mu_gamma)
  }
  
  ## 
  if(n_dyad_pred == 0){
    Z <- matrix(0, nrow = nrow(Z), ncol = 1)
  }
  
  ##Initial Blockmodel
  if(is.null(ctrl$block_init_t)){
    ctrl$block_init_t <- array(rnorm(mu_block, mu_block, sqrt(var_block)), c(n.blocks[2], n.blocks[1]))
  }
  
  
  ##Initial Beta 1
  if(is.null(ctrl$beta1_init)){
    prot <- array(.1, dim(ctrl$mu_beta1)[-3], dimnames=dimnames(ctrl$mu_beta1)[-3])
    ctrl$beta1_init <- vapply(seq.int(n.hmmstates),
                              function(m){
                                lm.fit(X1, t(ctrl$mm_init_t[[1]]))$coefficients
                              }, prot)
  }
  ##Initial Beta 2
  if(bipartite){
    if(is.null(ctrl$beta2_init)){
      prot <- array(.1, dim(ctrl$mu_beta2)[-3], dimnames=dimnames(ctrl$mu_beta2)[-3])
      ctrl$beta2_init <- vapply(seq.int(n.hmmstates),
                                function(m){
                                  lm.fit(X2,t(ctrl$mm_init_t[[2]]))$coefficients
                                }, prot)
    }
  }
  ## Create randomizer for order of updatePhis
  ctrl$phi_order <- rbinom(nrow(Z)[1],1,0.5) #ndyad
  #test print
  #print(ctrl$b_init_t)
  #print(ctrl$beta2_init)
  
  
  ## Estimate model
  if(ctrl$verbose){
    cat("Estimating model...\n");
  }
  X1_t <- t(X1)
  Z_t <- t(Z)
  
  ## Estimate model
  if(bipartite){
    X2_t <- t(X2)
    node_id_period1<-unlist(node_id_period1)
    node_id_period2<-unlist(node_id_period2)
    dimnames(ctrl$beta1_init) <- NULL
    dimnames(ctrl$beta2_init) <- NULL
    fit <- mmsbm_fitBi(Z_t, # matrix/array 1x100
                       X1_t,# matrix/array K1 x N1
                       X2_t,# matrix/array K2 x N2
                       Y, #length N_DYAD
                       t_id_d,#numeric: length N_DYAD
                       t_id_n1,#numeric: length N1
                       t_id_n2,#numeric: length N2
                       nodes_pp,#integer: length time*2
                       nodes_pp1,#integer: length time
                       nodes_pp2,#integer: length time
                       nt_id,# matrix/array 2 x N_DYAD
                       node_id_period1,#integer: time x N1
                       node_id_period2,#integer: time x N2
                       mu_block, # matrix/array K2 x K1
                       var_block,# matrix/array K2 x K1
                       ctrl$mu_beta1,#array predictors1+1 x K1 x time
                       ctrl$var_beta1,#array predictors1+1 x K1 x time
                       ctrl$mu_beta2,#array predictors2+1 x K2 x time
                       ctrl$var_beta2,#array predictors2+1 x K2 x time
                       as.vector(ctrl$mu_gamma),#vector ctrl$mu_gamma
                       as.vector(ctrl$var_gamma),#vector ctrl$var_gamma
                       ctrl$mm_init_t[[1]], # matrix/array K1 x N1
                       ctrl$mm_init_t[[2]], # matrix/array K2 x N2
                       ctrl$kappa_init_t, # matrix/array time x state
                       ctrl$block_init_t,#matrix/array K2 x K1
                       ctrl$beta1_init,#array predictors1+1 x K1 x time
                       ctrl$beta2_init,#array predictors2+1 x K2 x time
                       ctrl$gamma_init,#numeric vector 
                       ctrl
    )
    
  } else {
    dimnames(ctrl$beta1_init) <- NULL
    fit <- mmsbm_fit(Z_t,
                     X1_t,
                     Y,
                     t_id_d,
                     t_id_n1,
                     nodes_pp1,
                     nt_id,
                     node_id_period1,
                     mu_block,
                     var_block,
                     ctrl$mu_beta1,
                     ctrl$var_beta1,
                     ctrl$mu_gamma,
                     ctrl$var_gamma,
                     ctrl$mm_init_t[[1]],
                     ctrl$kappa_init_t,
                     ctrl$block_init_t,
                     ctrl$beta1_init,
                     ctrl$gamma_init,
                     ctrl)
  }
  if(!fit[["converged"]]){
    warning(paste("Model did not converge after", fit[["niter"]] - 1, "iterations.\n"))
  } else if (ctrl$verbose){
    cat("done after", fit[["niter"]] - 1, "iterations.\n")
  }
  ##Return transposes 
  fit[["TransitionKernel"]] <- t(fit[["TransitionKernel"]])
  fit[["BlockModel"]] <- t(fit[["BlockModel"]])
  
  ## Rescale and name coefficients
  fit[["DyadCoef"]] <- fit[["DyadCoef"]] / Z_sd
  if(length(fit[["DyadCoef"]])>0){
    Z <- t(t(Z) * Z_sd + Z_mean) #unscale
    fit[["BlockModel"]] <- fit[["BlockModel"]] - c(Z_mean %*% fit[["DyadCoef"]]) #recenter block model fixed effects
    names(fit[["DyadCoef"]]) <- colnames(Z) 
  }
  X1 <- t(t(X1) * X1_sd + X1_mean) #unscale
  # tmp1 <- .transfBeta(fit[["MonadCoef1"]], n.hmmstates,
  #                    X1_mean, X1_sd, n.blocks[1], colnames(X1))
  # tmp1[1,,] <- 1
  # tmp_beta1 <- maxLik::maxNR(alphaLBound,
  #                       start = tmp1,
  #                       tot_nodes = fit[["TotNodes1"]],
  #                       c_t=t(fit[["CountMatrix1"]]),
  #                       x_t=t(X1),
  #                       s_mat=fit[["Kappa"]],
  #                       t_id=t_id_n1,
  #                       var_beta=ctrl$var_beta1,
  #                       mu_beta=ctrl$mu_beta1, 
  #                       control=list(iterlim=100))
  # fit[["MonadCoef1"]] <- array(tmp_beta1$estimate, dim(fit[["MonadCoef1"]]))
  # rownames(fit[["MonadCoef1"]]) <- colnames(X1)
  # colnames(fit[["MonadCoef1"]]) <- paste("Group", 1:n.blocks[1])
  fit[["MonadCoef1"]] <- .transfBeta(fit[["MonadCoef1"]], n.hmmstates,
                                     X1_mean, X1_sd, n.blocks[1], colnames(X1))
  
  if(bipartite){
    X2 <- t(t(X2) * X2_sd + X2_mean) #unscale
    # tmp2 <- .transfBeta(fit[["MonadCoef2"]], n.hmmstates,
    #                     X2_mean, X2_sd, n.blocks[2], colnames(X2))
    # tmp2[1,,] <- 1
    # tmp_beta2 <- maxLik::maxNR(alphaLBound,
    #                            start = tmp2,
    #                            tot_nodes = fit[["TotNodes2"]],
    #                            c_t=t(fit[["CountMatrix2"]]),
    #                            x_t=t(X2),
    #                            s_mat=fit[["Kappa"]],
    #                            t_id=t_id_n2,
    #                            var_beta=ctrl$var_beta2,
    #                            mu_beta=ctrl$mu_beta2,
    #                            control=list(iterlim=1))
    # fit[["MonadCoef2"]] <- array(tmp_beta2$estimate, dim(fit[["MonadCoef2"]]))
    # rownames(fit[["MonadCoef2"]]) <- colnames(X2)
    # colnames(fit[["MonadCoef2"]]) <- paste("Group", 1:n.blocks[2])
    fit[["MonadCoef2"]] <- .transfBeta(fit[["MonadCoef2"]], n.hmmstates,
                                       X2_mean, X2_sd, n.blocks[2], colnames(X2))
  }
  
  
  ## Add other names
  colnames(fit[["Kappa"]]) <- unique(mfm1[,"(tid)"])
  dimnames(fit[["BlockModel"]]) <- c(replicate(1,paste("1 Group",1:n.blocks[1]), simplify = FALSE),replicate(1,paste("2 Group",1:n.blocks[2]), simplify = FALSE))
  dimnames(fit[["TransitionKernel"]]) <- replicate(2,paste("State",1:n.hmmstates), simplify = FALSE)
  colnames(fit[["MixedMembership1"]]) <- ntid1
  if(bipartite){
    colnames(fit[["MixedMembership2"]]) <- ntid2 
  }
  
  if(ctrl$hessian){
    if(ctrl$verbose){
      cat("Computing approximate vcov. matrices...\n")
    }
    ## Compute approximate standard errors
    ## for monadic coefficients
    # kappa_mat1 <- t(fit[["Kappa"]][,t_id_n1+1, drop=FALSE])
    # all_phi1 <- (fit[["CountMatrix1"]])
    fit$vcov_monad1 <- .vcovBeta(fit[["MonadCoef1"]],
                                 tot_nodes = fit[["TotNodes1"]],
                                 c_t=t(fit[["CountMatrix1"]]),
                                 x_t=t(X1),
                                 s_mat=fit[["Kappa"]],
                                 t_id=t_id_n1,
                                 var_beta=ctrl$var_beta1,
                                 mu_beta=ctrl$mu_beta1)
    
    if(bipartite){
      # kappa_mat2 <- t(fit[["Kappa"]][,t_id_n2+1, drop=FALSE])
      # all_phi2 <- (fit[["CountMatrix2"]])
      fit$vcov_monad2 <- .vcovBeta(fit[["MonadCoef2"]],
                                   tot_nodes = fit[["TotNodes2"]],
                                   c_t=t(fit[["CountMatrix2"]]),
                                   x_t=t(X2),
                                   s_mat=fit[["Kappa"]],
                                   t_id=t_id_n2,
                                   var_beta=ctrl$var_beta2,
                                   mu_beta=ctrl$mu_beta2)
      
    } 
    
    ## and for dyadic coefficients
    if(any(Z_sd > 0)){
      edge_eta <- Z %*% fit[["DyadCoef"]]
      z_map <- apply(fit[["MixedMembership1"]], 2, which.max) 
      w_map <- z_map
      if(bipartite){
        w_map <- apply(fit[["MixedMembership2"]], 2, which.max)
      }
      hessTheta_list <- lapply(1,
                               function(i, eta, z, w, B, ind){
                                 offset_bm <- B[cbind(z_map[ind[,1]], w_map[ind[,2]])]
                                 pred_edges <- plogis(offset_bm + eta)
                                 return(vcovGamma_ext(Z, pred_edges, c(ctrl$var_gamma)))
                               }, eta = edge_eta, z=z_samples, w=w_samples, B=fit[["BlockModel"]], ind = as.matrix(mfd[,c("(sid)","(rid)")]))
      fit$vcov_dyad <- as.matrix(hessTheta_list[[1]])
      colnames(fit$vcov_dyad) <- rownames(fit$vcov_dyad) <- names(fit[["DyadCoef"]])
    }
    
    
    
    if(ctrl$verbose){
      cat("done.\n")
    }
    
  }#end Hessian portion
  
  
  #Include used data
  attr(mfm1, "terms") <- NULL
  fit$monadic.data <- list(mfm1) #nodes are in nid1,nide2 naming conventions; time in tid convention
  attr(mfd, "terms") <- NULL
  fit$dyadic.data <- mfd
  fit$Y <- Y
  #
  if(bipartite){
    attr(mfm2, "terms") <- NULL
    fit$monadic.data[[2]] <- mfm2
  }
  
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
                    formula.dyad = formula.dyad,
                    formula.monad = formula.monad)
  
  ## Include used seed
  fit$seed <- ctrl$seed
  
  ## Include original call
  fit$call <- match.call()
  
  fit$bipartite <- bipartite
  
  ##Assign class for methods
  if(fit$bipartite){
    class(fit) <- c("mmsbmB", "mmsbm")
  }else{
    class(fit) <- "mmsbm"
  }
  return(fit)
}

