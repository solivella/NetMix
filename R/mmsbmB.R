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
#'        \item{spectral}{Boolean. Type of initialization algorithm for mixed-membership vectors in static case. If \code{TRUE} (default),
#'                    use spectral clustering with degree correction; otherwise, use kmeans algorithm}
#'        \item{init.dyn.gibbs}{Boolean. Should a collapsed Gibbs sampler of non-regression mmsbm be used to initialize
#'                    each time period? Setting to \code{TRUE} will be result in faster estimation that is very sensitive to
#'                    choice of alpha (see below)}            
#'        \item{alpha}{Numeric positive value. Concentration parameter for collapsed Gibbs sampler to find initial
#'                     mixed-membership values in dynamic case when \code{init.dyn.gibbs=TRUE}. Defaults to 0.5}            
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
                  missing="indicator method",
                  true.B_t = NULL,
                  true.beta1 = NULL,
                  true.beta2 = NULL,
                  true.gamma = NULL,
                  true.phi1 = NULL,
                  true.phi2 = NULL,
                  true.theta = NULL,
                  true.kappa = NULL,
                  nodes2 = NULL,
                  npred2 = NULL,
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
               nodes2 = NULL,
               states = n.hmmstates,
               times = 1,
               directed = TRUE,
               seed = sample(500,1),
               bipartite = TRUE,
               phi1_init_t = NULL,
               phi2_init_t = NULL,
               kappa_init_t = NULL,
               b_init_t = NULL,
               beta1_init = NULL,
               beta2_init = NULL,
               gamma_init = NULL,
               theta_init = NULL,
               spectral = TRUE,
               alpha = 0.5,
               missing = "indicator method",
               em_iter = 50,
               hessian=TRUE,
               se_sim=10,
               dyad_vcov_samp=1000,
               opt_iter = 10e3,
               
               init.dyn.gibbs = if (n.hmmstates > 1) TRUE else FALSE,
               forget_rate = 0.75,
               delay = 1.0,
               batch_size1 = 0.25,
               batch_size2 = 0.25,
               
               
               mu_b = c(1.0, 0.0),
               var_b = c(1.0, 1.0),
               mu_beta1 = 0.0,
               var_beta1 = 1.0,
               mu_beta2 = 0.0,
               var_beta2 = 1.0,
               mu_gamma = 0.0,
               var_gamma = 1.0,
               eta = 1.3,
               permute = TRUE,
               threads = 1,
               conv_tol = 1e-2,
               verbose = FALSE)
  
  # #pass truth if available: gamma
  if(!is.null(true.gamma)){
    ctrl$gamma_init<-true.gamma
  }
  # #pass truth if available: beta
  if(!is.null(true.beta1)){
    ctrl$beta1_init<-true.beta1
  }
  if(!is.null(true.beta2)){
    ctrl$beta2_init<-true.beta2
  }
  # #pass truth if available: B
  if(!is.null(true.B_t)){
    ctrl$b_init_t<-true.B_t
  }
  
  # #pass truth if available: phi
  if(!is.null(true.phi1)){
    ctrl$phi1_init_t <- t(true.phi1)
  }
  if(!is.null(true.phi2)){
    ctrl$phi2_init_t <- t(true.phi2)
  }
  
  directed<-ctrl$directed ## INSERTED HACK 
  
  ctrl[names(mmsbm.control)] <- mmsbm.control
  if(ctrl$bipartite == FALSE){
    ctrl$nodes2 <- 0
    ctrl$npred2 <- 0
  }else{
    ctrl$nodes2 <- nodes2
    ctrl$npred2 <- npred2
  }
  if(!is.null(ctrl$seed)){ set.seed(ctrl$seed)}
  mu_b <- var_b <- array(NA, c(n.blocks2, n.blocks1))
  diag(mu_b) <- ctrl[["mu_b"]][1]
  mu_b[upper.tri(mu_b)|lower.tri(mu_b)] <- ctrl[["mu_b"]][2]
  diag(var_b) <- ctrl[["var_b"]][1]
  var_b[upper.tri(var_b)|lower.tri(var_b)] <- ctrl[["var_b"]][2]
  
  if(ctrl$verbose){
    cat("Pre-processing data...\n")
  }
  
  stopifnot(class(formula.dyad) == "formula",
            class(formula.monad1) == "formula",
            class(formula.monad2) == "formula",
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
  
  ## Add time variable if only one period
  if(is.null(timeID)){
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
  if(missing=="indicator method"){
    if(ctrl$bipartite){
      .missing(formula.dyad, formula.monad1, formula.monad2, 
               data.dyad, data.monad1, data.monad2, missing=="indicator method")
    }else{
      .missing(formula.dyad, formula.monad1, formula.monad2=NULL, 
               data.dyad, data.monad1, data.monad2=NULL, missing=="indicator method")
    }
  }
  if(missing=="listwise deletion"){
    if(ctrl$bipartite){
      .missing(formula.dyad, formula.monad1, formula.monad2, 
               data.dyad, data.monad1, data.monad2, missing=="listwise deletion")
    }else{
      .missing(formula.dyad, formula.monad1, formula.monad2=NULL, 
               data.dyad, data.monad1, data.monad2=NULL, missing=="listwise deletion")
    }
  }
  
  mfd <- do.call(model.frame, list(formula = formula.dyad,
                                   data = data.dyad,
                                   drop.unused.levels = TRUE,
                                   tid = as.name(timeID),
                                   sid = as.name(senderID),
                                   rid = as.name(receiverID))) #ok for differential names
  mfd[,"(sid)"] <- as.character(mfd[,"(sid)"])
  mfd[,"(rid)"] <- as.character(mfd[,"(rid)"])
  #dyadic_order <- with(mfd, order(`(tid)`, `(sid)`, `(rid)`))
  #mfd <- mfd[dyadic_order, ]
  
  ut <- unique(mfd[["(tid)"]])
  periods <- length(ut)
  if(periods > 1){
    ctrl$times <- periods
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
  
  Y <- model.response(mfd)
  X1 <- scale(model.matrix(terms(mfm1), mfm1))
  X1_mean <- attr(X1, "scaled:center")
  X1_sd <- attr(X1, "scaled:scale")
  if(any(X1_sd==0)){
    constx <- which(X1_sd==0)
    X1[,constx] <- 1
    X1_sd[constx]<-1##ADDED to remove X1_sd elements that ==0 and set 1
  }
  n_monad1_pred <- ncol(X1)
  X2 <- scale(model.matrix(terms(mfm2), mfm2))
  X2_mean <- attr(X2, "scaled:center")
  X2_sd <- attr(X2, "scaled:scale")
  if(any(X2_sd==0)){
    constx <- which(X2_sd==0)
    X2[,constx] <- 1
    X2_sd[constx]<-1##ADDED to remove X2_sd elements that ==0 and set 1
  }
  n_monad2_pred <- ncol(X2)
  
  #sub_mfd <- mfd[,-match(all.vars(update(formula.dyad, .~0)) , names(mfd))] #all.vars(update(formula.dyad, .~0)) is name of outcome Y
  Z <- scale(model.matrix(terms(mfd), mfd))
  Z_mean <- attr(Z, "scaled:center")
  Z_sd <- attr(Z, "scaled:scale")
  if(any(Z_sd==0)){
    constz <- which(Z_sd==0)
    Z <- Z[,-constz, drop = FALSE]
  }
  n_dyad_pred <- ncol(Z)
  
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
  
  ## Modify prior means and variances to match transformed model matrix
  ctrl$mu_gamma <- .transf_muvar(ctrl$mu_gamma, FALSE, FALSE, Z)
  ctrl$var_gamma <- .transf_muvar(ctrl$var_gamma, TRUE, FALSE, Z, devs=Z_sd[-which(Z_sd==0)])
  ctrl$mu_beta1 <- .transf_muvar(ctrl$mu_beta1, FALSE, TRUE, X1, n.blocks1, n.hmmstates) 
  ctrl$mu_beta2 <- .transf_muvar(ctrl$mu_beta2, FALSE, TRUE, X2, n.blocks2, n.hmmstates) 
  ctrl$var_beta1 <- .transf_muvar(ctrl$var_beta1, TRUE, TRUE, X1, n.blocks1, n.hmmstates, c(1,X1_sd[-1]))
  ctrl$var_beta1[1,,]<-0.5 #### just for the intercepts, you could define a standard deviation of 0.5
  ctrl$var_beta2 <- .transf_muvar(ctrl$var_beta2, TRUE, TRUE, X2, n.blocks2, n.hmmstates, c(1,X2_sd[-1]))
  ctrl$var_beta2[1,,]<-0.5 #### just for the intercepts, you could define a standard deviation of 0.5
  
  ## Translate batch size to number of nodes
  ctrl$batch_size1 = max(1, floor(ctrl$batch_size1 * sum(nodes_pp1)))
  cat("batch_size1 is:",ctrl$batch_size1,"\n")
  ctrl$batch_size2 = max(1, floor(ctrl$batch_size2 * sum(nodes_pp2)))
  cat("batch_size2 is:",ctrl$batch_size2,"\n")
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
  if(ctrl$bipartite){
    soc_mats <- .createSocioB(dyads,all.nodes1,all.nodes2, ctrl$directed)
  }else{
    soc_mats <- .createSocioB(dyads,all.nodes,all.nodes, ctrl$directed)
  }
  
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
      state_init <- fitted(kmeans(dyad_time,
                                  n.hmmstates,
                                  nstart = 15), "classes")
      # state_init <- cluster::clara(dyad_time,
      #                              n.hmmstates,
      #                              samples = 15,
      #                              sampsize = min(nrow(dyad_time), 100 + 2 * n.hmmstates),
      #                              rngR = TRUE,
      #                              correct.d = TRUE)$clustering
      kappa_internal <- model.matrix(~ as.factor(state_init) - 1)
      kappa_internal <- .transf(kappa_internal)
      ctrl$kappa_init_t <- t(kappa_internal)
    } else {
      ctrl$kappa_init_t <- t(matrix(1, nrow = periods))
      state_init <- apply(ctrl$kappa_init_t, 2, which.max)
    }
  } else {
    state_init <- apply(ctrl$kappa_init_t, 2, which.max)
  } 
  if(n.hmmstates==1){
    names(state_init) = 1
  }
  
  
  #null Phi 1 and 2
  ## Adjusted this for when edges are never formed for a given node/year
  if(ctrl$bipartite){
    if((is.null(ctrl$phi1_init_t)|is.null(ctrl$phi2_init_t)) & (periods == 1)){
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
      ctrl$phi1_init_t <- do.call(cbind, lapply(phi_init_temp, `[[`, 1))#1st matrix of each element of big list
      ctrl$phi2_init_t <- do.call(cbind, lapply(phi_init_temp, `[[`, 2))#2nd matrix of each element of big list
    } 
  }else{
    if(is.null(ctrl$phi1_init_t)){
      if(is.null(ctrl$phi1_init_t) & (periods == 1)){
        phi_init_temp <- lapply(soc_mats, function(mat){
          if(ctrl$directed){
            D_o <- 1/sqrt(rowSums(mat) + 1)
            D_i <- 1/sqrt(colSums(mat) + 1)
            C_o <- t(D_o * mat)
            C_i <- t(D_i * t(mat))
            U <- t(C_o * D_i) %*% C_o +
              t(C_i * D_o) %*% C_i
          } else {
            D <- 1/sqrt(rowSums(mat) + 1)
            U <- t(D * mat) * D
          }
          if(ctrl$spectral) {#ignore spectral portion here
            n_elem <- n.blocks1 + 1
            res <- RSpectra::eigs_sym(U, n_elem) 
            eta <- res$vectors[,1:n_elem] %*% diag(res$values[1:n_elem]) #n_elem largest eigeenvalues
            target <- eta[,2:n_elem] / (eta[,1] + 1e-8)
            sig <- 1 - res$values[n_elem] / (res$values[n.blocks1]) #family1: n.blocks1
            sig <- ifelse(is.finite(sig), sig, 0)
            if(sig > 0.1){
              target <- target[,1:(n_elem - 2)]
            }
          } else {
            target <- U
          } #Kmeans -- detecting communities in family 1
          init_c <- sample(1:nrow(target), n.blocks1, replace = FALSE) #family1: n.blocks1
          clust_internal <- fitted(kmeans(target,
                                          n.blocks1, #family1: n.blocks1
                                          centers = target[init_c,],
                                          #algorithm = "Lloyd",
                                          nstart = 10), "classes")
          
          phi_internal <- model.matrix(~ as.factor(clust_internal) - 1)
          phi_internal <- .transf(phi_internal)
          rownames(phi_internal) <- rownames(mat)
          colnames(phi_internal) <- 1:n.blocks1 #family1: n.blocks1
          t(phi_internal)
        })
        #phi_init_orig <- do.call(cbind, phi_init_temp)
        #ctrl$phi_list <- phi_init_temp
        ctrl$phi1_init_t <- do.call(cbind, phi_init_temp) #family1: return ctrl$phi1_init_t
      }
      else if (is.null(ctrl$phi1_init_t) & (periods > 1)){
        temp_res <- vector("list", periods)
        for(i in 1:periods){
          if (ctrl$init.dyn.gibbs) {
            n_prior <- (dyads_pp[i] - nodes_pp[i]) * .05 #need to shift nodes_pp by length periods for family 2; no need for family 1
            a <- plogis(mu_b) * n_prior
            b <- n_prior - a
            lda_beta_prior <- lapply(list(b,a),
                                     function(prior){
                                       mat <- matrix(prior[2], n.blocks1, n.blocks1) # family1: n.blocks1/n.blocks1
                                       diag(mat) <- prior[1]
                                       return(mat)
                                     })
            ret <- lda::mmsb.collapsed.gibbs.sampler(network = soc_mats[[i]],
                                                     K = n.blocks1, # family1: n.blocks1
                                                     num.iterations = 100L,
                                                     burnin = 50L,
                                                     alpha = ctrl$alpha, #alpha can be 2x
                                                     beta.prior = lda_beta_prior)
            BlockModel <- with(ret, blocks.pos / (blocks.pos + blocks.neg + 1))
            MixedMembership <- prop.table(ret$document_expects, 2)
            colnames(MixedMembership) <- colnames(soc_mats[[i]])
            MixedMembership <- MixedMembership[,colnames(MixedMembership) %in% ntid1] # family1: ntid1
            temp_res[[i]] <- list(BlockModel = BlockModel,
                                  MixedMembership = MixedMembership)
          } 
          else {
            mfd_list <- split(mfd, mfd[,c("(tid)")])
            mfm_list <- split(mfm1, mfm1[,c("(tid)")]) #family1: mfm1
            temp_res[[i]] <- mmsbm(update(formula.dyad, .~1),
                                   formula.monad1 = ~ 1, #family1: formula.monad1
                                   senderID = "(sid)",
                                   receiverID = "(rid)",
                                   nodeID = "(nid)",
                                   timeID = "(tid)",
                                   data.dyad = mfd_list[[i]],
                                   data.monad1 = mfm_list[[i]][,c("(nid)", "(tid)")], #family1: data.monad1
                                   n.blocks1 = n.blocks1,
                                   n.blocks2 = n.blocks2,
                                   n.hmmstates = 1,
                                   directed = ctrl$directed,
                                   missing = missing,
                                   mmsbm.control = list(em_iter = ctrl$em_iter,
                                                        seed = ctrl$seed,
                                                        mu_b = ctrl$mu_b,
                                                        var_b = ctrl$var_b,
                                                        spectral = ctrl$spectral,
                                                        conv_tol = ctrl$conv_tol,
                                                        threads = ctrl$threads,
                                                        verbose = FALSE))
            temp_res[[i]]$MixedMembership <- temp_res[[i]]$MixedMembership[,colnames(temp_res[[i]]$MixedMembership) %in% ntid1] #family1: ntid1
          }
        }
        
        temp_res <- lapply(split(temp_res, state_init),
                           function(mods){
                             target <- t(mods[[1]]$MixedMembership)
                             rownames(target) <- sapply(strsplit(rownames(target), "@", fixed = TRUE, useBytes = TRUE), function(x)x[1])
                             res <- lapply(mods,
                                           function(mod, target_mat = target){
                                             split_names <- strsplit(colnames(mod$MixedMembership), "@", fixed = TRUE, useBytes = TRUE)
                                             mod_names <-  sapply(split_names, function(x)x[1])
                                             mod_time <- split_names[[1]][2]
                                             shared_nodes <- intersect(mod_names,
                                                                       rownames(target_mat))
                                             shared_nodes_mod <- paste(shared_nodes, mod_time, sep="@")
                                             cost_mat <- mod$MixedMembership[,shared_nodes_mod] %*% target_mat[shared_nodes,]
                                             perm <- clue::solve_LSAP(t(cost_mat), TRUE)
                                             mod$MixedMembership <- mod$MixedMembership[perm,]
                                             mod$BlockModel <- mod$BlockModel[perm, perm]
                                             return(mod)
                                           })
                             return(res)
                           })
        block_models <- lapply(temp_res, 
                               function(mods){
                                 Reduce("+", Map(function(x)x$BlockModel, mods)) / length(mods)
                               }) 
        perms_temp <- .findPerm(block_models, use.perms = ctrl$permute)
        phis_temp <- Map(function(x)x$MixedMembership, unlist(temp_res, recursive = FALSE))
        ctrl$phi1_init_t <- do.call(cbind,mapply(function(phi,perm){perm %*% phi},
                                                 phis_temp, perms_temp[state_init], SIMPLIFY = FALSE)) #this is now nblock1 x (n1 * nstate)
        ctrl$phi1_init_t<-ctrl$phi1_init_t[,match(ntid1,colnames(ctrl$phi1_init_t))] #order cols by ntid1
        rownames(ctrl$phi1_init_t) <- 1:n.blocks1 #family1: n.blocks1
      }
      
    }#end of nullPhi family1
    #null Phi 2  
    ctrl$phi2_init_t  <- ctrl$phi1_init_t #if not bipartite, then just let nullPhi 2 be nullPhi 1
  }
  
  
  #null Gamma
  if(is.null(ctrl$gamma_init)){
    if(ncol(Z) > 0){
      nsamp <- ifelse(nrow(Z) > 100e3, 100e3, nrow(Z))
      ind_dyad <- sample(1:nrow(Z), nsamp)
      ctrl$gamma_init <- coef(glm(Y[ind_dyad] ~ Z[ind_dyad,] - 1, family = binomial()))
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
    ctrl$b_init_t<-qlogis(approxB(Y, nt_id, ctrl$phi1_init_t, ctrl$phi2_init_t))
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
      phi1<-t(ctrl$phi1_init_t) #t(blk x node) -> node x blk
      phi2<-t(ctrl$phi2_init_t) 
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
    phi_state <- split.data.frame(t(ctrl$phi1_init_t), state_init[t_id_n1 + 1])
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
    phi_state <- split.data.frame(t(ctrl$phi2_init_t), state_init[t_id_n2 + 1])
    ctrl$beta2_init <- mapply(function(X_sub, phi_sub){
      phi_temp <- .transf(phi_sub)
      lm.fit(X_sub, log(phi_temp))$coefficients},
      X_state, phi_state,SIMPLIFY = "array")
  }
  if(anyNA(ctrl$beta2_init)){
    stop("Nearly singular design matrix in family 2; check monadic predictors.")
  }
  
  # null Theta
  if(!is.null(true.theta)){
    ctrl$theta_init<-true.theta
  }else{
    if(ctrl$bipartite){
      ctrl$theta_init<-array(rep(0, n.blocks1*n.blocks2*length(Y)), c(n.blocks2,n.blocks1,length(Y))) #dim theta: h,g,dyad
    }else{
      ctrl$theta_init<-rep(0, n.blocks1*n.blocks1*length(Y), c(n.blocks1,n.blocks1,length(Y)))
    }
  }
  
  ## Create randomizer for order of updatePhis
  ctrl$phi_order<-rbinom(nrow(Z)[1],1,0.5) #ndyad
  #test print
  #print(ctrl$b_init_t)
  #print(ctrl$beta2_init)
  
  ## Estimate model
  fit <- mmsbm_fit(t(Z),#N_DYAD_PRED x N_DYAD=t(Z)
                   t(X1),#2 x 50; (mpred1+1) x n1 = t(X1[,-1])
                   t(X2),#2 x 50; (mpred2+1) x n2 = t(X2[,-1])
                   Y, #length n1*n2, N_DYAD=Y
                   t_id_d,#length N_DYAD*TIME filled with 0s = t_id_d
                   t_id_n1,#length n1*time filled with 0s = t_id_n1
                   t_id_n2,#
                   nodes_pp, #stacked first family 1 nodes for each period, then family 2 nodes for each period
                   nodes_pp1,# length TIME, number of nodes per period in family 1
                   nodes_pp2,# length TIME, number of nodes per period
                   nt_id, #cbind family 1 and family 2
                   mu_b, #should be constant matrix; note should be passing it blk2 x blk1
                   var_b, #should be constant matrix; note should be passing it blk2 x blk1
                   ctrl$mu_beta1, #cube
                   ctrl$var_beta1, #cube
                   ctrl$mu_beta2,
                   ctrl$var_beta2,
                   ctrl$mu_gamma,
                   ctrl$var_gamma,#constant vector
                   ctrl$phi1_init_t, #for phi_init1
                   ctrl$phi2_init_t, #for phi_init2
                   ctrl$kappa_init_t, ##
                   ctrl$b_init_t,
                   ctrl$beta1_init, #for beta_init1
                   ctrl$beta2_init, #for beta_init2
                   ctrl$gamma_init,
                   ctrl$theta_init,
                   ctrl$phi_order,
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
                    formula.monad2 = formula.monad2,
                    nodes2 = nodes2,
                    npred2 = npred2
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

