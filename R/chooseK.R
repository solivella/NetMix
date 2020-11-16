#' Find best values for K latent groups based on auc
#'
#' The function produces a pair of K values that produce the largest auc value among a set of K options 
#' introduced by the user.
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
#' @param k1 Vector of integer values. How many latent groups for Family 1 should be considered in the search?
#' @param k2 Vector of integer values. How many latent groups for Family 2 should be considered in the search?
#' @param n.hmmstates Integer value. How many hidden Markov state should be used in the HMM? Defaults 
#'    to 1 (i.e. no HMM).  
#' @param directed Boolean. Is the network directed? Defaults to \code{TRUE}.
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
#'  @return List with named components:
#'     \describe{
#'       \item{bestk}{Row of best combination of k1 and k2.}
#'       \item{ktry}{matrix of all combinations of k1 and k2 values searched.}
#'       \item{train_models}{List of returned \code{mmsbmB} trained models of length all combinations of k1 and k2 values searched.
#'       Index corresponds to rows of \code{ktry}.}
#'       \item{pred_response}{List of returned predicted responses for each of the models run on testing data. Index corresponds to rows of \code{ktry}.}
#'       \item{model_auc}{Vector of auc values for each combination of k1 and k2. Index corresponds to rows of \code{ktry}.}
#'     }
#' @author Kosuke Imai (imai@@harvard.edu), Adeline Lo (aylo@@wisc.edu), Santiago Olivella (olivella@@unc.edu), Tyler Pratt (tyler.pratt@@yale.edu) 

chooseK<-function(formula.dyad,
                  formula.monad1=~1,
                  formula.monad2=~1,
                  senderID, 
                  receiverID,
                  nodeID1 = NULL,
                  nodeID2 = NULL,
                  timeID = NULL,
                  train.data.dyad, test.data.dyad,
                  train.data.monad1 = NULL, test.data.monad1 = NULL,
                  train.data.monad2 = NULL, test.data.monad2 = NULL,
                  k1=NULL,k2=NULL,
                  n.hmmstates,
                  directed = TRUE,
                  #nodes2 = NULL,
                  #npred2 = NULL,
                  mmsbm.control = list(),seed=123){
  set.seed(seed)
  require(pROC)
  if(is.null(k1)|is.null(k2)|any(2>k1)|any(2>k2)){stop("Error: mmsbmB object needs vectors of k1 AND k2 (greater than 1) values to try!")}
  ktry<-expand.grid(k1,k2) #dataframe with rows of combinations of k to try
  colnames(ktry)<-c("k1","k2")
  train_models<-pred<-vector("list",nrow(ktry))
  model_auc<-rep(NA,nrow(ktry))
  train_nodes2<-ifelse(is.null(train.data.monad2),nodes2,nrow(train.data.monad2))
  for(i in 1:nrow(ktry)){
    ## Training
    train_models[[i]] <- mmsbmB(formula.dyad = formula.dyad, formula.monad1 = formula.monad1,formula.monad2 = formula.monad2,
                                senderID = senderID, receiverID = receiverID,nodeID1 = nodeID1, nodeID2 = nodeID2, timeID = timeID,
                                data.dyad = train.data.dyad, data.monad1 = train.data.monad1, data.monad2 = train.data.monad2, #training data
                                n.blocks1 = ktry[i,1], n.blocks2 = ktry[i,2], #k to try
                                n.hmmstates = n.hmmstates, directed = FALSE, #nodes2 = nrow(train.data.monad2), npred2 = npred2,
                                mmsbm.control = mmsbm.control)
    ## Testing
    if (!is.null(test.data.monad1)){new.data.monad1<-test.data.monad1
    }else{new.data.monad1<-train.data.monad1}
    if (!is.null(test.data.monad2)){new.data.monad2<-test.data.monad2
    }else{new.data.monad2<-train.data.monad2}
    pred[[i]]<-predict.mmsbmB(object=train_models[[i]], 
                              new.data.dyad = test.data.dyad,
                              new.data.monad1  = new.data.monad1, 
                              new.data.monad2  = new.data.monad2,
                              parametric_mm = TRUE,
                              forecast = FALSE,
                              type = "response")
    
    # Syntax (response, predictor):
    model_auc[i]<-auc(test.data.dyad$Y, pred[[i]])
  }
  
  return(list(bestk=ktry[which.max(model_auc),],k=ktry, train_models=train_models, pred_response=pred, model_auc=model_auc  ))
}