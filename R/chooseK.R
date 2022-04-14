#' Find best values for K latent groups based on auc
#'
#' The function produces a pair of K values that produce the largest auc value among a set of K options 
#' introduced by the user.
#' @param formula.dyad A \code{formula} object. See \code{\link{mmsbm}}. 
#' @param formula.monad An optional \code{formula} object or list of \code{formula}s. See \code{\link{mmsbm}}.   
#' @param bipartite Boolean. See \code{\link{mmsbm}}. 
#' @param k1 Vector of integer values. How many latent groups for Family 1 (or whole network, if \code{bipartite = FALSE}) should be considered in the search? Defaults to  \code{c(2,3)}.
#' @param k2 Vector of integer values. How many latent groups for Family 2 should be considered in the search? Only required if \code{bipartite = TRUE}.
#' @param senderID Character string. See \code{\link{mmsbm}}. 
#' @param receiverID Character string. See \code{\link{mmsbm}}. 
#' @param nodeID Character string or list of two character strings. See \code{\link{mmsbm}}. 
#' @param timeID Optional character string or list of two character strings. See \code{\link{mmsbm}}. 
#' @param train.data.dyad Data frame. Dyadic data used to train the MMSBM model. See \code{data.dyad} in \code{\link{mmsbm}}. 
#' @param test.data.dyad Data frame. Dyadic data used to test the fit of MMSBM model. See \code{data.dyad} in \code{\link{mmsbm}}. 
#' @param train.data.monad Data frame or list of two data frames used to train the MMSBM model. See \code{data.monad} \code{\link{mmsbm}}.
#' @param test.data.monad Data frame or list of two data frames used to test the fit of MMSBM model. See \code{data.monad} \code{\link{mmsbm}}.
#' @param n.hmmstates Integer value. See \code{\link{mmsbm}}.
#' @param directed Boolean. See \code{\link{mmsbm}}.
#' @param mmsbm.control A named list of optional algorithm control parameters. See \code{\link{mmsbm}}.
#'
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
                  formula.monad=~1,
                  bipartite = FALSE,
                  k1 = c(2,3), k2 = NULL,
                  senderID, 
                  receiverID,
                  nodeID,
                  timeID = NULL,
                  train.data.dyad,
                  test.data.dyad,
                  train.data.monad = NULL,
                  test.data.monad = NULL,
                  n.hmmstates = 1,
                  directed = TRUE,
                  mmsbm.control = list(),
                  seed=123){
  set.seed(seed)
  require(pROC)
  if(bipartite){
    if(is.null(k2) | any(2>k2){
      stop("mmsbm object needs both k1 AND k2 vectors of values greater than 1.")
     }
    if(is.null(train.data.monad)){
      train.data.monad <- list(NULL, NULL)
    }
    if(is.null(test.data.monad)){
      test.data.monad <- list(NULL, NULL)
    }   
  }
  if(any(2>k1){
      stop("mmsbm object needs all values in k1 to be greater than 1.")
  }
  if(bipartite){
      ktry <- expand.grid(k1=k1,k2=k2) #dataframe with rows of combinations of k to try
  } else {
      ktry <- data.frame(k1 = k1)
  }
  train_models<-pred<-vector("list",nrow(ktry))
  model_auc<-rep(NA,nrow(ktry)) 
  for(i in 1:nrow(ktry)){
    ## Training
    train_models[[i]] <- mmsbm(formula.dyad = formula.dyad, formula.monad = formula.monad, 
                                senderID = senderID, receiverID = receiverID, nodeID = nodeID, timeID = timeID,
                                data.dyad = train.data.dyad, data.monad = train.data.monad, #training data
                                n.blocks = ktry[i, ], #k to try
                                n.hmmstates = n.hmmstates, directed = directed, 
                                mmsbm.control = mmsbm.control)
    ## Testing
    if(bipartite){
        pred[[i]]<-predict.mmsbmB(object=train_models[[i]], 
                                 new.data.dyad = test.data.dyad,
                                 new.data.monad1  = test.data.monad[[1]],
                                 new.data.monad2  = test.data.monad[[2]], 
                                 parametric_mm = TRUE,
                                 type = "response")
      
    } else {
        pred[[i]]<-predict.mmsbm(object=train_models[[i]], 
                                 new.data.dyad = test.data.dyad,
                                 new.data.monad  = test.data.monad, 
                                 parametric_mm = TRUE,
                                 type = "response")
    }
    
    # Syntax (response, predictor):
    model_auc[i]<-auc(test.data.dyad$Y, pred[[i]])
  }
  
  return(list(bestk=ktry[which.max(model_auc),],k=ktry, train_models=train_models, pred_response=pred, model_auc=model_auc ))
}
