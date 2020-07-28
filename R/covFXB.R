#' Generate estimated monadic covariate effects for estimated mmsbmB model
#'
#' The function estimates the effect of a shift in monadic covariate values (for one family) on the probability of edge formation in the network. 
#'
#' @param fm An object of class \code{mmsbmB}, a result of a call to \code{mmsbmB}. 
#' @param cov Character string identifying the monadic covariate to be shifted.  
#' @param shift Numeric value specifying the desired increase or decrease in the monadic covariate.  The monadic predictor will be shifted by this value for all nodes and time periods.
#' @param max.val An optional numeric value specifying the maximum possible value for the monadic covariate.
#'
#'     
#' @return List with named components:
#'     \describe{
#'       \item{Overall Avg. Effect}{Overall average effect of the covariate shift on the predicted probability of edge formation.}
#'       \item{Avg. Effect by Time}{Vector of average effects of the covariate shift on the predicted probability of edge formation for each time period.}
#'       \item{Avg. Effect by Node}{Vector of average effects of the covariate shift on the predicted probability of edge formation for each node.}
#'       \item{Avg. Effect by Dyad}{Vector of average effects of the covariate shift on the predicted probability of edge formation for each node dyad.}
#'       \item{Avg. Effect Dyad-Time}{Vector of estimated effects of the covariate shift on the predicted probability of edge formation for each node dyad-time unit.}
#'     }
#' 
#' @author Santiago Olivella (olivella@@unc.edu), Adeline Lo (aylo@@wisc.edu), Tyler Pratt (tyler.pratt@@yale.edu), Kosuke Imai (imai@@harvard.edu)
#' 

covFXB <- function(fm, cov1=NULL, cov2=NULL, shift, max.val=FALSE){
  #set up IDs
  if(!is.null(fm$call$senderID)){ sid<- fm$call$senderID }else{cat("Error: no senderID information in model object.")}
  if(!is.null(fm$call$receiverID)){ rid<- fm$call$receiverID }else{cat("Error: no receiverID information in model object.")}
  if(!is.null(fm$call$timeID)){ tid<- fm$call$timeID }else{cat("Error: no timeID information in model object.")}
  
  if(is.null(cov1)&is.null(cov2)){cat("Error: no monadic covariate provided for either family in bipartite network!\n")}
  if(!is.null(cov1)&!is.null(cov2)){cat("Error: only one family's monadic covariate can be provided!\n")}
  predict.ties <- predict.mmsbmB(fm, parametric_mm = TRUE, type="response")  ## check correct predict called
  #predict with shift in Family 1 covariate
  if(!is.null(cov1)){ 
    cov<-cov1
    monadic.data2 <- fm$monadic1.data
    monadic.data2[,cov] <- fm$monadic1.data[,cov] + shift
    if(!isFALSE(max.val)){
      monadic.data2[which(fm$monadic1.data[,cov] == max(fm$monadic1.data[,cov])),cov] <- max.val
    }
    predict.ties2 <- predict.mmsbmB(fm, new.data.monad1=monadic.data2, parametric_mm = TRUE, type="response") #chnage to just predict()?
    FX <- list(mean(predict.ties2 - predict.ties), #avg
               tapply(predict.ties2-predict.ties, fm$dyadic.data[,tid], mean), #time
               sapply(unique(fm$monadic1.data[,sid]), function(x){ #node
                 mean((predict.ties2-predict.ties)[fm$dyadic.data[,sid]==x])}),#always sender
               tapply(predict.ties2-predict.ties, paste(fm$dyadic.data[,sid], fm$dyadic.data[,rid], sep="_"), mean),#dyad
               predict.ties2 - predict.ties#dyad-time
               ,sapply(unique(fm$monadic1.data[,sid]), function(x){ #node
                 mean((predict.ties2)[fm$dyadic.data[,sid]==x])}) #node-newpredicted
               ,sapply(unique(fm$monadic1.data[,sid]), function(x){ #node
                 mean((predict.ties)[fm$dyadic.data[,sid]==x])})) #node-oldpredicted
    names(FX[[3]]) <- unique(fm$monadic1.data[,sid])
    names(FX[[5]]) <- paste(fm$dyadic.data[,sid], fm$dyadic.data[,rid], sep="_")
    names(FX) <- c(paste("Overall Avg. Effect of", cov), paste("Avg. Effect of", cov, "by Time"),
                   paste("Avg. Effect of", cov, "by Node (Family 1)"), paste("Avg. Effect of", cov, "by Dyad"),
                   paste("Effect of", cov, "by Dyad-Time")
                   ,"Node New Avg Predicted"
                   ,"Node Orig Avg Predicted")
  }
  #predict with shift in Family 2 covariate
  if(!is.null(cov2)){ 
    cov<-cov2
    monadic.data2 <- fm$monadic2.data
    monadic.data2[,cov] <- fm$monadic2.data[,cov] + shift
    if(!isFALSE(max.val)){
      monadic.data2[which(fm$monadic2.data[,cov] == max(fm$monadic2.data[,cov])),cov] <- max.val
    }
    predict.ties2 <- predict.mmsbmB(fm, bipartite=TRUE, new.data.monad2=monadic.data2, parametric_mm = TRUE, type="response")
    FX <- list(mean(predict.ties2 - predict.ties), #avg
               tapply(predict.ties2-predict.ties, fm$dyadic.data[,tid], mean), #time
               sapply(unique(fm$monadic2.data[,rid]), function(x){ #node
                 mean((predict.ties2-predict.ties)[fm$dyadic.data[,rid]==x])}),#always receiver
               tapply(predict.ties2-predict.ties, paste(fm$dyadic.data[,sid], fm$dyadic.data[,rid], sep="_"), mean),#dyad
               predict.ties2 - predict.ties#dyad-time
               ,sapply(unique(fm$monadic2.data[,rid]), function(x){ #node
                 mean((predict.ties2)[fm$dyadic.data[,rid]==x])}) #node-newpredicted
               ,sapply(unique(fm$monadic2.data[,rid]), function(x){ #node
                 mean((predict.ties)[fm$dyadic.data[,rid]==x])})) #node-oldpredicted
    names(FX[[3]]) <- unique(fm$monadic2.data[,rid])
    names(FX[[5]]) <- paste(fm$dyadic.data[,sid], fm$dyadic.data[,rid], sep="_")
    names(FX) <- c(paste("Overall Avg. Effect of", cov), paste("Avg. Effect of", cov, "by Time"),
                   paste("Avg. Effect of", cov, "by Node (Family 2)"), paste("Avg. Effect of", cov, "by Dyad"),
                   paste("Effect of", cov, "by Dyad-Time")
                   ,"Node New Avg Predicted"
                   ,"Node Orig Avg Predicted")
  }
  
  return(FX)
}