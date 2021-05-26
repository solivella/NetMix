#to go with plot.mmsbm
## Adapted from `fields`` package under GPL
#' @rdname auxfuns
.bar.legend <- function(colPalette, range){
  col <- rgb(colPalette(seq(0,1, length.out = 100)), maxColorValue = 255)
  zlim <- c(0, 1)
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar), add = TRUE)
  stick <- TRUE
  nlevel <- length(col)
  midpoints <- seq(zlim[1], zlim[2], length.out=nlevel)
  delta <- (midpoints[2] - midpoints[1])/2
  breaks <- c(midpoints[1] - delta, midpoints + delta)
  legend.mar <- 3.0
  char.size <- par()$cin[1]/par()$din[1]
  offset <- char.size * par()$mar[4]
  legend.width <- char.size * 1.2
  legend.mar <- legend.mar * char.size
  smallplot <- opar$plt
  smallplot[2] <- 1 - legend.mar
  smallplot[1] <- smallplot[2] - legend.width
  pr <- (smallplot[4] - smallplot[3]) * ((1 - 0.9)/2)
  smallplot[4] <- smallplot[4] - pr
  smallplot[3] <- smallplot[3] + pr
  bigplot <- opar$plt
  bigplot[2] <- min(bigplot[2], smallplot[1] - offset)
  dp <- smallplot[2] - smallplot[1]
  smallplot[1] <- min(bigplot[2] + offset, smallplot[1])
  smallplot[2] <- smallplot[1] + dp
  
  if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
    stop("plot region too small to add legend\n")
  }
  ix <- 1:2
  iy <- breaks
  nBreaks <- length(breaks)
  midpoints <- (breaks[1:(nBreaks - 1)] + breaks[2:nBreaks])/2
  iz <- matrix(midpoints, nrow = 1, ncol = length(midpoints))
  
  par(new = TRUE, pty = "m", plt = smallplot, err = -1)
  graphics::image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
                  ylab = "", col = col, breaks = breaks)
  axis.args <- c(list(side =  4, at = seq(0,1, length.out = 10),
                      labels = round(seq(range[1], range[2], length.out=10),2), 
                      cex.axis=0.75,
                      mgp = c(3, 1, 0), las = 2))
  
  do.call(graphics::axis, axis.args)
  graphics::box()
  legend.args <- list(text = "Edge\nProbability", side = 3, line = 0.5, cex = 1, adj = 0)
  do.call(mtext, legend.args)
}

#find best k: object is object of mmsbm or mmsbmB class; 
#for former just need k1 to have vector of integer values to try
#for latter need k1 and k2 to have vectors of integer values to try
set.seed(123)
sample_split_index<-sample(1:node2,node2,replace=FALSE) #sample split based on Family 2
train.data.dyad<- net3$dyad.data[which(net3$dyad.data$node2%in%sample_split_index[1:round(node2*0.8)]),]
test.data.dyad<- net3$dyad.data[which(net3$dyad.data$node2%in%sample_split_index[(round(node2*0.8)+1):length(sample_split_index)]),]
train.data.monad1<- net3$monad1.data
train.data.monad2<- net3$monad2.data[which(net3$monad2.data$node2%in%sample_split_index[1:round(node2*0.8)]),]
test.data.monad2<- net3$monad2.data[which(net3$monad2.data$node2%in%sample_split_index[81:length(sample_split_index)]),]


formula.dyad = Y ~ V1
formula.monad1 = ~ V2
formula.monad2 = ~ V2
senderID = "node1"
receiverID = "node2"
nodeID1 = "node1"
nodeID2 = "node2"
timeID = "time"
data.dyad = net3$dyad.data
data.monad1 = net3$monad1.data
data.monad2 = net3$monad2.data
n.blocks1 = net3$BLK1
n.blocks2 = net3$BLK2
n.hmmstates = net3$STATE
directed = FALSE
nodes2 = net3$NODE2
npred2 = mpred2+1
mmsbm.control = list(mu_b = c(5,-5)*-1,var_b = c(1, 1),spectral = TRUE,verbose=TRUE,em_iter = 100,conv_tol = 1e-5,bipartite = TRUE,batch_size1 = .80, batch_size2 = .80
)
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
                  nodes2 = NULL,
                  npred2 = NULL,
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
                         n.hmmstates = n.hmmstates, directed = FALSE, nodes2 = nrow(train.data.monad2), npred2 = npred2,
                         mmsbm.control = list(mu_b = c(5,-5)*-1,
                                              var_b = c(1, 1),
                                              spectral = TRUE,
                                              verbose=TRUE,
                                              em_iter = 100,
                                              conv_tol = 1e-5,
                                              bipartite = TRUE,
                                              batch_size1 = .80, batch_size2 = .80
                         ))
      ## Testing
      if (!is.null(test.data.monad1)){new.data.monad1<-test.data.monad1
      }else{new.data.monad1<-train.data.monad1}
      if (!is.null(test.data.monad2)){new.data.monad2<-test.data.monad2
      }else{new.data.monad2<-train.data.monad2}
      pred[[i]]<-predict.mmsbmB(object=train_models[[i]], new.data.dyad = test.data.dyad,
                     new.data.monad1  = new.data.monad1, 
                     new.data.monad2  = new.data.monad2,
                     parametric_mm = ifelse(!is.null(test.data.monad1)|!is.null(test.data.monad2),TRUE,FALSE),
                     forecast = FALSE,
                     type = "response")
      
      # Syntax (response, predictor):
      model_auc[i]<-auc(test.data.dyad$Y, pred[[i]])
    }
    
    return(list(bestk=ktry[which.max(model_auc),],k=ktry, train_models=train_models, pred_response=pred, model_auc=model_auc  ))
}

mmsbmB(formula.dyad = Y ~ V1,#1,#V1, #+ V2,
       formula.monad1 = ~ V2, 
       formula.monad2 = ~ V2,
       senderID = "node1",
       receiverID = "node2",
       nodeID1 = "node1",
       nodeID2 = "node2",
       timeID = "time",
       data.dyad = net3$dyad.data,
       data.monad1 = net3$monad1.data,
       data.monad2 = net3$monad2.data,
       n.blocks1 = net3$BLK1,
       n.blocks2 = net3$BLK2,
       n.hmmstates = net3$STATE,
       directed = FALSE,
       #pass true B, beta, gamma
       #true.B_t = t(net3$B),
       #true.beta1 = true.beta1,
       #true.beta2 = true.beta2,
       #true.gamma = net3$gamma_mat,
       #true.theta = theta,#unlist(theta_v),#unlist(net1$theta_full),
       #pass true phis
       #true.phi1 = net3$pi1_vecs,
       #true.phi2 = net3$pi2_vecs,
       nodes2 = net3$NODE2,
       npred2 = mpred2+1,#+1 for intercept
       mmsbm.control = list(mu_b = c(5,-5)*-1,
                            var_b = c(1, 1),
                            spectral = TRUE,
                            verbose=TRUE,
                            em_iter = 100,
                            conv_tol = 1e-5,
                            bipartite = TRUE,
                            batch_size1 = .80, batch_size2 = .80
       ))

#to go with predict.mmsbm
#' @rdname auxfuns
.pi.hat <- function(X, beta){
  if(dim(beta)[3]==1){
    mu <- exp(X %*% beta[,,1])
    pi.states <- list(t(mu))
  } else {
    pi.states <- lapply(1:dim(beta)[3], function(m){
      mu <- exp(X %*% beta[,,m])
      return(t(mu))
    })
  }
  return(pi.states)
}
#to go with predict.mmsbm
#' @rdname auxfuns
.mpower <- function(mat, p){
  orig <- mat
  while(p > 1){
    mat <- mat %*% orig
    p <- p - 1 
  }
  return(mat)
}
#to go with predict.mmsbm
#' @rdname auxfuns
.e.pi <- function(pi_l, kappa){
  if(is.null(dim(kappa))){
    return(pi_l[[1]])
  } else {
    n_states <- nrow(kappa)
    pi.states <- lapply(1:n_states,
                        function(m){
                          pi_l[[m]] * rep(kappa[m,], each=nrow(pi_l[[m]])) 
                        })
    return(Reduce("+", pi.states))
  }
}

### Bipartite degree calculations
# Raw.
degree.bipartite <- function(bipartite.net){
  degree.n <- rowSums(as.matrix(bipartite.net)) # degree for the first node set.
  degree.m <- colSums(as.matrix(bipartite.net)) # degree for the second node set.
  return(list=c(degree.n=degree.n,degree.m=degree.m)) # put them both together.
}

# Standardized.
degree.s.bipartite <- function(bipartite.net){
  degree.s.n <- (rowSums(as.matrix(bipartite.net)))/dim(as.matrix(bipartite.net))[2] # degree for the first node set.
  degree.s.m <- (colSums(as.matrix(bipartite.net)))/dim(as.matrix(bipartite.net))[1] # degree for the second node set.
  return(c(degree.s.n,degree.s.m)) # put them both together.
}


