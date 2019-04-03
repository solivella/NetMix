### Function for simulating networks ###

NetSim2 <- function(BLK = 4,
                    NODE = 20,
                    STATE = 1,
                    TIME = 1,
                    DIRECTED = TRUE,
                    N_PRED = c(2, 2),
                    B = NULL,
                    sVec = NULL,
                    beta_arr = NULL,
                    gamma_vec = NULL,
                    X = NULL,
                    Z = NULL){
  require(MCMCpack)
  
  stopifnot(TIME >= STATE)
  
  if(is.null(sVec)){
    s1v <- rnorm(STATE)
    s1 <- which.max(rmultinom(n=1, size=1, exp(s1v)/sum(exp(s1v))))
    if(STATE>1){
      A <- rdirichlet(STATE, rep(0.25, STATE))
      sVec <- c(s1, rep(NA, TIME-1))
      for(i in 2:length(sVec)){
        sVec[i] <- which.max(rmultinom(1, 1, A[sVec[i-1],]))
      }
    }
    if(STATE==1){sVec <- rep(s1, TIME)}
  }
  
  if(is.null(beta_arr)){
    beta_arr <- lapply(1:STATE, function(x){
      sapply(1:(N_PRED[1]+1), function(y){
        rnorm(BLK, 0, 3)})
    })
  }
  
  if(is.null(X)){
    X <-  if(N_PRED[1] > 0){
      rbind(rep(1, NODE*TIME), 
               sapply(1:(NODE*TIME), function(x){
                 rnorm(N_PRED[1], 0, 2)}))
    } else {
      matrix(1, nrow=1, ncol=NODE*TIME)
    }
    colnames(X) <- paste(rep(1:NODE, TIME), rep(1:TIME, each=NODE), sep="_")
  }
  
  nodes <- 1:NODE
  pi <- lapply(nodes, function(x){
    Xsub <- as.data.frame(X[,paste(x, 1:TIME, sep="_")])
    debetas <- sapply(1:ncol(Xsub), function(y){
      ebeta <- exp(t(Xsub[,y]) %*% t(beta_arr[[sVec[y]]]))
      return(rdirichlet(1, ebeta))
    })
    colnames(debetas) <- paste(x, 1:TIME, sep="_")
    return(debetas)
  })
  pi.mat <- t(do.call(cbind, pi))
  
  if(is.null(gamma_vec)){
    gamma_vec <- rnorm(N_PRED[2], -1, 1)
  }
  
  if(is.null(Z)){
    if(DIRECTED){dy <- expand.grid(1:NODE, 1:NODE)}
    if(!DIRECTED){dy <- t(combn(1:NODE, 2))}
    Z <- do.call(rbind, replicate(TIME, dy, simplify=FALSE))
    Z <- cbind(Z, rep(1:TIME, each = nrow(Z)/TIME))
    for(d in 1:N_PRED[2]){
      Z <- cbind(Z, rnorm(nrow(Z), 1, 1))
    }
    colnames(Z) <- c("node1", "node2", "time", paste("V", 1:N_PRED[2], sep=""))
  }
  
  z <- apply(Z, 1, function(x){
    which.max(rmultinom(1, 1, pi.mat[paste(x[1], x[3], sep="_"),]))
  })
  
  w <- apply(Z, 1, function(x){
    which.max(rmultinom(1, 1, pi.mat[paste(x[2], x[3], sep="_"),]))
  })
  
  if(is.null(B)){
    B <- matrix(rnorm(BLK*BLK), BLK, BLK)
    if(!DIRECTED){
      require(Matrix)
      B <- forceSymmetric(B)}
  }
  
  dgam <- as.matrix(Z[,-c(1:3)]) %*% as.matrix(gamma_vec)
  theta <- lapply(B, function(x){exp(x + dgam) / (1 + exp(x + dgam))}) #element1 for B[1,1], 2 for B[2,1], etc.
  
  prob.edge <- mapply(function(a, b, c){theta[[which(B==B[a,b])]][c]}, a=z, b=w, c=1:nrow(Z))
  Y <- sapply(prob.edge, function(x){rbinom(n=1, size=1, prob=x)})
  
  dyad.data <- as.data.frame(Z)
  dyad.data$Y <- Y
  dyad.data$grp1 <- z
  dyad.data$grp2 <- w
  
  monad.data <- as.data.frame(t(X))
  monad.data$node <- as.numeric(unlist(lapply(strsplit(rownames(monad.data), "_"), "[[", 1)))
  monad.data$time <- as.numeric(unlist(lapply(strsplit(rownames(monad.data), "_"), "[[", 2)))
  pi.mat2 <- pi.mat[colnames(X),]
  for(i in 1:BLK){
    n <- paste("pi",i,sep="")
    monad.data[,n] <- pi.mat2[,i]
  }
  
  return(list(BLK = BLK
              ,NODE = NODE
              ,STATE = STATE
              ,TIME = TIME
              ,DYAD_PRED = N_PRED
              ,MONAD_PRED = N_PRED
              ,DIRECTED = DIRECTED
              ,Y = Y
              ,sVec = sVec
              ,B = B
              ,beta_arr = beta_arr
              ,gamma_mat = gamma_vec
              ,dyad.data = dyad.data
              ,monad.data = monad.data
              ,X = X
              ,Z = Z
              ,pi_vecs = pi.mat2
              #,alpha = alpha
              #,grp = sapply(theta_grp, function(x)x$grp)
              ,theta = prob.edge
              ,theta_full = theta
  ))

  
}
                    
      