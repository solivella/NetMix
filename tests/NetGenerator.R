#######################
## Data generation 1:
## dynamic network.
#######################


NetSim <- function(BLK = 3,
                   NODE = 5,
                   STATE = 2,
                   TIME = 10,
                   DIRECTED = TRUE,
                   N_PRED = c(1,1),
                   B_t = NULL, 
                   sVec = NULL,
                   beta_arr,
                   gamma_vec,
                   X = NULL,
                   Z = NULL){
  #stopifnot(nrow(beta_arr)==N_PRED+1&ncol(beta_arr)==BLK)
  # if(N_PRED>0){
  #   stopifnot(length(gamma_vec)==N_PRED)
  # }
  NDYADS <- ifelse(DIRECTED, NODE^2, (1 + NODE) * NODE / 2)
  # if(!is.null(X)|!is.null(Z)){
  #   if(!is.list(X)|length(X)!=TIME|any(sapply(X,function(x)ncol(x)!=N_PRED+1))|any(sapply(X,function(x)nrow(x)!=NODE))){
  #     stop("X must be a list of length TIME, each element with N_PRED+1 columns (constant leading) and NODE rows.")
  #   }
  #   if(!is.list(Z)|length(Z)!=TIME|any(sapply(Z,function(x)ncol(x)!=N_PRED))|any(sapply(Z,function(x)nrow(x)!=NDYADS)))
  #     stop("Z must be a list of length TIME, each element with N_PRED columns (no constant) and NDYADS rows.")
  # }
  
  ## Transition Mat
  if(is.null(sVec)){
    A_raw <- matrix(runif(STATE^2, 0.5, 0.5), ncol=STATE)
    #A_raw <- as.matrix(band(A_raw, 0, 1))
    A_orig <- prop.table(A_raw,1)
    A <- A_orig
    
    stopifnot(ncol(B_t)==BLK)
    
    ## Simulate TIME steps
    piInit <- rep(0,STATE)
    piInit[1] <- 1
    sVec <- array(NA,TIME)
    sVec[1] <- 1
    kappa_mat <- array(NA, c(TIME, STATE))
    kappa_mat[1, ] <- piInit
    if(TIME > 1){
      for(i in 2:(TIME)){
        kappa_mat[i, ] <- piInit %*% SocNet:::.mpower(A_orig, i) 
        sVec[i] <- which(rmultinom(1,1,kappa_mat[i, ])==1)
      }
    }
  }
  
  
  ## Simulate N_PRED ind. covariates per node per time
  if(is.null(X)){
    if(N_PRED[1]==0){
      baseMatInd <- array(1, c(NODE, 1))
    } else{
      baseMatInd <- replicate(N_PRED[1], rnorm(NODE))
    }
    # baseMateInd <- scale(baseMatInd)
    X <- lapply(1:TIME,function(x){cbind(1, apply(baseMatInd, 2, function(x)rnorm(length(x), x, 0.5)))})
  }
  ## Simulate N_PRED dyadic covariates per dyad per time
  if(is.null(Z)){
    baseMatDyad <- 
      if(N_PRED[2]==0){
        array(0, c(NDYADS, 1))
      } else {
        replicate(N_PRED[2], rnorm(NDYADS))
      }
    #baseMatDyad <- scale(baseMatDyad)
    Z  <- lapply(1:TIME,function(x){apply(baseMatDyad, 2, function(x)rnorm(length(x), x, 0.5))})
  }
  
  
  ## Simulate blockmodel
  if(is.null(B_t)){
    B_t <- if(DIRECTED){
      matrix(rnorm(BLK^2, 0, 1), ncol = BLK)
    } else {
      forceSymmetric(matrix(rnorm(BLK^2, 0, 1), ncol = BLK))
    }
  }
  
  
  ## Form alpha
  alpha <- lapply(1:TIME
                  ,function(t){
                    a <- exp(X[[t]] %*% t(beta_arr[[sVec[t]]]))
                  })
  
  ## Simulate pi
  pi_vecs <- lapply(alpha
                    , function(a){
                      prop.table(matrix(rgamma(NODE*BLK,c(a)),ncol=BLK),1)
                    })
  
  ## Simulate theta
  theta_all <- lapply(1:TIME,
                      function(x){
                        part <- Z[[x]] %*% gamma_vec
                        sapply(part,function(i)1/(1 + exp(-B_t - i)))
                      })
  
  theta_grp <- lapply(1:TIME
                      ,function(x){
                        cats <- BLK^2
                        grp_cdf <- apply(pi_vecs[[x]]%x%pi_vecs[[x]],1,cumsum)
                        u <- rep(runif(NODE^2), rep(cats,NODE^2))
                        grp <- c(1:cats)[1 + apply(u > grp_cdf, 2, sum)]
                        if(!DIRECTED){
                          grp <- grp[which(lower.tri(diag(NODE), diag = TRUE))]
                        }
                        list(theta = 1/(1 + exp(- B_t[grp] - (Z[[x]] %*% gamma_vec)))
                             ,grp = grp)
                      }
  )
  
  ## Simulate Y
  Y <- lapply(theta_grp
              ,function(th){
                u <- if(DIRECTED){
                  runif(NDYADS)
                } else {
                  runif(NDYADS)
                }
                as.numeric(th$theta > u)
              })
  dyad.data <- as.data.frame((do.call(rbind, Z)))
  dyad.data$Y <- unlist(Y)
  
  year.names <- 1:TIME
  dyad.data$time <- rep(year.names, each = nrow(dyad.data)/TIME)
  dyad.data$grp <- unlist(lapply(theta_grp,function(x)x$grp))
  
  node.names <- 1:NODE
  if(DIRECTED){
    dyad.data$node1 <- rep(node.names, each = NODE, times = TIME)
    dyad.data$node2 <- rep(node.names, times = TIME * NODE)
    
  } else {
    sr <- t(outer(node.names,node.names, FUN=paste, sep=":"))
    sr_und <- do.call(rbind,strsplit(sr[lower.tri(sr, diag = TRUE)],":"))
    dyad.data$node1 <- rep(sr_und[,1], times = TIME)
    dyad.data$node2 <- rep(sr_und[,2], times = TIME)
  }
  if(N_PRED[1]!=0){
    monad.data <- as.data.frame((do.call(rbind, X))[,-1])
    names(monad.data) <- paste("V",1:ncol(monad.data), sep="")
    monad.data$node <- rep(node.names, times = TIME)
    monad.data$time <- rep(year.names, each = NODE)
    all_pi <- do.call(rbind, pi_vecs)
    for(i in 1:BLK){
      n <- paste("pi",i,sep="")
      monad.data[,n] <- all_pi[,i]
    }
  } else {
    monad.data  <- NULL
  }
  
  return(list(BLK = BLK
              ,NODE = NODE
              ,STATE = STATE
              ,TIME = TIME
              ,DYAD_PRED = N_PRED[2]
              ,MONAD_PRED = N_PRED[1]
              ,DIRECTED = DIRECTED
              ,Y = Y
              #,kappa_mat = kappa_mat
              ,sVec = sVec
              ,B = t(B_t)
              ,beta_arr = beta_arr
              ,gamma_mat = gamma_vec
              ,dyad.data = dyad.data
              ,monad.data = monad.data
              ,X = X
              ,Z = Z
              ,pi_vecs = do.call(rbind, pi_vecs)
              #,A_orig = A_orig
              ,alpha = alpha
              ,grp = sapply(theta_grp, function(x)x$grp)
              ,theta = lapply(theta_grp,function(x)x$theta)
              ,theta_full = theta_all#array(unlist(theta_all),c(BLK,BLK,(NODE^2)*TIME))
  ))
}
