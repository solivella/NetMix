#######################
## Data generation 1:
## dynamic network.
#######################


NetSim <- function(BLK = 3, NODE = 5, STATE = 2, TIME = 10, DIRECTED = TRUE, N_PRED = c(0,3),
                      B_t = NULL, A_orig = NULL, beta_arr, gamma_vec, alpha_conc = 0){
  library(expm, quietly=TRUE, warn.conflicts=FALSE)
  stopifnot()
  stopifnot(nrow(beta_arr)==N_PRED+1&ncol(beta_arr)==BLK)
  beta_arr[,1,] <- rep(0, (N_PRED+1)*STATE) 
  if(N_PRED>0){
    stopifnot(length(gamma_vec)==N_PRED)
  }
  ## Transition Mat
  if(is.null(A_orig)){
    A_raw <- matrix(runif(STATE^2, 0.5, 0.5), ncol=STATE)
    A_raw <- as.matrix(band(A_raw, 0, 1))
    A_orig <- prop.table(A_raw,1)
  }
  stopifnot((nrow(A_orig)==ncol(A_orig))&&(ncol(A_orig)==STATE))
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
      kappa_mat[i, ] <- A[sVec[i-1],] 
      sVec[i] <- which(rmultinom(1,1,A_orig[sVec[i-1],])==1)
      A <- A%*%A
    }
  }

  ## Simulate N_PRED ind. covariates per node per time
  if(N_PRED==0){
    baseMatInd <- array(1, c(NODE, 1))
  } else{
    baseMatInd <- cbind(1, runif(NODE), runif(NODE), runif(NODE))
  }
  #baseMateInd <- scale(baseMatInd)
  X <- replicate(TIME,baseMatInd, simplify = FALSE)
  ## Simulate N_PRED dyadic covariates per dyad per time
  baseMatDyad <- if(DIRECTED){
    if(N_PRED==0){
      array(0, c(NODE^2, 1))
    } else {
      cbind(runif(NODE^2), runif(NODE^2), runif(NODE^2))
    }
  } else {
    if(N_PRED==0){
      array(0, c((1 + NODE) * NODE / 2, 1))
    } else {
      cbind(runif((1 + NODE) * NODE / 2), runif((1 + NODE) * NODE / 2),  runif((1 + NODE) * NODE / 2))
    }
  }
  #baseMatDyad <- scale(baseMatDyad)
  Z <- if(DIRECTED){
    replicate(TIME, baseMatDyad, simplify = FALSE)
  } else {
    replicate(TIME, baseMatDyad, simplify = FALSE)
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
                    a <- exp(X[[t]] %*% beta_arr[,,sVec[t]])
                    prop.table(a, 1) * alpha_conc
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
  
  theta <- lapply(1:TIME
                  ,function(x){
                    cats <- BLK^2
                    grp_cdf <- apply(pi_vecs[[x]]%x%pi_vecs[[x]],1,cumsum)
                    u <- rep(runif(NODE^2), rep(cats,NODE^2))
                    grp <- c(1:cats)[1 + apply(u > grp_cdf, 2, sum)]
                    if(!DIRECTED)
                      grp <- grp[which(lower.tri(diag(NODE), diag = TRUE))]
                    1/(1 + exp(- B_t[grp] - (Z[[x]] %*% gamma_vec)))
                  }
  )
  
  ## Simulate Y
  Y <- lapply(theta
              ,function(th){
                u <- if(DIRECTED){
                  runif(NODE^2)
                } else {
                  runif((1 + NODE)*NODE/2)
                }
                as.numeric(th > u)
              })
  dyad.data <- as.data.frame((do.call(rbind, Z)))
  dyad.data$Y <- unlist(Y)
  
  year.names <- 1:TIME
  dyad.data$year <- rep(year.names, each = nrow(dyad.data)/TIME)
  
  node.names <- 1:NODE
  if(DIRECTED){
    dyad.data$sender <- rep(node.names, each = NODE, times = TIME)
    dyad.data$receiver <- rep(node.names, times = TIME * NODE)
    
  } else {
    sr <- t(outer(node.names,node.names, FUN=paste, sep=":"))
    sr_und <- do.call(rbind,strsplit(sr[lower.tri(sr, diag = TRUE)],":"))
    dyad.data$sender <- rep(sr_und[,1], times = TIME)
    dyad.data$receiver <- rep(sr_und[,2], times = TIME)
  }
  if(N_PRED!=0){
  monad.data <- as.data.frame((do.call(rbind, X))[,-1])
  names(monad.data) <- paste("V",1:ncol(monad.data), sep="")
  monad.data$node <- rep(node.names, times = TIME)
  monad.data$year <- rep(year.names, each = NODE)
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
              ,DYAD_PRED = N_PRED
              ,MONAD_PRED = N_PRED
              ,DIRECTED = DIRECTED
              ,Y = Y
              ,kappa_mat = kappa_mat
              ,sVec = sVec
              ,B = t(B_t)
              ,beta_arr = beta_arr
              ,gamma_mat = gamma_vec
              ,dyad.data = dyad.data
              ,monad.data = monad.data
              ,X = X
              ,pi_vecs = pi_vecs
              ,A_orig = A_orig
              ,alpha = alpha
              ,theta = theta
              ,theta_full = theta_all#array(unlist(theta_all),c(BLK,BLK,(NODE^2)*TIME))
  ))
}
