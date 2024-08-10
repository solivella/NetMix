rm(list=ls())
set.seed(02138)
require(pacman)
p_load(LaplacesDemon,MCMCpack,igraph,ggplot2,tidyverse,gridExtra,NetMix)



###########################
## initialize parameters
###########################

## All cases (we will use the easy case here)
## BM_list <- list(BM_easy = matrix(qlogis(c(0.85, 0.01, 0.01, 0.99)), ncol = 2),
##                BM_med = matrix(qlogis(c(0.65, 0.2, 0.35, 0.75)), ncol = 2),
##                BM_hard = matrix(qlogis(c(0.65, 0.5, 0.4, 0.45)), ncol = 2))


## beta_list <- list(beta_easy = list(array(c(-4.5, -4.5, ##Intercepts
##                                           0.0, 0.0), ## Predictor coefficients
##                                         c(2, 2)),
##                                   array(c(-4.5, -4.5,
##                                           0.0, 0.0),
##                                         c(2, 2))),
##                  beta_med =list(array(c(0.05, 0.75,
##                                         -0.75,  -1.0),
##                                       c(2, 2)),
##                                 array(c(0.05, 0.75,
##                                         -0.75,  -1.0),
##                                       c(2, 2))),
##                  beta_hard =list(array(c(0.0, 0.0,
##                                          -0.75,  -1.0),
##                                        c(2, 2)),
##                                  array(c(0.0, 0.0,
##                                          -0.75,  -1.0),
##                                        c(2, 2))))

#BM_easy = matrix(qlogis(c(0.9, 0.2, 0.05, 0.35)), ncol = 2) #2*2

#BM_easy = matrix(qlogis(c(0.9, 0.2, 0.05, 0.35)), ncol = 2) #easy
#BM_easy = matrix(qlogis(c(0.9, 0.45, 0.20, 0.35)), ncol = 2) #medium
#BM_easy = matrix(qlogis(c(0.9, 0.75, 0.60, 0.35)), ncol = 2) #hard

#BM_easy = matrix(qlogis(c(0.7, 0.2, 0.05, 0.35)), ncol = 2) #ok
 BM_easy = matrix(qlogis(c(0.6, 0.10, 0.01, 0.5)), ncol = 2) 
#BM_easy =matrix(qlogis(c(0.9, 0.2, 0.05, 0.35,0.02,0.8,0.75,0.3,0.65)), ncol =3) #3*3 BM used
#BM_easy =matrix(qlogis(c(0.9, 0.2, 0.05, 0.35,0.02,0.8,0.75,0.3,0.5)), ncol =3)

#BM_easy = matrix(qlogis(c(0.85, 0.01, 0.01, 0.99)), ncol = 2) #easy case

#BM_easy = matrix(qlogis(c(0.65, 0.2, 0.35, 0.75)), ncol = 2) #2*2
#BM_easy = matrix(qlogis(c(0.9, 0.2,0.5, 0.05, 0.35,0.8)), ncol = 2) #2*3
#BM_easy = matrix(qlogis(c(0.85, 0.05, 0.10, 0.35)), ncol = 2) #m3


#BM_easy = matrix(qlogis(c(0.65, 0.35, 0.20, 0.75)), ncol = 2)
#beta_easy = list(array(c(-2.5, -2.5 ##Intercepts
#), ## Predictor coefficients
#c(2, 1)),
#array(c(1, -1),
#      c(2, 1)))


#beta_easy = list(array(c(-2.5, -2.5, ##Intercepts
#                         0.5, -0.5), ## Predictor coefficients
#                       c(2, 2)),
#                 array(c(1, -1,
#                         1, -1),
#                       c(2, 2)))


beta_easy = list(array(c(0.05, 0.75, 
                         -0.75, -1.0), 
                       c(2, 2)),
                 array(c(-0.05, 0.55,
                         -0.75, 0.75),
                       c(2, 2)))




# easier beta with covariates
#beta_easy = list(array(c(-2.5,-2.5,-2.5,#easy beta V7
#                         -1.75,-2.5,1), 
#                       c(2, 3)),
#                 array(c(-1.75,-0.35,-1.5,0.5,-2,-1.75),
#                       c(2, 3)))

#sVec<-rep(c(1),c(10))
#sVec<-rep(c(1),c(1))
sVec<-rep(c(1,2),c(5,5))
nbS <- 2 # senators have 2 groups
nbB <- 2 # bills have 3 groups 

########################################
## Generate a dyadic bipartite network
########################################

## Sample monadic predictors
nnS <- 5
nnB <-5
TIME <- 10


XS <- matrix(NA, nrow = TIME, ncol = nnS)
XS[1,] <- rnorm(nnS,0,1)
if(TIME > 1){
  for(t in 2:TIME){
    #  XS[t,] <- rnorm(1, XS[t-1,], 1)
    XS[t,] <- rnorm(nnS,0,1) #no time dependency
  }
}
colnames(XS) <- paste0("S", 1:nnS)


XS_array <- array(0, dim = c(2, ncol(XS), nrow(XS)))
for (i in 1: nrow(XS)) {
  XS_array[1,,i] <- 1
  XS_array[2,,i] <- XS[i,]
}


XB <- matrix(NA, nrow = TIME, ncol = nnB)
XB[1,] <- rnorm(nnB,0,1)

if(TIME > 1){
  for(t in 2:TIME){
    # XB[t,] <- rnorm(1, XB[t-1,], 1)
    XB[t,] <- rnorm(nnB,0,1) #no time dependency
  }
}
colnames(XB) <- paste0("B", 1:nnB)
XB_array <- array(0, dim = c(2, ncol(XB), nrow(XB)))
for (i in 1: nrow(XB)) {
  XB_array[1,,i] <- 1
  XB_array[2,,i] <- XB[i,]
}

## Sample mixed membership for each family
piS <- lapply(1:nnS, function(x){
  Xsub <- as.data.frame(XS_array[,x,])#[1,] remove # to get no monadic cov
  debetas <- sapply(1:ncol(Xsub), function(y){
    ebeta <- exp(t(Xsub[,y]) %*% (beta_easy[[sVec[y]]])) #note: I removed t before beta_easy
    pi_vec <- rdirichlet(1, ebeta)
    if(anyNA(pi_vec)){
      cat_ind <- which.max(ebeta)
      pi_vec <- rep(0, length(ebeta))
      pi_vec[cat_ind] <- 1
    } 
    return(pi_vec)
  })
  colnames(debetas) <- paste(paste0("S", x), 1:TIME, sep="_")
  return(debetas)
})
piS <- t(do.call(cbind, piS))

piB <- lapply(1:nnB, function(x){
  Xsub <- as.data.frame(XB_array[,x,])#[1,] remove # to get no monadic cov
  debetas <- sapply(1:ncol(Xsub), function(y){
    ebeta <- exp(t(Xsub[,y]) %*% (beta_easy[[sVec[y]]])) #note: I removed t before beta_easy
    pi_vec <- rdirichlet(1, ebeta)
    if(anyNA(pi_vec)){
      cat_ind <- which.max(ebeta)
      pi_vec <- rep(0, length(ebeta))
      pi_vec[cat_ind] <- 1
    } 
    return(pi_vec)
  })
  colnames(debetas) <- paste(paste0("B", x), 1:TIME, sep="_")
  return(debetas)
})
piB <- t(do.call(cbind, piB))


## Sample dyadic predictor

Z <- do.call(rbind, replicate(TIME, expand.grid(1:nnS, 1:nnB), simplify=FALSE))
Z <- cbind(Z, rep(1:TIME, each = nrow(Z)/TIME), matrix(NA, nrow(Z), 1))
t1 <- list(rnorm(nrow(Z)/TIME, 0, 2))  #only one predictor for now 
#t1 <- rep(0,rnorm(nrow(Z)/TIME))
if(TIME > 1){
  for(t in 2:TIME){
    t1[[t]] <- t1[[t-1]] + rnorm(length(t1[[t-1]]), 0, 1)
  }
}
Z[,1+3] <- do.call(c, t1)
colnames(Z) <- c("S", "B", "year", paste("V", 1:1, sep=""))

## Sample dyad-specific block memberships 
z <- apply(Z, 1, function(x){
  which.max(rmultinom(1, 1, piS[paste0("S",paste(x[1], x[3], sep="_")),])) #S
})

w <- apply(Z, 1, function(x){
  which.max(rmultinom(1, 1, piB[paste0("B",paste(x[2], x[3], sep="_")),])) #B
})

## Sample network edges

gamma_vec = c(0.1)
dgam <- (as.matrix(Z[,-c(1:3)])) %*% as.matrix(gamma_vec) 
#dgam<-as.matrix(rep(0,nrow(Z))) #use this for no dyadic pred
theta <- lapply(BM_easy, function(x){exp(x + dgam) / (1 + exp(x + dgam))}) #element1 for B[1,1], 2 for B[2,1], etc.
prob.edge <- mapply(function(a, b, c){theta[[which(BM_easy==BM_easy[a,b])[1]]][c]}, a=z, b=w, c=1:nrow(Z))
Y <- sapply(prob.edge, function(x){rbinom(n=1, size=1, prob=x)}) 


## List Generator

df_dyad_1 <- data.frame(Y = Y,
                        var1 = Z[,4],  # comment out for no dyad cov
                        id1 = paste0("S",Z[,1]),
                        id2 = paste0("B", Z[,2]),
                        year=Z[,3])
df_monad_S0 <- as.data.frame(XS)
df_monad_B0 <- as.data.frame(XB)

df_monad_S <- df_monad_S0 %>% 
  rownames_to_column(var = "year") %>%  
  pivot_longer(1:nnS+1, names_to = "id", values_to = "VarS1")#%>%
#select(-VarS1) remove # to get no monadic cov


df_monad_B <- df_monad_B0 %>% 
  rownames_to_column(var = "year") %>%  
  pivot_longer(1:nnB+1, names_to = "id", values_to = "VarB1")#%>%
#select(-VarB1) remove # to get no monadic cov

piS <- piS[paste(df_monad_S$id, df_monad_S$year, sep = "_"),]
piB <- piB[paste(df_monad_B$id, df_monad_B$year, sep = "_"),]

for(i in 1:nbS){
  n <- paste0("piS", i)
  df_monad_S[,n] <- piS[,i]
}

for(i in 1:nbB){
  n <- paste("piB",i,sep="")
  df_monad_B[,n] <- piB[,i]
}

netSim <- list(piS = piS, piB = piB,
               df_dyad_1=df_dyad_1,
               df_monad_B=df_monad_B,df_monad_S=df_monad_S,
               nbS = nbS, nbB = nbB, nnS = nnS, nnB = nnB,
               sVec=sVec,TIME=TIME)

#beta_var_array <- array(c(0.0001, 1, ##Intercepts
#                          0.0001, 1,
#                          0.0001, 1,
#                          0.0001, 1), ## Predictor coefficients
#                        c(2, 2, 2))
beta_var_array <- array(c(0.0001, 1, ##Intercepts
                         # 0.0001, 1,
                          0.0001, 1,
                        0.0001, 1,
                       #   0.0001, 1, ##Intercepts
                          0.0001, 1), ## Predictor coefficients
                        c(2, 2, 2)) #ncov, ngroup, nstate

#beta_var_array <- array(c(0.0001, ##Intercepts
#                          0.0001, 
#                          0.0001,
#                          0.0001, 
#                          0.0001, ##Intercepts
#                          0.0001 ), ## Predictor coefficients
#                        c(1, 3, 2))
##############################
## Fitting the simulated data
##############################

res_dynbi <- mmsbm(formula.dyad = Y~var1, #Y~1 for no monadic cov
                   formula.monad = list(~VarS1, ~VarB1), #comment out for no dyadic cov
                   #formula.monad = list(~1, ~1), # no monadic covariates
                   timeID="year",
                   senderID = "id1",
                   receiverID = "id2",
                   nodeID = list("id","id"),
                   bipartite= TRUE,
                   data.dyad = netSim[["df_dyad_1"]],
                   data.monad = list(netSim[["df_monad_S"]],netSim[["df_monad_B"]]),
                   n.blocks = c(2,2),n.hmmstates =2,
                   mmsbm.control = list(verbose = TRUE,
                                        threads=1,
                                        svi = TRUE,
                                        vi_iter = 10000,
                                        # batch_size = 1.0,
                                        conv_tol = 1e-3,
                                        #var_beta=list(c(0.01),c(0.01)),
                                        var_beta=list(beta_var_array, 
                                                     beta_var_array),
                                        hessian = TRUE))

vcov_monad1<-res_dynbi$vcov_monad1
vcov_monad2<-res_dynbi$vcov_monad2

coef_monad1<-res_dynbi$MonadCoef1
coef_monad2<-res_dynbi$MonadCoef2
coef_monad1
coef_monad2
