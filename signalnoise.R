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
bm<- matrix(qlogis(c(0.6, 0.10, 0.01, 0.5)), ncol = 2) 
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

beta_easy = list(array(c(-2.5, -2.5, ##Intercepts
                         0.5, -0.5), ## Predictor coefficients
                       c(2, 2)),
                 array(c(1, -1,
                         1, -1),
                       c(2, 2)))

# three blocks?
#beta_easy = list(array(c(-4.5, -4.5,-4.5 ##Intercepts
#), ## Predictor coefficients
#c(3, 1)),
#array(c(1, -1,-1),
#      c(3, 1)))


sVec<-rep(c(1,2),c(25,25))

nbS <- 2 # senators have 2 groups
nbB <- 2

########################################
## Generate a dyadic bipartite network
########################################

## Sample monadic predictors
nnS <- 400
nnB <-200
TIME <- 50


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
    ebeta <- exp(t(Xsub[,y]) %*% t(beta_easy[[sVec[y]]]))
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

split_strings <- strsplit(rownames(piS), "_")
S_period <- lapply(split_strings, function(x) x[2]) %>% unlist() %>% as.numeric()
piS <- piS %>% data.frame() %>% mutate(period = S_period)


piB <- lapply(1:nnB, function(x){
  Xsub <- as.data.frame(XB_array[,x,])#[1,] remove # to get no monadic cov
  debetas <- sapply(1:ncol(Xsub), function(y){
    ebeta <- exp(t(Xsub[,y]) %*% t(beta_easy[[sVec[y]]]))
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

split_strings <- strsplit(rownames(piB), "_")
B_period <- lapply(split_strings, function(x) x[2]) %>% unlist() %>% as.numeric()
piB <- piB %>% data.frame() %>% mutate(period = B_period)

#------------------------------------------
## Signal-to-Noise Ratio
#----------------------------------------

ratio_calc <- function(bm, piS, piB, nnS, nnB){
  
 p11_1 <- piS %>% filter(period <= 25) %>% select(1) %>% colMeans()
 p12_1 <- piS %>% filter(period <= 25) %>% select(2) %>% colMeans()

  P1_1 <- sqrt(diag(c(nnS*p11_1, nnS*p12_1)))
  
  p21_1 <- piS %>% filter(period <= 25) %>% select(1) %>% colMeans()
  p22_1 <- piS %>% filter(period <= 25) %>% select(2) %>% colMeans()
  P2_1 <- sqrt(diag(c(nnB*p21_1, nnB*p22_1)))
  
  #info1 <- svd(P1_1 %*% bm %*% P2_1)$d[2]^2/svd(P1_1 %*% bm %*% P2_1)$d[1]
  PQ1<-P1_1 %*% bm %*% t(P2_1)
  eigen_res1 <- eigen(PQ1)
  eigenval1 <- eigen_res1$values
  eigenval1 <- unique(eigenval1)
  eigenval1 <- sort(eigenval1, decreasing = TRUE, method = "radix")
  info1<-(eigenval1[2])^2/eigenval1[1]
  
  p11_2 <- piS %>% filter(period > 25) %>% select(1) %>% colMeans()
  p12_2 <- piS %>% filter(period > 25) %>% select(2) %>% colMeans()
  P1_2 <- sqrt(diag(c(nnS*p11_2, nnS*p12_2)))
  
  p21_2 <- piS %>% filter(period > 25) %>% select(1) %>% colMeans()
  p22_2 <- piS %>% filter(period > 25) %>% select(2) %>% colMeans()
  P2_2 <- sqrt(diag(c(nnB*p21_2, nnB*p22_2)))
  
  #info2 <- svd(P1_2 %*% bm %*% P2_2)$d[2]^2/svd(P1_2 %*% bm %*% P2_2)$d[1]
  PQ2<-P1_2 %*% bm %*% t(P2_2)
  eigen_res2 <- eigen(PQ2)
  eigenval2 <- eigen_res2$values
  eigenval2 <- unique(eigenval2)
  eigenval2 <- sort(eigenval2, decreasing = TRUE, method = "radix")
  info2<-(eigenval2[2])^2/eigenval2[1]
  
  return(c(info1, info2))
}

bm<-mat
ratio_calc(bm,piS, piB, nnS, nnB)
#[1] 93.02992 66.61547


#-------------------------------------------
# Second network
#---------------------------------------
bm<- matrix(qlogis(c(0.75, 0.10, 0.01, 0.5)), ncol = 2) 
beta_easy = list(array(c(-2.5, -2.5, ##Intercepts
                         0.5, -0.5), ## Predictor coefficients
                       c(2, 2)),
                 array(c(1, -1,
                         1, -1),
                       c(2, 2)))


sVec<-rep(c(1,2),c(25,25))

nbS <- 2 # senators have 2 groups
nbB <- 2

########################################
## Generate a dyadic bipartite network
########################################

## Sample monadic predictors
nnS <- 400
nnB <-200
TIME <- 50


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
    ebeta <- exp(t(Xsub[,y]) %*% t(beta_easy[[sVec[y]]]))
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

split_strings <- strsplit(rownames(piS), "_")
S_period <- lapply(split_strings, function(x) x[2]) %>% unlist() %>% as.numeric()
piS <- piS %>% data.frame() %>% mutate(period = S_period)


piB <- lapply(1:nnB, function(x){
  Xsub <- as.data.frame(XB_array[,x,])#[1,] remove # to get no monadic cov
  debetas <- sapply(1:ncol(Xsub), function(y){
    ebeta <- exp(t(Xsub[,y]) %*% t(beta_easy[[sVec[y]]]))
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

split_strings <- strsplit(rownames(piB), "_")
B_period <- lapply(split_strings, function(x) x[2]) %>% unlist() %>% as.numeric()
piB <- piB %>% data.frame() %>% mutate(period = B_period)

ratio_calc(bm,piS, piB, nnS, nnB)
#[1]66.77685 34.91699

#Easy 
bm<-matrix((c(0.85, 0.01, 0.01, 0.99)), ncol = 2)
c<-ratio_calc(bm,piS, piB, nnS, nnB)
c[2]^2/c[1]
#Medium
bm<-matrix((c(0.65, 0.20, 0.35, 0.75)), ncol = 2)
c<-ratio_calc(bm,piS, piB, nnS, nnB)
c[2]^2/c[1]
#Hard
bm<-matrix((c(0.65, 0.5, 0.4, 0.45)), ncol = 2)
c<-ratio_calc(bm,piS, piB, nnS, nnB)
c[2]^2/c[1]



bm_list <- list()
interpolate <- function(start, end, length.out = 50) {
  seq(from = start, to = end, length.out = length.out)
}

a <- interpolate(0.85, 0.65)
b <- interpolate(0.01, 0.5)
c <- interpolate(0.01, 0.4)
d <- interpolate(0.99, 0.45)


for (i in 1:50) {
  bm <- matrix(c(a[i], b[i], c[i], d[i]), ncol = 2)
  bm_list[[i]] <- bm
}


bm_list[[1]]  # The first BM with most different entries
bm_list[[50]] # The last BM with almost identical entries

range<-1:50
results <- lapply(1:length(range), function(x) ratio_calc(bm_list[[x]], 
                                                          piS = piS, piB = piB, nnS = 400, nnB = 200) ) %>% 
  unlist() %>% matrix(ncol = 2, byrow = TRUE) %>% data.frame() %>% 
  mutate(entry = range) %>% 
  rename(info1 = X1, info2 = X2)

results %>% ggplot(aes(x = entry)) +
  geom_line(aes(y = info1, color = "state1")) +
  geom_line(aes(y = info2, color = "state2"))
