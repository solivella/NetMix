#########################
## Bipartite simulations
## Two monadic predictors
## One dyadic predictor.
########################

##Workspace setup
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


sVec<-rep(c(1,2),c(5,5))

nbS <- 2 # senators have 2 groups
nbB <- 2

########################################
## Generate a dyadic bipartite network
########################################

## Sample monadic predictors
nnS <- 100
nnB <-50
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



##############################
## Fitting the simulated data
##############################

res_dynbi <- mmsbm(formula.dyad = Y~var1, #Y~1 for no monadic cov
                   formula.monad = list(~VarS1, ~VarB1), #ocmment out for no dyadic cov
                   timeID="year",
                   senderID = "id1",
                   receiverID = "id2",
                   nodeID = list("id","id"),
                   bipartite= TRUE,
                   data.dyad = netSim[["df_dyad_1"]],
                   data.monad = list(netSim[["df_monad_S"]],netSim[["df_monad_B"]]),
                   n.blocks = c(2,2),n.hmmstates = 2,
                   mmsbm.control = list(verbose = TRUE,
                                        threads=1,
                                        svi = FALSE,
                                        vi_iter = 2000,
                                        batch_size = 1.0,
                                        conv_tol = 1e-3,
                                        var_beta=list(c(0.01),
                                                      c(0.01)),
                                        hessian = FALSE))

## statistics for S
loss.mat.phi<- res_dynbi$MixedMembership1 %*% netSim$piS
phi_ord <- clue::solve_LSAP(t(loss.mat.phi), TRUE)
orig_mm_names <- rownames(res_dynbi$MixedMembership1)
res_dynbi$MixedMembership1 <- res_dynbi$MixedMembership1[phi_ord, ]
rownames(res_dynbi$MixedMembership1) <- orig_mm_names
orig_g_names <- dimnames(res_dynbi$BlockModel)
res_dynbi$BlockModel <- res_dynbi$BlockModel[phi_ord, phi_ord]
dimnames(res_dynbi$BlockModel) <- orig_g_names #comment out following lines if no covariates
loss.mat.kappa <- res_dynbi$Kappa %*% model.matrix(~-1+as.factor(netSim$sVec))
kappa_ord <- clue::solve_LSAP((loss.mat.kappa), TRUE)
orig_coef_names <- dimnames(res_dynbi$MonadCoef1)
res_dynbi$MonadCoef1 <- res_dynbi$MonadCoef1[, phi_ord, kappa_ord, drop = FALSE]
dimnames(res_dynbi$MonadCoef1) <- orig_coef_names
orig_k_names <- rownames(res_dynbi$Kappa)
res_dynbi$Kappa <- res_dynbi$Kappa[kappa_ord, , drop = FALSE]
rownames(res_dynbi$Kappa) <- orig_k_names
orig_tk_names <- dimnames(res_dynbi$TransitionKernel)
res_dynbi$TransitionKernel <- res_dynbi$TransitionKernel[kappa_ord, kappa_ord, drop = FALSE]
dimnames(res_dynbi$TransitionKernel) <- orig_tk_names
pred_piS <- data.frame(#netSimwork = rep(c(type), each=with(SynthnetSims[[i]],NODE*BLK*TIME)),
  NodeType = "Senator",
  Group = factor(rep(c(1:netSim$nbS), each = with(netSim,nnS*TIME))),
  Year = rep(paste("Year",1:netSim$TIME), each=netSim$nnS, times=netSim$nbS),
  State = rep(c("State 1","State 2"),c(25*netSim$nnS,25*netSim$nnS)),
  V2 = rep(netSim$df_monad_S[,"VarS1"], netSim$nbS), #comment out if no covariates
  Truth = c(netSim$piS),
  Pred = c(t(res_dynbi$MixedMembership1[phi_ord,])))

#pred_piS <-pred_piS%>%
#  mutate(switched_Pred=1-Pred)%>%
#  mutate(best_Pred=ifelse((switched_Pred-Truth)^2<= (Pred-Truth)^2, switched_Pred, Pred ))

## statistics for B
loss.mat.phi<- res_dynbi$MixedMembership2 %*% netSim$piB
phi_ord <- clue::solve_LSAP(t(loss.mat.phi), TRUE)
orig_mm_names <- rownames(res_dynbi$MixedMembership2)
res_dynbi$MixedMembership2 <- res_dynbi$MixedMembership2[phi_ord, ]
rownames(res_dynbi$MixedMembership2) <- orig_mm_names
orig_g_names <- dimnames(res_dynbi$BlockModel)
res_dynbi$BlockModel <- res_dynbi$BlockModel[phi_ord, phi_ord]
dimnames(res_dynbi$BlockModel) <- orig_g_names #comment out following lines if no covariates
loss.mat.kappa <- res_dynbi$Kappa %*% model.matrix(~-1+as.factor(netSim$sVec))
kappa_ord <- clue::solve_LSAP((loss.mat.kappa), TRUE)
orig_coef_names <- dimnames(res_dynbi$MonadCoef2)
res_dynbi$MonadCoef2 <- res_dynbi$MonadCoef2[, phi_ord, kappa_ord, drop = FALSE]
dimnames(res_dynbi$MonadCoef2) <- orig_coef_names
orig_k_names <- rownames(res_dynbi$Kappa)
res_dynbi$Kappa <- res_dynbi$Kappa[kappa_ord, , drop = FALSE]
rownames(res_dynbi$Kappa) <- orig_k_names
orig_tk_names <- dimnames(res_dynbi$TransitionKernel)
res_dynbi$TransitionKernel <- res_dynbi$TransitionKernel[kappa_ord, kappa_ord, drop = FALSE]
dimnames(res_dynbi$TransitionKernel) <- orig_tk_names
pred_piB <- data.frame(#netSimwork = rep(c(type), each=with(SynthnetSims[[i]],NODE*BLK*TIME)),
  NodeType = "Bill",
  Group = factor(rep(c(1:netSim$nbB), each = with(netSim,nnB*TIME))),
  Year = rep(paste("Year",1:netSim$TIME), each=netSim$nnB, times=netSim$nbB),
  State = rep(c("State 1","State 2"),c(25*netSim$nnB,25*netSim$nnB)),
  V2 = rep(netSim$df_monad_B[,"VarB1"], netSim$nbB), #comment out if no covariates
  Truth = c(netSim$piB),
  Pred = c(t(res_dynbi$MixedMembership2[phi_ord,])))

#pred_piB <-pred_piB%>%
#  mutate(switched_Pred=1-Pred)%>%
#  mutate(best_Pred=ifelse((switched_Pred-Truth)^2<= (Pred-Truth)^2, switched_Pred, Pred ))

plot(x = pred_piB$Truth, y = pred_piB$Pred, xlab="True Value", ylab="Predicted Membership")
plot(x=pred_piS$Truth, y=pred_piS$Pred, xlab="True Value", ylab="Predicted Membership")

#pred_pi <- rbind(pred_piS, pred_piB)
#plot(x = pred_piB[,5], y = pred_piB[,8], xlab="True Value", ylab="Switched Predicted Membership")
#plot(pred_piS[,5], pred_piS[,8], xlab="True Value", ylab="Switched Predicted Membership")

# How good is the initialization after by-period switching?
mm_init<-res_dynbi$mm_init
init_b<-c(mm_init[[2]][1,],mm_init[[2]][2,])
plot(y=init_b, x=c(netSim$piB[,1],netSim$piB[,2]), xlab="True Value", ylab="Initialized Membership")

init_s<-c(mm_init[[1]][1,],mm_init[[1]][2,])
plot(y=init_s, x=pred_piS[,5], xlab="True Value", ylab="Initialized Membership")

# How good is the initialization before switching?
mm_init<-res_dynbi$mm_orig
init_b<-c(mm_init[[2]][1,],mm_init[[2]][2,])
plot(y=init_b, x=c(netSim$piB[,1],netSim$piB[,2]), xlab="True Value", ylab="Initialized Membership")

init_s<-c(mm_init[[1]][1,],mm_init[[1]][2,])
plot(y=init_s, x=c(netSim$piS[,1],netSim$piS[,2]), xlab="True Value", ylab="Initialized Membership")

#Initialization compared to prediction
#init_b<-c(mm_init[[2]][1,],mm_init[[2]][2,])
#plot(x=init_b, y=pred_piB[,6], ylab="Predicted Membership", xlab="Initialized Membership")

#init_s<-c(mm_init[[1]][1,],mm_init[[1]][2,])
#plot(x=init_s, y=pred_piS[,6], ylab="Predicted Membership", xlab="Initialized Membership")

## Model performance by year
pred_piB$year<-parse_number(pred_piB$Year)
pred_piS$year<-parse_number(pred_piS$Year)

ggplot(pred_piB, aes(x=Truth, y=Pred))+
  geom_point(alpha=1, size=0.15,color="red") +
  facet_wrap(~year) +
  theme_bw() +
  ylab("Estimate") +
  xlab("True mixed-membership")
ggplot(pred_piS, aes(x=Truth, y=Pred))+
  geom_point(alpha=1, size=0.15,color="red") +
  facet_wrap(~year) +
  theme_bw() +
  ylab("Estimate") +
  xlab("True mixed-membership")

# initialization performance by year 
pred_piB$year<-parse_number(pred_piB$Year)
pred_piS$year<-parse_number(pred_piS$Year)

mm_init<-res_dynbi$mm_init
init_b<-c(mm_init[[2]][1,],mm_init[[2]][2,])
init_s<-c(mm_init[[1]][1,],mm_init[[1]][2,])

pred_piB2<-pred_piB%>%mutate(init_b=init_b)
pred_piS2<-pred_piS%>%mutate(init_b=init_s)

ggplot(pred_piB2, aes(x=Truth, y=init_b))+
  geom_point(alpha=1, size=0.2,color="red") +
  facet_wrap(~year) +
  theme_bw() +
  ylab("Initialization with realignment") +
  xlab("True mixed-membership")

ggplot(pred_piS2, aes(x=Truth, y=init_s))+
  geom_point(alpha=1, size=0.2,color="red") +
  facet_wrap(~year) +
  theme_bw() +
  ylab("Initialization with realignment") +
  xlab("True mixed-membership")


#badyear_S<-c(2,3,16,17,23,26,27,29,30,36,40,47,49)
#bad<-pred_piS%>%filter(year%in%badyear_S)
#good<-pred_piS%>%filter(!year%in%badyear_S)

#plot(x = bad[,5], y = bad[,6], xlab="True Value", ylab="Switched Predicted Membership")
#plot(x = good[,5], y = good[,8], xlab="True Value", ylab="Switched Predicted Membership")

#near0<-pred_piS%>%filter(Truth>0.01 & Truth<0.99)
#plot(x = near0[,5], y = near0[,8], xlab="True Value", ylab="Predicted Membership")

#d<-pred_piS%>%filter(year==7)
#plot(x = d[,5], y = d[,8], xlab="True Value", ylab="Predicted Membership")



## Plot simulated network
Z<-as.data.frame(Z)
Z2<-Z%>%mutate(groupS=z,groupB=w)
Z2$link<-Y
Z3<-Z2%>%select(groupS,groupB,link,year)
d11<-sum(Z3%>%filter(groupS==1,groupB==1)%>%select(link))/nrow(Z3%>%filter(groupS==1,groupB==1)%>%select(link))
d12<-sum(Z3%>%filter(groupS==1,groupB==2)%>%select(link))/nrow(Z3%>%filter(groupS==1,groupB==2)%>%select(link))
d21<-sum(Z3%>%filter(groupS==2,groupB==1)%>%select(link))/nrow(Z3%>%filter(groupS==2,groupB==1)%>%select(link))
d22<-sum(Z3%>%filter(groupS==2,groupB==2)%>%select(link))/nrow(Z3%>%filter(groupS==2,groupB==2)%>%select(link))
m<-matrix(c(d11,d12,d21,d22),
          nrow=2,ncol=2,
          byrow=T)
Z3$groupS<-as.factor(Z3$groupS)
Z3$groupB<-as.factor(Z3$groupB)
Z3$d<-ifelse(Z3$groupS==1 & Z3$groupB==1,d11,
             ifelse(Z3$groupS==1 & Z3$groupB==2,d12,
                    ifelse(Z3$groupS==2 &Z3$groupB==1,d21,d22)))
Z3<-cbind(Z3,Z2[,1:2])
#ggplot(Z3, aes(x = groupS, y = groupB, fill = link)) +
#  geom_tile() +
#  scale_fill_gradient(low = "white", high = "red") +
#  geom_text(aes(label = d), color = "black", size = 4)+
#  labs(x = "Senator", y = "Bill", title = "Network Generation Heatmap")

Z3$link<-as.factor(Z3$link)


ggplot(Z3, aes(x = S, y = B, color = link)) +
  geom_point(size=0.5)+ scale_colour_grey(start = 1, end = 0)+
  facet_grid(groupS ~ groupB) +
  labs(x = "Bills", y = "Senators", title = "Network Generation Heatmap")
#state 1
Z3_1<-Z3%>%filter(year<=25)
ggplot(Z3_1, aes(x = S, y = B, color = link)) +
  geom_point(size=0.5)+ scale_colour_grey(start = 1, end = 0)+
  facet_grid(groupS ~ groupB) +
  labs(x = "Bills", y = "Senators", title = "Network Generation Heatmap")
#state 2
Z3_2<-Z3%>%filter(year<=25)
ggplot(Z3_2, aes(x = S, y = B, color = link)) +
  geom_point(size=0.5)+ scale_colour_grey(start = 1, end = 0)+
  facet_grid(groupS ~ groupB) +
  labs(x = "Bills", y = "Senators", title = "Network Generation Heatmap")

# mixed membership in the simulated network
plot1<-ggplot(data = pred_piB2, aes(x = Truth)) +
     geom_density(fill = "blue", alpha = 0.5) +
       labs(title = "Bills")+theme_bw()+xlab("Membership in Group 1")
plot2<-ggplot(data = pred_piS2, aes(x = Truth)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Senators")+theme_bw()+xlab("Membership in Group 1")
grid.arrange(plot1, plot2, ncol = 2)

## Or: plot using data before running the model ##
piB2<-as.data.frame(netSim$piB)
piS2<-as.data.frame(netSim$piS)
piB2$year<- rep(1:TIME, each = nnB)
piS2$year<- rep(1:TIME, each = nnS)

plot1<-ggplot(data = piB2, aes(x = V1)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Bills")+theme_bw()+xlab("Membership in Group 1")
plot2<-ggplot(data = piS2, aes(x = V1)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Senators")+theme_bw()+xlab("Membership in Group 1")
grid.arrange(plot1, plot2, ncol = 2)

#State 1
plot1<-ggplot(data = piB2[piB2$year<=25,], aes(x = V1)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Bills")+theme_bw()+xlab("Membership in Group 1")
plot2<-ggplot(data = piS2[piS2$year<=25,], aes(x = V1)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Senators")+theme_bw()+xlab("Membership in Group 1")
grid.arrange(plot1, plot2, ncol = 2)
#State 2
plot1<-ggplot(data = piB2[piB2$year>25,], aes(x = V1)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Bills")+theme_bw()+xlab("Membership in Group 1")
plot2<-ggplot(data = piS2[piS2$year>25,], aes(x = V1)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Senators")+theme_bw()+xlab("Membership in Group 1")
grid.arrange(plot1, plot2, ncol = 2)

## What is happening in year x?
yeardata<-function(year){
  y<-year
Z2<-Z%>%mutate(groupS=z,groupB=w)
Z2$link<-Y
Z2<-Z2%>%filter(year==y)
Z3<-Z2%>%select(groupS,groupB,link)
d11<-sum(Z3%>%filter(groupS==1,groupB==1)%>%select(link))/nrow(Z3%>%filter(groupS==1,groupB==1)%>%select(link))
d12<-sum(Z3%>%filter(groupS==1,groupB==2)%>%select(link))/nrow(Z3%>%filter(groupS==1,groupB==2)%>%select(link))
d21<-sum(Z3%>%filter(groupS==2,groupB==1)%>%select(link))/nrow(Z3%>%filter(groupS==2,groupB==1)%>%select(link))
d22<-sum(Z3%>%filter(groupS==2,groupB==2)%>%select(link))/nrow(Z3%>%filter(groupS==2,groupB==2)%>%select(link))
m<-matrix(c(d11,d12,d21,d22),
          nrow=2,ncol=2,
          byrow=T)
Z3$groupS<-as.factor(Z3$groupS)
Z3$groupB<-as.factor(Z3$groupB)
Z3$d<-ifelse(Z3$groupS==1 & Z3$groupB==1,d11,
             ifelse(Z3$groupS==1 & Z3$groupB==2,d12,
                    ifelse(Z3$groupS==2 &Z3$groupB==1,d21,d22)))
Z3<-cbind(Z3,Z2[,1:2])
Z3$link<-as.factor(Z3$link)
p1<-ggplot(Z3, aes(x = S, y = B, color = link)) +
  geom_point(size=0.5)+ scale_colour_grey(start = 1, end = 0)+
  facet_grid(groupS ~ groupB) +theme_bw()+
  labs(x = "Bill", y = "Senator", title = "Network Generation Heatmap")
return(p1)
}
Z3$link2<-as.numeric(Z3$link)-1
probabilities <- Z3 %>%
  group_by(groupS, groupB) %>%
  summarize(probability_of_link = mean(link2)) %>%
  ungroup()

yeardens<-function(year){
  y<-year
  pred_piB2<-pred_piB2%>%filter(year==y)
  pred_piS2<-pred_piS2%>%filter(year==y)
  plot1<-ggplot(data = pred_piB2, aes(x = Truth)) +
    geom_density(fill = "blue", alpha = 0.5) +ylim(0,1.3)+
    labs(title = "Bills")+theme_bw()+xlab("Membership in Group 1")
  plot2<-ggplot(data = pred_piS2, aes(x = Truth)) +
    geom_density(fill = "blue", alpha = 0.5) +ylim(0,1.3)+
    labs(title = "Senators")+theme_bw()+xlab("Membership in Group 1")
  p<-grid.arrange(plot1, plot2, ncol = 2)
  return(p)
}

# Plot the wrong years
#p<-yeardens(1)
#ggsave(p,file="D:/Aprojects/NetMix/plots/1220/y1_dens.pdf",device="pdf",width = 9, height = 7)
#p<-yeardata(1)
##ggsave(p,file="D:/Aprojects/NetMix/plots/1220/y1_data.pdf",device="pdf",width = 9, height = 7)
# Compare with good years
#p<-yeardens(13)
#ggsave(p,file="D:/Aprojects/NetMix/plots/1220/y13_dens.pdf",device="pdf",width = 9, height = 7)
#p<-yeardata(13)
#ggsave(p,file="D:/Aprojects/NetMix/plots/1220/y13_data.pdf",device="pdf",width = 9, height = 7)


#ggplot(pred_piB2, aes(x=Truth))+
#  geom_density(fill = "blue", alpha = 0.5) +ylim(0,1.3)+
#  facet_wrap(~year) +
#  theme_bw() +
#  xlab("Membership in Group 1")
#ggplot(pred_piS2, aes(x=Truth))+
#  geom_density(fill = "blue", alpha = 0.5) +ylim(0,1.3)+
#  facet_wrap(~year) +
#  theme_bw() +
#  xlab("Membership in Group 1")

# niter and lowerbound in the wrong years
res_dynbi$init_niter[c(5,7,10,18,19,22)] #  9  4 84 17  4 17
res_dynbi$init_niter 
# [1] 192  39  28  48   9  29   4  47  41  84  84 111  45  37  24  55  18  17   4  37  70  17  14  52  63  90  52  66  71
#[30]  91 182  51  56  65 122 148  46  78 134  68 139  84  64  82  54  87  90 179  80  49
mean(res_dynbi$init_niter[c(5,7,10,18,19,22)]) # 22.5
mean(res_dynbi$init_niter) #67.94

res_dynbi$init_lb
#[1] -18.33485 -21.54925 -58.05289 -21.48659 -64.76354 -23.25325 -64.24366 -19.60838 -56.29238 -59.02733 -18.15181
#[12] -18.80010 -20.47846 -20.09121 -20.70268 -18.70800 -61.78846 -66.06576 -59.63365 -21.94468 -20.74039 -60.38903
#[23] -47.20506 -18.91685 -20.54324 -37.25913 -56.53977 -81.92180 -54.80411 -50.90592 -47.85753 -37.40558 -55.64426
#[34] -46.54527 -54.50066 -50.71628 -65.66365 -53.44938 -53.49331 -60.89337 -54.35163 -39.20375 -38.47908 -52.28368
#[45] -65.25825 -48.56214 -76.79909 -54.86837 -76.13405 -55.85707
res_dynbi$init_lb[c(5,7,10,18,19,22)] # -64.76354 -64.24366 -59.02733 -66.06576 -59.63365 -60.38903
mean(res_dynbi$init_lb[c(5,7,10,18,19,22)]) #-62.35383
mean(res_dynbi$init_lb) #-45.40337



#Year 3
i<-6
dy<-netSim[["df_dyad_1"]]%>%filter(year==i)
sdf<-netSim[["df_monad_S"]]%>%filter(year==i)
bdf<-netSim[["df_monad_B"]]%>%filter(year==i)
m_3<-mmsbm(formula.dyad = Y~var1,
           formula.monad = list(~VarS1, ~VarB1),
           # timeID="year",
           senderID = "id1",
           receiverID = "id2",
           nodeID = list("id","id"),
           bipartite= TRUE,
           data.dyad = dy,
           data.monad = list(sdf,bdf),
           n.blocks = c(2,2), 
           mmsbm.control = list(verbose = TRUE,
                                threads=1,
                                svi = TRUE,
                                vi_iter = 5000,
                                batch_size = 1.0,
                                conv_tol = 1e-3,
                                var_beta=list(c(0.01),
                                              c(0.01)),
                                hessian = FALSE,seed=1891
                                ))
bm3<-plogis(m_3$BlockModel)

## Year 1
i<-8
dy<-netSim[["df_dyad_1"]]%>%filter(year==i)
sdf<-netSim[["df_monad_S"]]%>%filter(year==i)
bdf<-netSim[["df_monad_B"]]%>%filter(year==i)
m_1<-mmsbm(formula.dyad = Y~1,
           #formula.monad = list(~VarS1, ~VarB1),
           # timeID="year",
           senderID = "id1",
           receiverID = "id2",
           nodeID = list("id","id"),
           bipartite= TRUE,
           data.dyad = dy,
           data.monad = list(sdf,bdf),
           n.blocks = c(2,2), 
           mmsbm.control = list(verbose = TRUE,
                                threads=1,
                                svi = TRUE,
                                vi_iter = 5000,
                                batch_size = 1.0,
                                conv_tol = 1e-2,
                                var_beta=list(c(0.01),
                                              c(0.01)),
                                hessian = FALSE,seed=3374 
           ))

bm1<-plogis(m_1$BlockModel)

netSim$piB<-netSim$piB[1:50,1]
netSim$piS<-netSim$piS[1:100,1]

plot(netSim$piB,m_3$MixedMembership2[1,])

## statistics for S
loss.mat.phi<- m_s$MixedMembership1 %*% netSim$piS
phi_ord <- clue::solve_LSAP(t(loss.mat.phi), TRUE)
orig_mm_names <- rownames(m_s$MixedMembership1)
m_s$MixedMembership1 <- m_s$MixedMembership1[phi_ord, ]
rownames(m_s$MixedMembership1) <- orig_mm_names
orig_g_names <- dimnames(m_s$BlockModel)
m_s$BlockModel <- m_s$BlockModel[phi_ord, phi_ord]
dimnames(m_s$BlockModel) <- orig_g_names #comment out following lines if no covariates
loss.mat.kappa <- m_s$Kappa %*% model.matrix(~-1+as.factor(netSim$sVec))
kappa_ord <- clue::solve_LSAP((loss.mat.kappa), TRUE)
orig_coef_names <- dimnames(m_s$MonadCoef1)
m_s$MonadCoef1 <- m_s$MonadCoef1[, phi_ord, kappa_ord, drop = FALSE]
dimnames(m_s$MonadCoef1) <- orig_coef_names
orig_k_names <- rownames(m_s$Kappa)
m_s$Kappa <- m_s$Kappa[kappa_ord, , drop = FALSE]
rownames(m_s$Kappa) <- orig_k_names
orig_tk_names <- dimnames(m_s$TransitionKernel)
m_s$TransitionKernel <- m_s$TransitionKernel[kappa_ord, kappa_ord, drop = FALSE]
dimnames(m_s$TransitionKernel) <- orig_tk_names
pred_piS <- data.frame(#netSimwork = rep(c(type), each=with(SynthnetSims[[i]],NODE*BLK*TIME)),
  NodeType = "Senator",
  Group = factor(rep(c(1:netSim$nbS), each = with(netSim,nnS*TIME))),
  Year = rep(paste("Year",1:netSim$TIME), each=netSim$nnS, times=netSim$nbS),
  State = rep(c("State 1","State 2"),c(25*netSim$nnS,25*netSim$nnS)),
  V2 = rep(netSim$df_monad_S[,"VarS1"], netSim$nbS), #comment out if no covariates
  Truth = c(netSim$piS),
  Pred = c(t(m_s$MixedMembership1[phi_ord,])))

#pred_piS <-pred_piS%>%
#  mutate(switched_Pred=1-Pred)%>%
#  mutate(best_Pred=ifelse((switched_Pred-Truth)^2<= (Pred-Truth)^2, switched_Pred, Pred ))

## statistics for B
loss.mat.phi<- m_s$MixedMembership2 %*% netSim$piB
phi_ord <- clue::solve_LSAP(t(loss.mat.phi), TRUE)
orig_mm_names <- rownames(m_s$MixedMembership2)
m_s$MixedMembership2 <- m_s$MixedMembership2[phi_ord, ]
rownames(m_s$MixedMembership2) <- orig_mm_names
orig_g_names <- dimnames(m_s$BlockModel)
m_s$BlockModel <- m_s$BlockModel[phi_ord, phi_ord]
dimnames(m_s$BlockModel) <- orig_g_names #comment out following lines if no covariates
loss.mat.kappa <- m_s$Kappa %*% model.matrix(~-1+as.factor(netSim$sVec))
kappa_ord <- clue::solve_LSAP((loss.mat.kappa), TRUE)
orig_coef_names <- dimnames(m_s$MonadCoef2)
m_s$MonadCoef2 <- m_s$MonadCoef2[, phi_ord, kappa_ord, drop = FALSE]
dimnames(m_s$MonadCoef2) <- orig_coef_names
orig_k_names <- rownames(m_s$Kappa)
m_s$Kappa <- m_s$Kappa[kappa_ord, , drop = FALSE]
rownames(m_s$Kappa) <- orig_k_names
orig_tk_names <- dimnames(m_s$TransitionKernel)
m_s$TransitionKernel <- m_s$TransitionKernel[kappa_ord, kappa_ord, drop = FALSE]
dimnames(m_s$TransitionKernel) <- orig_tk_names
pred_piB <- data.frame(#netSimwork = rep(c(type), each=with(SynthnetSims[[i]],NODE*BLK*TIME)),
  NodeType = "Bill",
  Group = factor(rep(c(1:netSim$nbB), each = with(netSim,nnB*TIME))),
  Year = rep(paste("Year",1:netSim$TIME), each=netSim$nnB, times=netSim$nbB),
  State = rep(c("State 1","State 2"),c(25*netSim$nnB,25*netSim$nnB)),
  V2 = rep(netSim$df_monad_B[,"VarB1"], netSim$nbB), #comment out if no covariates
  Truth = c(netSim$piB),
  Pred = c(t(m_s$MixedMembership2[phi_ord,])))

plot(x = pred_piB$Truth, y = pred_piB$Pred, xlab="True Value", ylab="Predicted Membership")
plot(x = pred_piS$Truth, y = pred_piS$Pred, xlab="True Value", ylab="Predicted Membership")


#the model did not swap for year 3
bm1<-plogis(m_1$BlockModel)
calculate_norm <- function(matrix1, matrix2) {
        return(base::norm(matrix1 - matrix2, type = "f"))
       }
calculate_norm(bm1,bm3) # 0.1950383
calculate_norm(bm1,bm3[c(2, 1), c(2, 1)]) #0.290135
bm1
bm3
bm3[c(2, 1), c(2, 1)]
# data in year 3:
y<-8
Z2<-Z%>%mutate(groupS=z,groupB=w)
Z2$link<-Y
Z2<-Z2%>%filter(year==y)
Z3<-Z2%>%select(groupS,groupB,link)
d11<-sum(Z3%>%filter(groupS==1,groupB==1)%>%select(link))/nrow(Z3%>%filter(groupS==1,groupB==1)%>%select(link))
d12<-sum(Z3%>%filter(groupS==1,groupB==2)%>%select(link))/nrow(Z3%>%filter(groupS==1,groupB==2)%>%select(link))
d21<-sum(Z3%>%filter(groupS==2,groupB==1)%>%select(link))/nrow(Z3%>%filter(groupS==2,groupB==1)%>%select(link))
d22<-sum(Z3%>%filter(groupS==2,groupB==2)%>%select(link))/nrow(Z3%>%filter(groupS==2,groupB==2)%>%select(link))
m<-matrix(c(d11,d12,d21,d22),
          nrow=2,ncol=2,
          byrow=T)
Z3$groupS<-as.factor(Z3$groupS)
Z3$groupB<-as.factor(Z3$groupB)
Z3$d<-ifelse(Z3$groupS==1 & Z3$groupB==1,d11,
             ifelse(Z3$groupS==1 & Z3$groupB==2,d12,
                    ifelse(Z3$groupS==2 &Z3$groupB==1,d21,d22)))
Z3<-cbind(Z3,Z2[,1:2])
Z3$link<-as.factor(Z3$link)

Z3$link2<-as.numeric(Z3$link)-1
probabilities <- Z3 %>%
  group_by(groupS, groupB) %>%
  summarize(probability_of_link = mean(link2)) %>%
  ungroup()
matrix(c(probabilities$probability_of_link),nrow=2,ncol=2,byrow=T)





## statistics for S
loss.mat.phi<- m_s$MixedMembership1 %*% netSim$piS
phi_ord <- clue::solve_LSAP(t(loss.mat.phi), TRUE)
orig_mm_names <- rownames(m_s$MixedMembership1)
m_s$MixedMembership1 <- m_s$MixedMembership1[phi_ord, ]
rownames(m_s$MixedMembership1) <- orig_mm_names
orig_g_names <- dimnames(m_s$BlockModel)
m_s$BlockModel <- m_s$BlockModel[phi_ord, phi_ord]
dimnames(m_s$BlockModel) <- orig_g_names #comment out following lines if no covariates
loss.mat.kappa <- m_s$Kappa %*% model.matrix(~-1+as.factor(netSim$sVec))
kappa_ord <- clue::solve_LSAP((loss.mat.kappa), TRUE)
orig_coef_names <- dimnames(m_s$MonadCoef1)
m_s$MonadCoef1 <- m_s$MonadCoef1[, phi_ord, kappa_ord, drop = FALSE]
dimnames(m_s$MonadCoef1) <- orig_coef_names
orig_k_names <- rownames(m_s$Kappa)
m_s$Kappa <- m_s$Kappa[kappa_ord, , drop = FALSE]
rownames(m_s$Kappa) <- orig_k_names
orig_tk_names <- dimnames(m_s$TransitionKernel)
m_s$TransitionKernel <- m_s$TransitionKernel[kappa_ord, kappa_ord, drop = FALSE]
dimnames(m_s$TransitionKernel) <- orig_tk_names
pred_piS <- data.frame(#netSimwork = rep(c(type), each=with(SynthnetSims[[i]],NODE*BLK*TIME)),
  NodeType = "Senator",
  Group = factor(rep(c(1:netSim$nbS), each = with(netSim,nnS*TIME))),
  Year = rep(paste("Year",1:netSim$TIME), each=netSim$nnS, times=netSim$nbS),
  State = rep(c("State 1","State 2"),c(25*netSim$nnS,25*netSim$nnS)),
  V2 = rep(netSim$df_monad_S[,"VarS1"], netSim$nbS), #comment out if no covariates
  Truth = c(netSim$piS),
  Pred = c(t(m_s$MixedMembership1[phi_ord,])))

#pred_piS <-pred_piS%>%
#  mutate(switched_Pred=1-Pred)%>%
#  mutate(best_Pred=ifelse((switched_Pred-Truth)^2<= (Pred-Truth)^2, switched_Pred, Pred ))

## statistics for B
loss.mat.phi<- m_s$MixedMembership2 %*% netSim$piB
phi_ord <- clue::solve_LSAP(t(loss.mat.phi), TRUE)
orig_mm_names <- rownames(m_s$MixedMembership2)
m_s$MixedMembership2 <- m_s$MixedMembership2[phi_ord, ]
rownames(m_s$MixedMembership2) <- orig_mm_names
orig_g_names <- dimnames(m_s$BlockModel)
m_s$BlockModel <- m_s$BlockModel[phi_ord, phi_ord]
dimnames(m_s$BlockModel) <- orig_g_names #comment out following lines if no covariates
loss.mat.kappa <- m_s$Kappa %*% model.matrix(~-1+as.factor(netSim$sVec))
kappa_ord <- clue::solve_LSAP((loss.mat.kappa), TRUE)
orig_coef_names <- dimnames(m_s$MonadCoef2)
m_s$MonadCoef2 <- m_s$MonadCoef2[, phi_ord, kappa_ord, drop = FALSE]
dimnames(m_s$MonadCoef2) <- orig_coef_names
orig_k_names <- rownames(m_s$Kappa)
m_s$Kappa <- m_s$Kappa[kappa_ord, , drop = FALSE]
rownames(m_s$Kappa) <- orig_k_names
orig_tk_names <- dimnames(m_s$TransitionKernel)
m_s$TransitionKernel <- m_s$TransitionKernel[kappa_ord, kappa_ord, drop = FALSE]
dimnames(m_s$TransitionKernel) <- orig_tk_names
pred_piB <- data.frame(#netSimwork = rep(c(type), each=with(SynthnetSims[[i]],NODE*BLK*TIME)),
  NodeType = "Bill",
  Group = factor(rep(c(1:netSim$nbB), each = with(netSim,nnB*TIME))),
  Year = rep(paste("Year",1:netSim$TIME), each=netSim$nnB, times=netSim$nbB),
  State = rep(c("State 1","State 2"),c(25*netSim$nnB,25*netSim$nnB)),
  V2 = rep(netSim$df_monad_B[,"VarB1"], netSim$nbB), #comment out if no covariates
  Truth = c(netSim$piB),
  Pred = c(t(m_s$MixedMembership2[phi_ord,])))

plot(x = pred_piB$Truth, y = pred_piB$Pred, xlab="True Value", ylab="Predicted Membership")
plot(x = pred_piS$Truth, y = pred_piS$Pred, xlab="True Value", ylab="Predicted Membership")




