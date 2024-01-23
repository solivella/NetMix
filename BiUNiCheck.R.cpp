#########################
## Bipartite simulations
## Two monadic predictors
## No dyadic predictors.
########################

##Clear worskapce
rm(list=ls())

library(NetMix)
library(MCMCpack)
library(igraph)
library(ggplot2)
library(doMC)
library(tidyverse)



## Net Creator
net_maker <- function(XS, XB, beta, BM){
  nnS <- nrow(XS)
  nnB <- nrow(XB)
  nbS <- nrow(BM)
  nbB <- ncol(BM)
  
  alphaS <- exp(XS %*% beta)
  alphaB <- exp(XB %*% beta)
  
  piS <- t(apply(alphaS, 1, function(x)MCMCpack::rdirichlet(1, x)))
  piB <- t(apply(alphaB, 1, function(x)MCMCpack::rdirichlet(1, x)))
  
  ##Dyadic 
  Z1 <- rnorm(nnS*nnB)
  ##Network
  y <- rbinom(nnS*nnB, 1, plogis(c(piS %*% BM %*% t(piB)))) 
  
  ##Dyadic
  df_dyad_1 <- data.frame(Y = y,
                          var1 = Z1,  
                          id1 = paste0("S",rep(1:nnS, times = nnB)),
                          id2 = paste0("B", rep(1:nnB, each = nnS)))
  ##Monadic
  df_monad_S <- data.frame(varS1 = XS[,2],
                           id = paste0("S", 1:nnS))
  df_monad_B <- data.frame(varB1 = XB[,2],
                           id = paste0("B", 1:nnB))
  return(list(df_dyad_1 = df_dyad_1,
              df_monad_S = df_monad_S,
              df_monad_B = df_monad_B, 
              nbS = nbS, nbB = nbB, nnS = nnS, nnB = nnB,
              piS = piS, piB = piB ))
}
## Full sim funtion
sim_fun <- function(XS, XB, beta, BM, SetSize){
  net <- net_maker(XS, XB, beta, BM)
  
  ## Estimate bipartite MMSBM
  system.time(res_bipart <- mmsbm(formula.dyad = Y~1,
                                  formula.monad = list(~varS1, ~varB1),
                                  senderID = "id1",
                                  receiverID = "id2",
                                  nodeID = list("id","id"),
                                  bipartite= TRUE,
                                  data.dyad = net[["df_dyad_1"]],
                                  data.monad = list(net[["df_monad_S"]],net[["df_monad_B"]]),
                                  n.blocks = c(net[["nbS"]],net[["nbB"]]),
                                  mmsbm.control = list(verbose = FALSE,
                                                       threads=1,
                                                       svi = FALSE,
                                                       vi_iter = 1000,
                                                       batch_size = 1.0,
                                                       conv_tol = 1e-3,
                                                       var_beta=list(c(0.01),
                                                                     c(0.01)),
                                                       hessian = TRUE)))
  
  
  ## Plot results
  #plot(res_bipart$LowerBound_full, type='l')
  
  ##Type S
  loss.mat.netS<- res_bipart$MixedMembership1 %*% net[["piS"]]
  netS_order <- clue::solve_LSAP(t(loss.mat.netS), TRUE)
  pred_piS <- data.frame(NodeType = "Senator",
                         Size = SetSize,
                         Group = factor(rep(c(1:net[["nbS"]]), each = net[["nnS"]], times = 1)),
                         Truth = c(c(net[["piS"]])),
                         Pred = c(c(t(res_bipart$MixedMembership1[netS_order,]))))
  #corS <- cor(pred_piS[,c("Truth","Pred")])[1,2]
  #Type B
  loss.mat.netB<- res_bipart$MixedMembership2 %*% net[["piB"]]
  netB_order <- clue::solve_LSAP(t(loss.mat.netB), TRUE)
  pred_piB <- data.frame(NodeType = "Bill",
                         Size = SetSize,
                         Group = factor(rep(c(1:net[["nbB"]]), each = net[["nnB"]], times = 1)),
                         Truth = c(c(net[["piB"]])),
                         Pred = c(c(t(res_bipart$MixedMembership2[netB_order,]))))
  pred_pi <- rbind(pred_piS, pred_piB)
  #corB <- cor(pred_piB[,c("Truth","Pred")])[1,2]
  
  
  ## Monad coef
  
  
  ##SEs
  #S
  se1 <- sqrt(diag(res_bipart$vcov_monad1))
  
  #B
  se2 <- sqrt(diag(res_bipart$vcov_monad2))
  
  
  return(list(Preds = pred_pi,
              model = res_bipart
  ))
  
}

##All cases
asim_nnB <- c(200, 2000)
asim_nnS <- floor(asim_nnB/2)

all_n <- rbind(cbind(nB=asim_nnB,nS=asim_nnS))
all_X <- lapply(seq.int(nrow(all_n)),
                function(i){
                  set.seed(831213)
                  nnB_this <- all_n[i,1]
                  nnS_this <- all_n[i,2]
                  list(XS = cbind(1, rnorm(nnS_this, 0, 1.5)),
                       XB = cbind(1, rnorm(nnB_this, 0, 1.5)))
                })

## Scenarios
beta_list <- list(beta_easy = array(c(-4.5, -4.5, ##Intercepts
                                      0.0, 0.0), ## Predictor coefficients
                                    c(2, 2)),
                  beta_med = array(c(0.05, 0.75,
                                     -0.75,  -1.0),
                                   c(2, 2)),
                  beta_hard = array(c(0.0, 0.0,
                                      -0.75,  -1.0),
                                    c(2, 2)))
BM_list <- list(BM_easy = matrix(qlogis(c(0.85, 0.01, 0.01, 0.99)), ncol = 2),
                BM_med = matrix(qlogis(c(0.65, 0.2, 0.35, 0.75)), ncol = 2),
                BM_hard = matrix(qlogis(c(0.65, 0.5, 0.4, 0.45)), ncol = 2))
## Indeces
cases_data <- rep(seq.int(nrow(all_n)), each = 25*3) 
cases_scen <- rep(seq.int(3), times = 25*2) 
sizes <- rep(c("Small","Large"), each=25*3)

## Sims
##Run sim 
library(doParallel)
registerDoMC(10)
set.seed(831213)
res_sim <- foreach(i=1:length(cases_data), .errorhandling='remove',.combine="rbind",
                   .export=c("all_X", "cases_data", "cases_scen", "BM_list","beta_list")) %dopar% {
                     sim_fun(all_X[[cases_data[i]]][["XS"]],
                             all_X[[cases_data[i]]][["XB"]], 
                             beta_list[[cases_scen[i]]],
                             BM_list[[cases_scen[i]]], 
                             sizes[i])
                   }

save(res_sim, file=paste0("tests/BipartSim.RData"))




# 
corr_res <- do.call(rbind, res_sim[c(76:78)])
corr_res$Scenario <- factor(rep(c("Easy","Medium","Hard"), each=6000), levels=c("Easy","Medium","Hard"))
# ## Plots
ggplot(filter(corr_res,Group == 1),aes(x=Truth)) +
  facet_grid(NodeType~Scenario)+
  geom_density(fill="gray60") +
  ylab("") +
  xlab("Membership in Group 1") +
  theme_bw()

# ##Corr
## Mod-list 
pdf("CorrPlotSim.pdf", height=4, width=7)
ggplot(corr_res, aes(x=Truth, y=Pred))+
  geom_point(alpha=1, size=2.5) +
  facet_grid(NodeType~Scenario) +
  theme_bw() +
  ylab("Estimate") +
  xlab("True mixed-membership")
dev.off()

##Error of SE
res_sim_large <- do.call(rbind,lapply(res_sim[151:225],
                                      function(x){
                                        data.frame(coef_est = c(x$MonadCoef1[2,2,1],
                                                                x$MonadCoef2[2,2,1]),
                                                   SE = c(sqrt(x$vcov_monad1[4,4]),
                                                          sqrt(x$vcov_monad2[4,4])),
                                                   NodeType = c("Bills","Senators"),
                                                   Size = "Large Network")
                                      }))
res_sim_small <- do.call(rbind,lapply(res_sim[226:300],
                                      function(x){
                                        data.frame(coef_est = c(x$MonadCoef1[2,2,1],
                                                                x$MonadCoef2[2,2,1]),
                                                   SE = c(sqrt(x$vcov_monad1[4,4]),
                                                          sqrt(x$vcov_monad2[4,4])),
                                                   NodeType = c("Bills","Senators"),
                                                   Size = "Small Network")
                                      }))
res_sim_df <- rbind(res_sim_small, res_sim_large)
res_sim_df$Size <- factor(res_sim_df$Size, levels=c("Small Network","Large Network"))

se_res <- res_sim_df %>% 
  group_by(NodeType, Size) %>% 
  summarize(est_se_bias = abs(SE-sd(coef_est)))
pdf("SEErrorPlot.pdf", height=3, width=6)
ggplot(se_res, aes(as.factor(NodeType), est_se_bias)) +
  facet_wrap(~Size) +
  geom_boxplot() +
  geom_hline(yintercept=0.0) +
  xlab("Type of Node") +
  ylab("Approximate bias of standard error") +
  theme_bw()
dev.off()

## Goodness of fit
sim_net <- res_sim[[151]]
class(sim_net) <- c("mmsbmB","mmsbm")
sim_net$forms$timeID <- "(tid)"
gof_obj <- gof(sim_net,
               gof_stat = c("Degree","Dyad Shared Partners", "Geodesics"),
               level=0.9)
pdf("GOFSim.pdf", height=3, width=11)
plot(gof_obj)
dev.off()


### Time tests

timed_sims <- function(XS, XB, beta, BM){
  net <- net_maker(XS, XB, beta, BM)
  
  propS <- 0.4 #min(5e2/net[["nnS"]], 0.6) 
  propB <- 0.4 #min(1e3/net[["nnB"]], 0.6) 
  
  ## Estimate bipartite MMSBM
  strt_time <- Sys.time()
  system.time(res_bipart <- mmsbm(formula.dyad = Y~1,
                                  formula.monad = list(~varS1, ~varB1),
                                  senderID = "id1",
                                  receiverID = "id2",
                                  nodeID = list("id","id"),
                                  bipartite= TRUE,
                                  data.dyad = net[["df_dyad_1"]],
                                  data.monad = list(net[["df_monad_S"]],net[["df_monad_B"]]),
                                  n.blocks = c(net[["nbS"]],net[["nbB"]]),
                                  mmsbm.control = list(verbose = FALSE,
                                                       threads=1,
                                                       svi = TRUE,
                                                       vi_iter = 70000,
                                                       batch_size = c(propS, propB),
                                                       conv_tol = 1e-3,
                                                       var_beta=list(c(0.01),
                                                                     c(0.01)),
                                                       hessian = FALSE)))
  end_time <- Sys.time()
  return(data.frame(n = net[["nnB"]] + net[["nnS"]],
                    time = (end_time-strt_time),
                    iter = res_bipart$niter,
                    converged = res_bipart$converged))
}

asim_nnB <- seq(200, 10000, length.out = 15)
asim_nnS <- floor(asim_nnB/2)

all_n <- rbind(cbind(nB=asim_nnB,nS=asim_nnS))
all_X <- lapply(seq.int(nrow(all_n)),
                function(i){
                  set.seed(831213)
                  nnB_this <- all_n[i,1]
                  nnS_this <- all_n[i,2]
                  list(XS = cbind(1, rnorm(nnS_this, 0, 1.5)),
                       XB = cbind(1, rnorm(nnB_this, 0, 1.5)))
                })

## Scenarios
beta_list <- list(beta_med = array(c(0.05, 0.75,##Intercepts
                                     -0.75,  -1.0), ## Predictor coefficients
                                   c(2, 2)))

BM_list <- list(BM_med = matrix(qlogis(c(0.65, 0.2, 0.35, 0.75)), ncol = 2))
##Run sim 
library(doParallel)
registerDoMC(15)
set.seed(831213)
res_sim_times <- foreach(i=1:length(all_X), .errorhandling='remove',.combine="rbind",
                   .export=c("all_X","BM_list","beta_list")) %dopar% {
                     timed_sims(all_X[[i]][["XS"]],
                                all_X[[i]][["XB"]], 
                                beta_list[[1]],
                                BM_list[[1]])
                   }

save(res_sim_times, file=paste0("tests/BipartSimTimes.RData"))
res_sim_times <- res_sim_times %>% ## Adjust time scales
  mutate(time=case_when(n < 1000 ~ time/60,
                        between(n, 1000, 6600) ~ time, 
                        between(n, 7000, 9000) ~ time*60,
                        n > 9000 ~ time *3600)) #%>% 
  #filter(n>4000)

pdf("TimeSim.pdf",height=3, width=4)
ggplot(res_sim_times, aes(x=n, y=(time/iter))) +
  geom_line() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(y="Seconds per iteration", x="Total number of vertices") +
  theme_bw()
dev.off()
  
