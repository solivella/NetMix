#####################################
## Testing versions of MID data on ##
## NetMix model #####################
#####################################

library(NetMix)
load("tests/MIDdata.Rdata") 

MID_monad$country[MID_monad$country=="Sao Tom\xe9-Principe"] <- "Sao Tome Principe"
MID_dyad$country1[MID_dyad$country1=="Sao Tom\xe9-Principe"] <- "Sao Tome Principe"
MID_dyad$country2[MID_dyad$country2=="Sao Tom\xe9-Principe"] <- "Sao Tome Principe"


## Simulate some random data to see what clusters are recovered
MID_dyad$sim2 <- rbinom(nrow(MID_dyad), 1, 0.01)

set.seed(999)
MID_sim2 <- mmsbm(formula.dyad = sim2 ~ ally,# + contiguity + dist + 
                  #peaceyrs + spline1 + spline2 + spline3,
                  #formula.monad =  ~ polity + logNMC,
                  data.dyad = MID_dyad, data.monad = MID_monad,
                  senderID = "country1", receiverID = "country2", 
                  nodeID = "country", 
                  timeID = "year",
                  n.blocks = 4, 
                  n.hmmstates = 2,
                  directed=FALSE,
                  missing="listwise deletion",
                  mmsbm.control = list(verbose = TRUE,
                                       #mu_b = c(-5, 5),
                                       em_iter = 500,
                                       conv_tol = 1e-4,
                                       threads = parallel::detectCores()
                  ))
head(MID_sim2)
lapply(1:4, function(x){sample(unique(MID_monad$country), 10)})

lapply(head(MID_sim2), function(x){
  s <- names(x)
  sapply(s, function(y){nrow(MID_sim2$monadic.data[MID_sim2$monadic.data$`(nid)`==y,])})
})


MID_dyad$ally_sim <- rbinom(nrow(MID_dyad), 1, 0.01)
set.seed(999)
MID_sim2b <- mmsbm(formula.dyad = sim2 ~ ally_sim,# + contiguity + dist + 
                  #peaceyrs + spline1 + spline2 + spline3,
                  #formula.monad =  ~ polity + logNMC,
                  data.dyad = MID_dyad, data.monad = MID_monad,
                  senderID = "country1", receiverID = "country2", 
                  nodeID = "country", 
                  timeID = "year",
                  n.blocks = 4, 
                  n.hmmstates = 2,
                  directed=FALSE,
                  missing="listwise deletion",
                  mmsbm.control = list(verbose = TRUE,
                                       #mu_b = c(-5, 5),
                                       em_iter = 500,
                                       conv_tol = 1e-4,
                                       threads = parallel::detectCores()
                  ))
head(MID_sim2b)
lapply(head(MID_sim2b), function(x){
  s <- names(x)
  sapply(s, function(y){nrow(MID_sim2b$monadic.data[MID_sim2b$monadic.data$`(nid)`==y,])})
})
lapply(1:4, function(x){sample(unique(MID_monad$country), 10)})



MID_dyad3 <- MID_dyad[MID_dyad$year %in% c(1948:2010) & !is.na(MID_dyad$MID_onset_DY),]
MID_monad3 <- MID_monad[MID_monad$year  %in% c(1948:2010),]
states <- sort(table(MID_monad3$country), decreasing=T)
k <- names(states)[states ==63]
MID_dyad4 <- MID_dyad3[MID_dyad3$country1 %in% k & MID_dyad3$country2 %in% k,]
MID_monad4 <- MID_monad3[MID_monad3$country  %in% k,]

set.seed(999)
MID_sim3 <- mmsbm(formula.dyad = sim2 ~ ally,# + contiguity + dist + 
                  #peaceyrs + spline1 + spline2 + spline3,
                  #formula.monad =  ~ polity + logNMC,
                  data.dyad = MID_dyad4, data.monad = MID_monad4,
                  senderID = "country1", receiverID = "country2", 
                  nodeID = "country", 
                  timeID = "year",
                  n.blocks = 4, 
                  n.hmmstates = 2,
                  directed=FALSE,
                  missing="listwise deletion",
                  mmsbm.control = list(verbose = TRUE,
                                       #mu_b = c(-5, 5),
                                       em_iter = 500,
                                       conv_tol = 1e-4,
                                       threads = parallel::detectCores()
                  ))
head(MID_sim3)


## Estimate using dataset with no missing nodes 
MID_dyad3 <- MID_dyad[MID_dyad$year %in% c(1948:2010) & !is.na(MID_dyad$MID_onset_DY),]
MID_monad3 <- MID_monad[MID_monad$year  %in% c(1948:2010),]

states <- sort(table(MID_monad3$country), decreasing=T)
k <- names(states)[states ==63]

MID_dyad4b <- MID_dyad3[MID_dyad3$country1 %in% k & MID_dyad3$country2 %in% k,]
MID_monad4b <- MID_monad3[MID_monad3$country  %in% k,]


plot(unique(MID_monad4$year), tapply(MID_monad4$country, MID_monad4$year, length))


set.seed(999)
MID_onset_4_2 <- mmsbm(formula.dyad = MID_DY ~ IGOmems_joint + 
                         ally + contiguity + dist + peaceyrs + 
                         spline1 + spline2 + spline3,
                       formula.monad =  ~ polity + logNMC, 
                       data.dyad = MID_dyad4b, data.monad = MID_monad4b,
                       senderID = "country1", receiverID = "country2", 
                       nodeID = "country", 
                       timeID = "year",
                       n.blocks = 4, 
                       n.hmmstates = 2,
                       directed=FALSE,
                       mmsbm.control = list(verbose = TRUE,
                                            em_iter = 500,
                                            mu_b = c(-5, 5),
                                            conv_tol = 1e-4,
                                            threads = parallel::detectCores()
                       ))
cluster.time(MID_onset_4_2)
summary(MID_onset_4_2)
plot(MID_onset_4_2)
exp(MID_onset_4_2$BlockModel) / (1 + exp(MID_onset_4_2$BlockModel)) # really high
mean(MID_onset_4_2$dyadic.data$MID_DY) # edges occur in 0.009% of cases
head(MID_onset_4_2) 
head(MID_onset_4_2, t=1955:1970, n=5)



## Same but with no peace years and splines
## Try with no covariates, without splines, listwise deletion
# try exact same specification as before



set.seed(999)
MID_onset_4_2b <- mmsbm(formula.dyad = MID_DY ~ IGOmems_joint + 
                          ally + contiguity + dist,# + peaceyrs + 
                        #spline1 + spline2 + spline3,
                        formula.monad =  ~ polity + logNMC, 
                        data.dyad = MID_dyad4b, data.monad = MID_monad4b,
                        senderID = "country1", receiverID = "country2", 
                        nodeID = "country", 
                        timeID = "year",
                        n.blocks = 4, 
                        n.hmmstates = 2,
                        directed=FALSE,
                        mmsbm.control = list(verbose = TRUE,
                                             em_iter = 500,
                                             mu_b = c(-5, 5),
                                             conv_tol = 1e-4,
                                             threads = parallel::detectCores()
                        ))
cluster.time(MID_onset_4_2b)
summary(MID_onset_4_2b)
plot(MID_onset_4_2b)
exp(MID_onset_4_2b$BlockModel) / (1 + exp(MID_onset_4_2b$BlockModel)) # this looks reasonable
colMeans(exp(MID_onset_4_2b$BlockModel) / (1 + exp(MID_onset_4_2b$BlockModel)) )
head(MID_onset_4_2b)
head(MID_onset_4_2b, t=1955:1990, n=20)



