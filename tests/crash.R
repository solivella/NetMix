library(NetMix)
load("crash.Rdata") 

set.seed(999)
test <- mmsbm(formula.dyad = MID_onset_DY ~ IGOmems_joint + 
                           ally + contiguity + dist + peaceyrs + 
                           spline1 + spline2 + spline3,
                         formula.monad =  ~ polity + logNMC, 
                         data.dyad = MID_dyad1820, data.monad = MID_monad1820,
                         senderID = "country1", receiverID = "country2", 
                         nodeID = "country", 
                         timeID = "year",
                         n.blocks = 4, 
                         n.hmmstates = 1,
                         directed=FALSE,
                         mmsbm.control = list(
                           verbose = TRUE,
                           em_iter = 3000,
                           threads = parallel::detectCores()
                         ))
