library(NetMix)
load("crash.Rdata") 


set.seed(999)
MID_onset_4_2 <- mmsbm(formula.dyad = MID_onset_DY ~ IGOmems_joint + 
                         ally + contiguity + dist + peaceyrs + 
                         spline1 + spline2 + spline3,
                       formula.monad =  ~ polity + logNMC, 
                       data.dyad = MID_dyad2, data.monad = MID_monad2,
                       senderID = "country1", receiverID = "country2", 
                       nodeID = "country", 
                       timeID = "year",
                       n.blocks = 4, 
                       n.hmmstates = 2,
                       directed=FALSE,
                       mmsbm.control = list(verbose = TRUE,
                                            em_iter = 500,
                                            conv_tol = 1e-4,
                                            threads = parallel::detectCores()
                       ))
