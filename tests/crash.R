load("tests/MIDdata.Rdata") 
library(NetMix)
set.seed(999)
MID_monad$ccode <- as.factor(as.numeric(paste(MID_monad$ccode)))
MID_dyad$ccode1 <- as.factor(as.numeric(paste(MID_dyad$ccode1)))
MID_dyad$ccode2 <- as.factor(as.numeric(paste(MID_dyad$ccode2)))
system.time(fit_onset_f <- mmsbm(formula.dyad = MID_onset ~ trade_dep_low + IGOmems_joint +
                       ally + contiguity + dist +
                       peaceyrs + spline1 + spline2 + spline3,
                     formula.monad =  ~ polity + logNMC + gdp_growth_67_a + region,
                     data.dyad = MID_dyad, 
                     data.monad = MID_monad,
                     senderID = "country1", 
                     receiverID = "country2", 
                     nodeID = "country", 
                     timeID = "year",
                     n.groups = 4, 
                     n.hmmstates = 2,
                     directed=FALSE,
                     mmsbm.control = list(var_b = c(1,1),
                                          var_beta = 1,
                                          var_xi = 1,
                                          var_gamma = 1,
                                          #mu_b = c(-5, 5),
                                          verbose = TRUE,
                                          em_iter = 3000,
                                          threads = 4
                     )))
