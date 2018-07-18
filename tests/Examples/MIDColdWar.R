#########################
## Example: MID network
## during the Cold War
#########################

library(NetMix)

data(MID_monad_CW)
data(MID_dyad_CW)


## Estimate model (takes about 5 minutes)
system.time(fit_onset <- mmsbm(formula.dyad = mid_onset ~ hially + dist + 
                                 peaceyrs + spline1 + spline2 + spline3,
                               formula.monad =  ~ polity + logcinc,
                               data.dyad = MID_dyad_CW, data.monad = MID_monad_CW,
                               senderID = "country1", receiverID = "country2", 
                               nodeID = "country", 
                               timeID = "year",
                               n.groups = 3, 
                               n.hmmstates = 2,
                               directed = FALSE,
                               mmsbm.control = list(em_iter = 3000)))

## Summarize model
summary(fit_onset) 

## Plot blockmodel
plot(fit_onset)





