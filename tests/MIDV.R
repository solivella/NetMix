########################################
## MID application for dynMMSBM paper ##
########################################


library(NetMix)
load("tests/MIDdata.Rdata") 

MID_dyad2 <- MID_dyad[MID_dyad$year <= 2010,]
MID_monad2 <- MID_monad[MID_monad$year <= 2010,]

start_time <- Sys.time()
set.seed(999)
MID_onset_6_2o <- mmsbm(formula.dyad = MID_onset_DY ~ IGOmems_joint + 
                          ally + contiguity + dist + peaceyrs + 
                          spline1 + spline2 + spline3,
                        formula.monad =  ~ polity + logNMC, 
                        data.dyad = MID_dyad2, data.monad = MID_monad2,
                        senderID = "country1", receiverID = "country2", 
                        nodeID = "country", 
                        timeID = "year",
                        n.blocks = 6, 
                        n.hmmstates = 2,
                        directed=FALSE,
                        mmsbm.control = list(verbose = TRUE,
                                             em_iter = 500,
                                             mu_b = c(-5, 5),
                                             conv_tol = 1e-4,
                                             threads = parallel::detectCores())
                        )
end_time <- Sys.time()
end_time - start_time

load("tests/fit_onset.Rdata")
fit_onset <- MID_onset_6_2o

fit_onset$iter

summary(fit_onset)

plot(fit_onset)
plot(fit_onset, type="membership")

mean(fit_onset$dyadic.data$MID_onset_DY)
mean(predict(fit_onset, type="response"))
mean(predict(fit_onset, type="response", parametric_mm=TRUE))



## plot some individual states over time
states <- c("USA", "Iceland", "UK", "Nicaragua", "Russia", "Iraq")
plots <- list()
ind <- 1
for(i in states){
  Nodes <- unlist(lapply(strsplit(colnames(fit_onset$MixedMembership), "@"), "[[", 1))
  avgmem <- lapply(1:nrow(fit_onset$MixedMembership), function(x){
    fit_onset$MixedMembership[x,which(Nodes==i)]
  })
  ts <- unlist(lapply(strsplit(names(avgmem[[1]]), "@"), "[[", 2))
  avgmem <- as.data.frame(cbind(rep(ts, nrow(fit_onset$MixedMembership)), unlist(avgmem),
                                rep(1:nrow(fit_onset$MixedMembership), each=length(ts))))
  colnames(avgmem) <- c("Time", "Membership", "Group")
  avgmem$Group <- factor(avgmem$Group, levels=length(unique(avgmem$Group)):1)
  if(class(avgmem$Membership) == "factor"){avgmem$Membership <- as.numeric(as.character(avgmem$Membership))}
  if(class(avgmem$Time) == "factor"){avgmem$Time <- as.numeric(as.character(avgmem$Time))}
  plots[[ind]] <- ggplot() + ggtitle(paste("Group Membership Over Time,", i)) + theme(plot.title = element_text(hjust = 0.5)) +
    geom_area(aes(y = Membership, x = Time, fill=Group), data = avgmem,
              stat="identity", position="stack")  + guides(fill=guide_legend(title="Group"))
  ind <- ind + 1
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
pdf("cluster_time_node_onset.pdf", height=8, width=15)
multiplot(plots[[1]], plots[[2]], plots[[3]], plots[[4]],
          plots[[5]], plots[[6]], cols=3)
dev.off()



## look at covariate effects
p <- covFX(fit_onset, "polity", sd(fit_onset$monadic.data$polity)) 
p[[1]]
head(sort(p[[3]])) # countries made more peaceful by pos shift in polity
head(sort(p[[3]], decreasing=T)) # countries made more violent 
hist(p[[4]])

plot(fit_onset, type="effect", FX=p)



## compare to logit - in sample fit

load("tests/MIDdata.Rdata")
MID_dyad2 <- MID_dyad

MID_dyad2$polity_missing <- ifelse((is.na(MID_dyad2$polity_1) | is.na(MID_dyad2$polity_2)), 1, 0)
MID_dyad2$polity_1[is.na(MID_dyad2$polity_1)] <- 0
MID_dyad2$polity_2[is.na(MID_dyad2$polity_2)] <- 0

MID_dyad2$IGOmems_joint_missing <- ifelse(is.na(MID_dyad2$IGOmems_joint), 1, 0)
MID_dyad2$IGOmems_joint[is.na(MID_dyad2$IGOmems_joint)] <- 0

MID_dyad2$dist_missing <- ifelse(is.na(MID_dyad2$dist), 1, 0)
MID_dyad2$dist[is.na(MID_dyad2$dist)] <- 0


MID_dyad2$dem <- ifelse(MID_dyad2$polity_2 > 6 & MID_dyad2$polity_1 > 6, 1, 0)
MID_dyad2$aut <- ifelse(MID_dyad2$polity_2 < -6 & MID_dyad2$polity_1 < -6, 1, 0)
MID_dyad2$mixed <- ifelse(MID_dyad2$dem==0 & MID_dyad2$aut==0, 1, 0)
MID_dyad2$cinc_ratio <- MID_dyad2$logNMC_1 / MID_dyad2$logNMC_2
MID_dyad2b <- MID_dyad2[MID_dyad2$year <= 2010,]


logit.MID <- glm(MID_onset_DY ~ IGOmems_joint + ally + contiguity + dist + 
                   dem + aut + cinc_ratio + 
                   IGOmems_joint_missing + dist_missing +  
                   polity_missing + 
                   peaceyrs + spline1 + spline2 + spline3,
                 data = MID_dyad2b, family = binomial(link = "logit"))
pred.logit <- predict(logit.MID, newdata = MID_dyad2b, type="response")
pred.MMSB <- predict(fit_onset, type="response")
pred.MMSBp <- predict(fit_onset, type="response", parametric_mm=TRUE)

library(ROCR)
MID.roc <- prediction(predictions=c(pred.MMSB), labels=fit_onset$dyadic.data$MID_onset_DY)
MID.roc2 <- performance(MID.roc, measure="tpr", x.measure="fpr")

MID.rocD <- prediction(predictions=c(pred.MMSBp), labels=fit_onset$dyadic.data$MID_onset_DY)
MID.roc2D <- performance(MID.rocD, measure="tpr", x.measure="fpr")

logit.roc <- prediction(predictions=c(pred.logit), labels=MID_dyad2b$MID_onset_DY)
logit.roc2 <- performance(logit.roc, measure="tpr", x.measure="fpr")


pdf("ROC_onset.pdf", width=7, height=7)
par(mfrow=c(1,1))
plot(unlist(attributes(MID.roc2)["x.values"]), 
     unlist(attributes(MID.roc2)["y.values"]),
     type="l",
     xlab="False Positive Rate",
     ylab="True Positive Rate", lwd=2)
lines(unlist(attributes(logit.roc2)["x.values"]),
      unlist(attributes(logit.roc2)["y.values"]),
      col="red", lwd=2)
lines(unlist(attributes(MID.roc2D)["x.values"]),
      unlist(attributes(MID.roc2D)["y.values"]),
      col="dark green", lwd=2)
lines(c(0, 1), c(0,1), lty=2)
legend(x=0.7, y=0.2, legend=c("hmmsb", "hmmsb w/ parametric prediction", "logit"),
       col=c("black", "dark green", "red"), lty=1, lwd=2)
dev.off()
 

# compare to logit - out of sample fit

MID.sub <- MID_dyad[MID_dyad$year <= 2008,]
MID.monad.sub <- MID_monad[MID_monad$year <= 2008,]

set.seed(999)
MID_onset_2008 <- mmsbm(formula.dyad = MID_onset_DY ~ IGOmems_joint + 
                          ally + contiguity + dist,
                        formula.monad =  ~ polity + logNMC, 
                        data.dyad = MID.sub, data.monad = MID.monad.sub,
                        senderID = "country1", receiverID = "country2", 
                        nodeID = "country", 
                        timeID = "year",
                        n.blocks = 6, 
                        n.hmmstates = 2,
                        directed=FALSE,
                        mmsbm.control = list(verbose = TRUE,
                                             em_iter = 500,
                                             mu_b = c(-5, 5),
                                             conv_tol = 1e-4,
                                             threads = parallel::detectCores())
)


MID_dyad2c <- MID_dyad2b[MID_dyad2b$year <= 2008,]
logit.2008 <- glm(MID_onset_DY ~ IGOmems_joint + ally + contiguity + dist + 
                   dem + aut + cinc_ratio + 
                   IGOmems_joint_missing + dist_missing +  
                   polity_missing,
                 data = MID_dyad2b, family = binomial(link = "logit"))


MID.new <- MID_dyad2b[MID_dyad2b$year %in% 2009:2010,]
MID.monad.new <- MID_monad[as.numeric(paste(MID_monad$year)) %in% 2009:2010,]
MID.monad.new$polity_missing <- ifelse(is.na(MID.monad.new$polity), 1, 0)
MID.monad.new$polity[is.na(MID.monad.new$polity)] <- 0

MID.monad.new2 <- MID_monad[MID_monad$year == 2008,]
MID.monad.new2$polity_missing <- ifelse(is.na(MID.monad.new2$polity), 1, 0)
MID.monad.new2$polity[is.na(MID.monad.new2$polity)] <- 0
MID_dyad.new2 <- MID_dyad2b[MID_dyad2b$year==2008,]

MMSB.predictOOS <- predict(MID_onset_2008, new.data.dyad=MID_dyad.new2, 
                           new.data.monad=MID.monad.new2, type="response",
                           parametric_mm = TRUE, forecast=TRUE)
logit.predictOOS <- predict(logit.2008, newdata=MID.new, type="response")

mean(predict(MID_onset_2008, new.data.dyad=MID_dyad.new2, new.data.monad=MID.monad.new2, type="response"))
mean(predict(MID_onset_2008, new.data.dyad=MID_dyad.new2, new.data.monad=MID.monad.new2, type="response", parametric_mm=TRUE))
mean(predict(MID_onset_2008, new.data.dyad=MID_dyad.new2, new.data.monad=MID.monad.new2, type="response", parametric_mm=TRUE, forecast=TRUE))



library(ROCR)
MMSB.roc <- prediction(predictions=c(MMSB.predictOOS), labels=MID.new$MID_onset_DY)
MMSB.roc2 <- performance(MMSB.roc, measure="tpr", x.measure="fpr")

logit.roc <- prediction(predictions=c(logit.predictOOS), labels=MID.new$MID_onset_DY)
logit.roc2 <- performance(logit.roc, measure="tpr", x.measure="fpr")

pdf("ROC_onsetOOS.pdf", height=8, width=8)
par(mfrow=c(1,1))
plot(unlist(attributes(MMSB.roc2)["x.values"]), 
     unlist(attributes(MMSB.roc2)["y.values"]),
     type="l",
     xlab="False Positive Rate",
     ylab="True Positive Rate", lwd=2)
lines(unlist(attributes(logit.roc2)["x.values"]),
      unlist(attributes(logit.roc2)["y.values"]),
      type="l",
      col="red", lwd=2)
lines(c(0, 1), c(0,1), lty=2)
legend(x=0.7, y=0.2, legend=c("hmmsb","logit"),
       col=c("black", "red"), lty=1, lwd=2)
dev.off()




