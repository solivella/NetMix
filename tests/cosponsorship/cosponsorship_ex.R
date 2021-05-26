# Cosponsorship : focus on 107th US Congress -- January 3, 2001 – January 3, 2003
rm(list=ls())
library(NetMix)
library(stringr)
#set up the data
#NOTES
#** #Dean Barkley (id=40106) missing ideology information; hard-code/ remove
#SH
SH<-read.csv(file="tests/cosponsorship/cosponsor_data/SH.csv")
#dates
dates<-read.csv(file="tests/cosponsorship/cosponsor_data/dates.txt",header=FALSE,col.names=c("bill_date"))
dates$bill_date<-as.Date(dates$bill_date,format=c("%Y-%m-%d"))
#ICPSR code of each bill sponsor
sponsors<-read.table(file="tests/cosponsorship/cosponsor_data/sponsors.txt",header=FALSE,col.names="sponsor")
cosponsors<-read.table(file="tests/cosponsorship/cosponsor_data/cosponsors.txt",header=FALSE,sep=""
                       ,col.names=paste0("V",seq_len(399)), fill = TRUE)
orig.bills<-read.table(file="tests/cosponsorship/cosponsor_data/bills.txt",header=FALSE,col.names="bill_name")
pvtbills<-read.table(file="tests/cosponsorship/cosponsor_data/pvtbills.txt",header=FALSE)

#combine data
data<-data.frame(date=dates$bill_date,billname=orig.bills$bill_name,private=pvtbills$V1,sponsor=sponsors$sponsor)
data<-cbind(data,cosponsors)
data.107<-subset(data,date>as.Date("2001-01-03") & date < as.Date("2003-01-03"))
saveRDS(data.107,file="tests/cosponsorship/cosponsor_data/data_107.rds")
#reduce data.107 to just bills from senate: SC, SE, SJ, SN, (remove SP amendements don't have as many sponsors, but still do have them)
ind<-c(which(str_detect(data.107$billname, "SC_")),which(str_detect(data.107$billname, "SE_")),which(str_detect(data.107$billname, "SJ_")),which(str_detect(data.107$billname, "SN_")))
data.107<-data.107[ind,]
#n1: unique senators (sponsors + cosponsors) from data.107 = result in 100 after remove dean barkley (no nominate scores because didn't have min roll call)
data.107<-data.107[which(data.107$sponsor!=40106&data.107$V1!=40106),]
legis.names=unique(sort(unlist(apply(data.107[,4:ncol(data.107)],2,unique))))
saveRDS(legis.names,file="tests/cosponsorship/cosponsor_data/legisnames.rds") #names associated with node1: 100
#n2: unique bills = 2633
node2=nrow(data.107)

# data.107 node1 names must match monadic data1 names
#monadic data 1: senators; party, node1, time
senate107<-subset(SH,congress==107&chamber=="S"&ids!=40106)
senate107<-with(senate107, senate107[order(labels),])#order by senator name alpha

senate107<-senate107[,c("ids","labels","seniority","ideol1","ideol2","party")]
#code gender from list:https://history.house.gov/Exhibitions-and-Publications/WIC/Historical-Data/Women-Representatives-and-Senators-by-Congress/
senate107$sex <- ifelse(senate107$labels=="Boxer, Barbara [CA]"|senate107$labels=="Cantwell, Maria [WA]"|
                          senate107$labels=="Carnahan, Jean [MO]"|senate107$labels=="Clinton, Hillary Rodham [NY]"|
                          senate107$labels=="Collins, Susan M. [ME]"|senate107$labels=="Feinstein, Dianne [CA]"|
                          senate107$labels=="Hutchison, Kay Bailey [TX]"|senate107$labels=="Landrieu, Mary [LA]"|
                          senate107$labels=="Lincoln, Blanche [AR]"|senate107$labels=="Mikulski, Barbara A. [MD]"|
                          senate107$labels=="Murray, Patty [WA]"|senate107$labels=="Snowe, Olympia J. [ME]"|
                          senate107$labels=="Stabenow, Debbie [MI]","f","m")
sum(is.na(match(senate107$ids,legis.names)))#check: sum(is.na(match(legis.names,congress$id)))
senate107$time<-1 #add time variable
node1.names<-senate107$ids
node1=length(node1.names)
senate107$node1<-1:nrow(senate107)

#monadic data 2: bills; private, type, node2, time
billnum<-billtype<-billsponsor<-NA
for(i in 1:length(data.107$billname)){
  tmp<-as.character(data.107$billname[i])
  billtype[i]<-str_split(tmp,"_")[[1]][1]
  billnum[i]<-str_split(tmp,".txt")[[1]][1]
  billsponsor[i]<-data.107$sponsor[i]
}
bills<-data.frame(private=data.107$private,billsponsor=billsponsor,billnum=billnum,billtype=billtype,node2=1:node2,time=1)
#pick up covariates of sponsor of that bill -- billsponsor_seniority, billsponsor_ideol1, billsponsor_ideol2, billsponsor_party,billsponsor_sex
ind<-match(bills$billsponsor,senate107$ids) #index in senate data
bills$billsponsor_seniority<-senate107$seniority[ind]
bills$billsponsor_ideol1<-senate107$ideol1[ind]
bills$billsponsor_ideol2<-senate107$ideol2[ind]
bills$billsponsor_party<-senate107$party[ind]
bills$billsponsor_sex<-senate107$sex[ind]

#picking up bills covariates: congressional bills project for bill categories
#Bill type ("hr" (house bill); "s" (senate bill); "hres" (House resolution); "sres" (Senate resolution); 
#"hcon" (House Concurrent Resolution); "scon" (Senate Concurrent Resolution); "hjres" (House Joint Resolution); 
#"sjres" (Senate Joint Resolution).
bills93_114<-read.csv(file="tests/cosponsorship/cosponsor_data/bills93-114.csv",header=TRUE,sep=";")
bills93_114$billtype<-NA
bills93_114$billtype[bills93_114$BillType=="hconres"]<-"HC"
bills93_114$billtype[bills93_114$BillType=="hres"]<-"HE"
bills93_114$billtype[bills93_114$BillType=="hjres"]<-"HJ"
bills93_114$billtype[bills93_114$BillType=="hr"]<-"HR"
#no amendments in congressional bills project
bills93_114$billtype[bills93_114$BillType=="scon"]<-"SC"
bills93_114$billtype[bills93_114$BillType=="sres"]<-"SE"
bills93_114$billtype[bills93_114$BillType=="sjres"]<-"SJ"
bills93_114$billtype[bills93_114$BillType=="s"]<-"SN"
saveRDS(bills93_114,file="tests/cosponsorship/cosponsor_data/bills93-114.rds")
bills107<-subset(bills93_114,Cong==107&Chamber=="1")#just Senators Congress 107
bills107$BillID_lastnum<-NA #last number to use for matching later
for(i in 1:nrow(bills107)){
  bills107$BillID_lastnum[i]<-str_split(bills107$BillID[i],"-")[[1]][3]
}
saveRDS(bills107,file="tests/cosponsorship/cosponsor_data/cov-bills107.rds")
#match with `bills` data via billtype, BillNum, and pull in information on "Major" for category
bills$Major<-bills$Minor<-NA #a lot of the Major are NA -- nearly all because of the "SP" amendments
for(i in 1:nrow(bills)){
  print(i)
  ind<-str_split(bills$billnum[i],"_")[[1]][3]
  tmp<-subset(bills107,BillID_lastnum==ind&billtype==bills$billtype[i])
  # if(nrow(tmp)>1) {print("STOP\n")
  #   break}
  if(nrow(tmp)>0){
    bills$Major[i]<-tmp$Major[1]
    bills$Minor[i]<-tmp$Minor[1]
  }else{
    next
  }
}
#code cruder ver Major 
#1=Economy (1=Macro, 4=Agriculture, 5=Labor,8=Energy,15=Domestic Commerce,17=Tech,18=Foreign Trade,19=International), 
#2=Legal (2=Civil Rights,9=Immigration,12=Law&Crime)
#3=Social programs/Public goods (3=Health,6=Education,7=Environment, 10=Transportation,13=Social Welfare,14=Housing, 21=Public lands,23=Culture)
#4=Security (16=Defense)
#5=Gov operations (20=Gov operations)
#6=Other (all others)
bills$Major2<-NA
bills$Major2<-ifelse(bills$Major==1|bills$Major==4|bills$Major==5|bills$Major==8|bills$Major==15|bills$Major==17|bills$Major==18|bills$Major==19,1,bills$Major2)
bills$Major2<-ifelse(bills$Major==2|bills$Major==9|bills$Major==12,2,bills$Major2)
bills$Major2<-ifelse(bills$Major==3|bills$Major==6|bills$Major==7|bills$Major==10|bills$Major==13|bills$Major==14|bills$Major==21|bills$Major==23,3,bills$Major2)
bills$Major2<-ifelse(bills$Major==16,4,bills$Major2)
bills$Major2<-ifelse(bills$Major==20,5,bills$Major2) #gov operations=5, other=6
bills$Major2[is.na(bills$Major)]<-6
#code cruder ver Major (ver 2)
#1=Economy (1=Macro, 4=Agriculture, 5=Labor,8=Energy,15=Domestic Commerce,17=Tech)
#2=International (18=Foreign Trade,19=International)
#3=Social programs/Public goods (2=Civil Rights,3=Health,6=Education,7=Environment,9=Immigration,10=Transportation,
      #13=Social Welfare,14=Housing, 21=Public lands,23=Culture)
#4=Security (12=Law&Crime,16=Defense)
#5=Gov operations (20=Gov operations)
#6=Other (all others)
bills$Major3<-NA
bills$Major3<-ifelse(bills$Major==1|bills$Major==4|bills$Major==5|bills$Major==8|bills$Major==15|bills$Major==17,1,bills$Major3)#economy
bills$Major3<-ifelse(bills$Major==18|bills$Major==19,2,bills$Major3)#international
bills$Major3<-ifelse(bills$Major==2|bills$Major==3|bills$Major==6|bills$Major==7|bills$Major==9|bills$Major==10|
                       bills$Major==13|bills$Major==14|bills$Major==21|bills$Major==23,3,bills$Major3)#social programs/pub goods
bills$Major3<-ifelse(bills$Major==12|bills$Major==16,4,bills$Major3)#defense
bills$Major3<-ifelse(bills$Major==20,5,bills$Major3)#gov operations
bills$Major3[is.na(bills$Major)]<-6#other
#Do some level fixing
senate107$party<-ifelse(senate107$party==100,"Democrat",ifelse(senate107$party==200,"Republican",NA)) #100=Democrat, 200=Republican
senate107$party<-as.factor(senate107$party)
senate107$sex<-as.factor(senate107$sex)
bills$Major2<-as.factor(bills$Major2)
bills$Major3<-as.factor(bills$Major3)
bills$billsponsor_party<-ifelse(bills$billsponsor_party==100,"Democrat",ifelse(bills$billsponsor_party==200,"Republican",NA)) #100=Democrat, 200=Republican
bills$billsponsor_party<-as.factor(bills$billsponsor_party)
bills$billsponsor_sex<-as.factor(bills$billsponsor_sex)

#save monadic data
saveRDS(node1.names,file="tests/cosponsorship/cosponsor_data/node1names.rds") #names associated with node1: 1-100
saveRDS(senate107,"tests/cosponsorship/cosponsor_data/senate-107.rds")
saveRDS(bills,"tests/cosponsorship/cosponsor_data/bills-107.rds")


#check that all node1.names are ONLY the senators who cosponsor
cosponsors<-unique((unlist(apply(data.107[,5:ncol(data.107)],2,unique))))
cosponsors<-cosponsors[which(!is.na(cosponsors))]
sum(is.na(match(node1.names,cosponsors)))
sum(is.na(match(cosponsors,node1.names)))

#dyadic data:
  #go through bill by bill, filling in sponsor info; let sponsor variable = 1 if sponsor (use this to remove rows after)
  #code cosponsors as Y=1, otherwise Y=0; 
dyad.data<-vector("list",length=node2)
for(i in 1:node2){
  cat("Bill no:  ",i,"\n")
  #Sponsor variable
  tmp.sponsor<-rep(0,node1)
  tmp.sponsor[match(data.107$sponsor[i],node1.names)]<-1
  #Sponsor seniority; removed - this is now @ monadic level
  #tmp.seniority<-rep(senate107$seniority[match(data.107$sponsor[i],senate107$ids)],node1)
  #Sponsor ideology1,2; removed - this is now @ monadic level
  #tmp.ideol1<-rep(senate107$ideol1[match(data.107$sponsor[i],senate107$ids)],node1)
  #tmp.ideol2<-rep(senate107$ideol2[match(data.107$sponsor[i],senate107$ids)],node1)
  #Sponsor party; removed - this is now @ monadic level
  #tmp.party<-rep(senate107$party[match(data.107$sponsor[i],senate107$ids)],node1)
  #Sponsor sex; removed - this is now @ monadic level
  #tmp.sex<-rep(senate107$sex[match(data.107$sponsor[i],senate107$ids)],node1)
  #Y Edge variable
  y<-rep(0,node1)
  tmp.vec<-data.107[i,5:ncol(data.107)][!is.na(data.107[i,5:ncol(data.107)])] #only cosponsors
  tmp.vec<-tmp.vec[which(tmp.vec!=40100)]
  y[match(tmp.vec,node1.names)]<-1
  dyad.data[[i]]<-data.frame(date=data.107$date[i],node1=1:node1, node2=i, Y=y, sponsor=tmp.sponsor
                             #, sponsor.seniority=tmp.seniority,sponsor.ideol1=tmp.ideol1
                             #,sponsor.ideol2=tmp.ideol2,sponsor.party=tmp.party,sponsor.sex=tmp.sex
                             )
}
dyad.data<-do.call(rbind,dyad.data)
dyad.data$time<-1 #only 107th session
#rm rows where the node has sponsored the bill (any rows with sponsor=1)
saveRDS(dyad.data,file="tests/cosponsorship/cosponsor_data/dyad-data-107-FULL.rds")
dyad.data<-subset(dyad.data,sponsor==0)

#save dyadic data
saveRDS(dyad.data,file="tests/cosponsorship/cosponsor_data/dyad-data-107.rds")

########### FINISH ##############

### Call in data ###
dyad.data<-readRDS("tests/cosponsorship/cosponsor_data/dyad-data-107.rds")
congress<-readRDS("tests/cosponsorship/cosponsor_data/congress-107.rds")
bills<-readRDS("tests/cosponsorship/cosponsor_data/bills-107.rds")
node1=nrow(congress)
node2=nrow(bills)

library(NetMix)

### Model 1: k1=k2=2, monad1.pred=monad2.pred=dyad.pred=1
co.model1 <- mmsbm(formula.dyad = Y ~ sponsor,
                    formula.monad1 = ~ party, 
                    formula.monad2 = ~ private,
                    senderID = "node1",
                    receiverID = "node2",
                    nodeID1 = "node1",
                    nodeID2 = "node2",
                    timeID = "time",
                    data.dyad = dyad.data,
                    data.monad1 = congress,
                    data.monad2 = bills,
                    n.blocks1 = 2,#blocks for legislators
                    n.blocks2 = 2,#blocks for bills
                    n.hmmstates = 1,
                    directed = TRUE,
                    nodes2 = node2,
                    npred2 = 2,#+1 for intercept
                    mmsbm.control = list(mu_b = c(3, -5),
                                         var_b = c(1, 1),
                                         var_beta = 1,
                                         eta = 1.3,
                                         seed = 1234,
                                         alpha = 0.5,
                                         spectral = TRUE,
                                         verbose=TRUE,
                                         em_iter = 1000,
                                         conv_tol = 1e-4,
                                         bipartite = TRUE
                                         ,permute = TRUE
                                         ,directed = TRUE
                    ))

save(co.model1,file="tests/cosponsorship/co.model1.RData")

### Model 2: k1=3,k2=2, monad1.pred=2, monad2.pred=dyad.pred=1
co.model2 <- mmsbm(formula.dyad = Y ~ sponsor,
                   formula.monad1 = ~ party + senate, 
                   formula.monad2 = ~ private,
                   senderID = "node1",
                   receiverID = "node2",
                   nodeID1 = "node1",
                   nodeID2 = "node2",
                   timeID = "time",
                   data.dyad = dyad.data,
                   data.monad1 = congress,
                   data.monad2 = bills,
                   n.blocks1 = 3,#blocks for legislators
                   n.blocks2 = 2,#blocks for bills
                   n.hmmstates = 1,
                   directed = TRUE,
                   nodes2 = node2,
                   npred2 = 2,#+1 for intercept
                   mmsbm.control = list(mu_b = c(3, -5),
                                        var_b = c(1, 1),
                                        var_beta = 1,
                                        eta = 1.3,
                                        seed = 1234,
                                        alpha = 0.5,
                                        spectral = TRUE,
                                        verbose=TRUE,
                                        em_iter = 1000,
                                        conv_tol = 1e-4,
                                        bipartite = TRUE
                                        ,permute = TRUE
                                        ,directed = TRUE
                   ))

save(co.model2,file="tests/cosponsorship/co.model2.RData")

###################################################################################################
################################### Plots #########################################################
###################################################################################################
## Plots of relevant parameters
# Libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

node1.names<-readRDS("tests/cosponsorship/cosponsor_data/node1names.rds")
congress<-readRDS("tests/cosponsorship/cosponsor_data/congress-107.rds")
bills<-readRDS("tests/cosponsorship/cosponsor_data/bills-107.rds")

###### Model 1 ########
load("tests/cosponsorship/co.model1.RData")
#Blockmodel:
co.model1$BlockModel

#Dyadic predictor: sponsor (rather than just cosponsor)
co.model1$DyadCoef

#Congress monadic predictor: party
SH<-read.csv(file="tests/cosponsorship/cosponsor_data/SH.csv") #congress info
SH.107<-subset(SH,congress==107)
co.model1$MonadCoef1

SH$labels[match(congress$id,SH$ids)] #congressperson name

#Bill monadic predictor: private (or not private bill); applies to individual or entity instead of to everyone
co.model1$MonadCoef2

###  Congressional latent groups
dat.m1<-data.frame(name=rep(SH$labels[match(congress$id,SH$ids)],2)
                   ,prob=c(co.model1$`MixedMembership 1`[1,],co.model1$`MixedMembership 1`[2,])
                   ,group=rep(c(1,2),each=length(SH$labels[match(congress$id,SH$ids)])))
dat.m1<-data.frame(name=SH.107$labels[match(congress$id,SH.107$ids)],party=SH.107$party[match(congress$id,SH.107$ids)]
                   ,group1=co.model1$`MixedMembership 1`[1,],group2=co.model1$`MixedMembership 1`[2,])
dat.m1$party<-ifelse(dat.m1$party==100,"Democrat",ifelse(dat.m1$party==200,"Republican","Other"))
  #sort by group1 probas
dat.m1<-dat.m1[order(-dat.m1$group1),]
dat.m1$name<-as.vector(dat.m1$name)
#sample from rows of congresspeople because way too many
tmp<-dat.m1[seq(1,nrow(dat.m1),3),]
tmp$name<-factor(tmp$name, levels = tmp$name)
axiscolors<-ifelse(tmp$party=="Democrat","royalblue1",ifelse(tmp$party=="Republican","firebrick3","gray"))
plot(1,type='n', xlim=c(0,2*nrow(tmp)), ylim=c(0,1), ylab="Probability in Group 1", xaxt='n', xlab="",las=2)
axis(1, at=seq(1,2*nrow(tmp),2),labels=tmp$name,las=2,cex.axis=0.5)
points(seq(1,2*nrow(tmp),2),tmp$group1,pch=16,cex=0.8
     ,col=axiscolors #color them based on party
     ,las=2,xlab="")
lines(seq(1,2*nrow(tmp),2),tmp$group1,col="darkgray")

### Bills latent groups
## Recall:
#HC   House Concurrent Resolutions
#HE 	House Resolutions
#HJ 	House Joint Resolutions
#HR 	House Bills
#HZ 	House Amendments
#SC 	Senate Concurrent Resolutions
#SE 	Senate Resolutions
#SJ 	Senate Joint Resolutions
#SN 	Senate Bills
#SP 	Senate Amendments 
data.107<-readRDS(data.107,file="tests/cosponsorship/cosponsor_data/data_107.rds")
dat.m2<-data.frame(name=data.107$billname,private=bills$private
                   ,group1=co.model1$`MixedMembership 2`[1,],group2=co.model1$`MixedMembership 2`[2,])
#sort by group1 probas
dat.m2<-dat.m2[order(-dat.m2$group1),]
dat.m2$name<-as.vector(dat.m2$name)
#sample from rows of bills because way too many
tmp<-dat.m2[seq(1,nrow(dat.m2),50),]
tmp$name<-factor(tmp$name, levels = tmp$name)
axiscolors<-ifelse(tmp$private==1,"darkgoldenrod1","gray")
plot(1,type='n', xlim=c(0,2*nrow(tmp)), ylim=c(0,1), ylab="Probability in Group 1", xaxt='n', xlab="",las=2)
axis(1, at=seq(1,2*nrow(tmp),2),labels=tmp$name,las=2,cex.axis=0.5)
points(seq(1,2*nrow(tmp),2),tmp$group1,pch=16,cex=0.8
       ,col=axiscolors #color them based on private bill or not
       ,las=2,xlab="")
lines(seq(1,2*nrow(tmp),2),tmp$group1,col="darkgray")

###### Model 2 ########
load("tests/cosponsorship/co.model2.RData")
#Blockmodel:
co.model2$BlockModel

#Dyadic predictor: sponsor (rather than just cosponsor)
co.model2$DyadCoef

#Congress monadic predictor: party; senate
#npred=2+1, blk=2, state=1
SH<-read.csv(file="tests/cosponsorship/cosponsor_data/SH.csv") #congress info
SH.107<-subset(SH,congress==107)
co.model2$MonadCoef1

SH$labels[match(congress$id,SH$ids)] #congressperson name

#Bill monadic predictor: private (or not private bill)
co.model2$MonadCoef2 #npred=1+1, blk=2, state=1

###  Congressional latent groups
dat.m1<-data.frame(name=rep(SH$labels[match(congress$id,SH$ids)],2)
                   ,prob=c(co.model2$`MixedMembership 1`[1,],co.model2$`MixedMembership 1`[2,])
                   ,group=rep(c(1,2),each=length(SH$labels[match(congress$id,SH$ids)])))
dat.m1<-data.frame(name=SH.107$labels[match(congress$id,SH.107$ids)],party=SH.107$party[match(congress$id,SH.107$ids)]
                   ,group1=co.model2$`MixedMembership 1`[1,],group2=co.model2$`MixedMembership 1`[2,]
                   ,group3=co.model2$`MixedMembership 1`[3,])
dat.m1$party<-ifelse(dat.m1$party==100,"Democrat",ifelse(dat.m1$party==200,"Republican","Other"))
#sort by group1 probas
dat.m1<-dat.m1[order(-dat.m1$group1),]
dat.m1$name<-as.vector(dat.m1$name)
#sample from rows of congresspeople because way too many
tmp<-dat.m1[seq(1,nrow(dat.m1),1),]
tmp$name<-factor(tmp$name, levels = tmp$name)
axiscolors<-ifelse(tmp$party=="Democrat","royalblue1",ifelse(tmp$party=="Republican","firebrick3","gray"))
plot(1,type='n', xlim=c(0,2*nrow(tmp)), ylim=c(0,1), ylab="Probability in Group 1", xaxt='n', xlab="",las=2)
axis(1, at=seq(1,2*nrow(tmp),2),labels=tmp$name,las=2,cex.axis=0.5)
#points(seq(1,2*nrow(tmp),2),tmp$group1,pch=16,cex=0.8
       #,col=axiscolors #color them based on party
       #,las=2,xlab="")
lines(seq(1,2*nrow(tmp),2),tmp$group1,col="darkgray")

#points(seq(1,2*nrow(tmp),2),tmp$group2,pch=16,cex=0.8
       #,col=axiscolors #color them based on party
       #,las=2,xlab="")
lines(seq(1,2*nrow(tmp),2),tmp$group2,col="lightgoldenrod1")

#points(seq(1,2*nrow(tmp),2),tmp$group3,pch=16,cex=0.8
       #,col=axiscolors #color them based on party
       #,las=2,xlab="")
lines(seq(1,2*nrow(tmp),2),tmp$group3,col="goldenrod1")

#Ternary plot
library(Ternary)

TernaryPlot(atip="Group 1", btip="Group 2", ctip="Group 3", point="up")

AddToTernary(points, tmp[,c("group1","group2","group3")]*100, pch=21, cex=1.5, 
             bg=apply(tmp[,c("group1","group2","group3")]*100, 1,
                      function (x) rgb(x[1], x[2], x[3], 128, maxColorValue=255)
                      #,character(1)
                      )
)

tmp2<-apply(tmp[,c("group1","group2","group3")]*100,1,jitter,amount=10)
AddToTernary(text, t(tmp2), tmp$name, cex=0.5, font=2)


### Bills latent groups
## Recall:
#HC   House Concurrent Resolutions
#HE 	House Resolutions
#HJ 	House Joint Resolutions
#HR 	House Bills
#HZ 	House Amendments
#SC 	Senate Concurrent Resolutions
#SE 	Senate Resolutions
#SJ 	Senate Joint Resolutions
#SN 	Senate Bills
#SP 	Senate Amendments 
data.107<-readRDS(file="tests/cosponsorship/cosponsor_data/data_107.rds")
dat.m2<-data.frame(name=data.107$billname,private=bills$private
                   ,group1=co.model2$`MixedMembership 2`[1,],group2=co.model2$`MixedMembership 2`[2,])
#sort by group1 probas
dat.m2<-dat.m2[order(-dat.m2$group1),]
dat.m2$name<-as.vector(dat.m2$name)
#sample from rows of bills because way too many
tmp<-dat.m2[seq(1,nrow(dat.m2),50),]
tmp$name<-factor(tmp$name, levels = tmp$name)
axiscolors<-ifelse(tmp$private==1,"darkgoldenrod1","gray")
plot(1,type='n', xlim=c(0,2*nrow(tmp)), ylim=c(0,1), ylab="Probability in Group 1", xaxt='n', xlab="",las=2)
axis(1, at=seq(1,2*nrow(tmp),2),labels=tmp$name,las=2,cex.axis=0.5)
points(seq(1,2*nrow(tmp),2),tmp$group1,pch=16,cex=0.8
       ,col=axiscolors #color them based on private bill or not
       ,las=2,xlab="")
lines(seq(1,2*nrow(tmp),2),tmp$group1,col="darkgray")


