##############################################################
## Estimation of the posterior probabilities (abc function) ##
##############################################################

require(abc)

##vector with the observed sumstats (FST_F:M, FST_HG:M, FST_HG:M)
stat.obs<-c(0.06,0.07,0.17)

# ##open BayeSSC stat files
# panmixia<-read.csv("panmixia_stat.csv",header=T)
# split<-read.csv("split_stat.csv",header=T)
# 
# ##matrix with the simulated sumstats: pm (panmixia model); sm (split model)
# pm.stat.sim<-panmixia[,c(16,23,39)]
# colnames(pm.stat.sim) <- c("fst.FvM","fst.HGvM","fst.HGvF")
# 
# sm.stat.sim<-split[,c(16,23,39)]
# colnames(sm.stat.sim) <- c("fst.FvM","fst.HGvM","fst.HGvF")
# 
# ##matrix with the simulated parameters: N (size at Neolithic (T=300)); UP (size at Upper Palaeolithic (T=1800))
# pm.par.sim<-panmixia[,63:64]
# colnames(pm.par.sim)<-c("N","UP")
# 
# sm.par.sim<-split[,67:68]
# colnames(sm.par.sim)<-c("N","UP")


##open modified BayeSSC stat files (FST_F:M; FST_HG:M; FST_HG:M; size at N; size at UP)
panmixia<-read.csv("panmixia_stat_mod.csv",header=T)
split<-read.csv("split_stat_mod.csv",header=T)

##matrix with the simulated sumstats: pm (panmixia model); sm (split model)
pm.stat.sim<-panmixia[,c(1,2,3)]
colnames(pm.stat.sim) <- c("fst.FvM","fst.HGvM","fst.HGvF") #names for the sumstats

sm.stat.sim<-split[,c(1,2,3)]
colnames(sm.stat.sim) <- c("fst.FvM","fst.HGvM","fst.HGvF") #names for the sumstats

##matrix with the simulated parameters: N (size at Neolithic (T=300)); UP (size at Upper Palaeolithic (T=1800))
pm.par.sim<-panmixia[,4:5]
colnames(pm.par.sim)<-c("N","UP") #names for the parameters

sm.par.sim<-split[,4:5]
colnames(sm.par.sim)<-c("N","UP") #names for the parameters



###############################################
##rejection algorithm (Pritchard et al.,1999)##
###############################################

pm_rej <- abc(target=stat.obs, param=pm.par.sim, sumstat=pm.stat.sim,tol=.1,method="rejection")
summary(pm_rej, intvl = .9)

sm_rej <- abc(target=stat.obs, param=sm.par.sim, sumstat=sm.stat.sim,tol=.1,method="rejection")
summary(sm_rej, intvl = .9)



##########################################################################################
##local linear regression (Beaumont et al., 2002) with correction for heteroscedascity  ##
##########################################################################################

## with logit transformation [ logit(p)= log(p/(1-p))]

#Panmixia model
pm.lin.logit <- abc(target=stat.obs, param=pm.par.sim, sumstat=pm.stat.sim,tol=.01,method="loclinear",transf=c("logit","logit"),logit.bounds = rbind(c(1000,100000),c(10,5000)))
summary(pm.lin.logit, intvl = .9)
plot(pm.lin.logit, param=pm.par.sim)
# plot(pm.lin.logit, param=pm.par.sim,file="pm.lin.logit") #save to pdf

# Split model
sm.lin.logit <- abc(target=stat.obs, param=sm.par.sim, sumstat=sm.stat.sim,tol=.01,method="loclinear",transf=c("logit","logit"),logit.bounds = rbind(c(1000,100000),c(10,5000)))
summary(sm.lin.logit , intvl = .9)
plot(sm.lin.logit, param=sm.par.sim)
# plot(sm.lin.logit, param=sm.par.sim,file="sm.lin.logit") #save to pdf

## with log transformation [ log(p)]

#Panmixia model
pm.lin.log <- abc(target=stat.obs, param=pm.par.sim, sumstat=pm.stat.sim,tol=.01,method="loclinear",transf=c("log","log"))
summary(pm.lin.log, intvl = .9)
plot(pm.lin.log, param=pm.par.sim)
# plot(pm.lin.log, param=pm.par.sim,file="pm.lin.log") #save to pdf

# Split model
sm.lin.log <- abc(target=stat.obs, param=sm.par.sim, sumstat=sm.stat.sim,tol=.01,method="loclinear",transf=c("log","log"))
summary(sm.lin.log , intvl = .9)
plot(sm.lin.log, param=sm.par.sim)
# plot(sm.lin.log, param=sm.par.sim,file="sm.lin.log") #save to pdf


###################################################################################
## neural networks with correction for heteroscedascity (Blum and François,2010) ##
###################################################################################

## with logit transformation [ logit(p)= log(p/(1-p))]

#Panmixia model
pm.net.logit <- abc(target=stat.obs, param=pm.par.sim, sumstat=pm.stat.sim,tol=.01,method="neuralnet",transf=c("logit","logit"),logit.bounds = rbind(c(1000,100000),c(10,5000)))
summary(pm.net.logit, intvl = .9)
plot(pm.net.logit, param=pm.par.sim)
#plot(pm.net.logit, param=pm.par.sim,file="pm.lin.logit") #save to pdf

# Split model
sm.net.logit <- abc(target=stat.obs, param=sm.par.sim, sumstat=sm.stat.sim,tol=.01,method="neuralnet",transf=c("logit","logit"),logit.bounds = rbind(c(1000,100000),c(10,5000)))
summary(sm.net.logit , intvl = .9)
plot(sm.net.logit, param=sm.par.sim)
#plot(sm.net.logit, param=sm.par.sim,file="sm.lin.logit") #save to pdf

## with log transformation [ log(p)]

#Panmixia model
pm.net.log <- abc(target=stat.obs, param=pm.par.sim, sumstat=pm.stat.sim,tol=.01,method="neuralnet",transf=c("log","log"))
summary(pm.net.log, intvl = .9)
plot(pm.net.log, param=pm.par.sim)
#plot(pm.net.log, param=pm.par.sim,file="pm.lin.log") #save to pdf

# Split model
sm.net.log <- abc(target=stat.obs, param=sm.par.sim, sumstat=sm.stat.sim,tol=.01,method="neuralnet",transf=c("log","log"))
summary(sm.net.log , intvl = .9)
plot(sm.net.log, param=sm.par.sim)
#plot(sm.net.log, param=sm.par.sim,file="sm.lin.log") #save to pdf



#The "abc" object, contains the following components:
# adj.values : The regression adjusted values, when method is "loclinear" or "neuralnet".
# unadj.values: The unadjusted values that correspond to "rejection" method.

#for more details see help (?abc)



###################################
## Posterior Probabilities Plots ##
###################################

#### Panmixia (loc linear) ####

#pdf("panmixia_posteriores.pdf", width=12,height=7)
par(mfrow=c(1,2), oma = c(0,0,2,0))

##Prior, and posterior probabilities for Neolithic size (N) parameter
den.lin.pm.logit_N<-density(pm.lin.logit$adj.values[,1]) #with logit transformation
den.lin.pm.log_N<-density(pm.lin.log$adj.values[,1]) #with log transformation

plot(den.lin.pm.logit_N, xlab=expression(paste(italic(N[N])," posteriors")),main="") 
lines(den.lin.pm.log_N,col="grey") 
lines(density(pm.lin.logit$unadj.values[,1]), lty=2) #rejection only
lines(density(pm.par.sim[,1]),lty=3) #prior
legend("topright",c("prior","rejection only","loclinear with logit","loclinear with log"), lty=c(3,2,1,1),col=c(1,1,1,"grey"),bty="n")


##Prior, and posterior probabilities for Palaeolithic size (UP) parameter
den.lin.pm.logit_UP<-density(pm.lin.logit$adj.values[,2]) #with logit transformation
den.lin.pm.log_UP<-density(pm.lin.log$adj.values[,2]) #with log transformation

plot(den.lin.pm.log_UP, xlab=expression(paste(italic(N[UP])," posteriors")),main="",col="grey") 
lines(den.lin.pm.logit_UP)  
lines(density(pm.lin.log$unadj.values[,2]), lty=2) #rejection only
lines(density(pm.par.sim[,2]),lty=3) #prior
legend("topright",c("prior","rejection only","loclinear with logit","loclinear with log"), lty=c(3,2,1,1), col=c("black","black","black","grey"), bty="n")

par(mfrow=c(1,1))
title(main="Panmixia")
#dev.off()


#### Split (loclinear) ####

#pdf("split_posteriores.pdf", width=12,height=7)
par(mfrow=c(1,2), oma = c(0,0,2,0))

##Prior, and posterior probabilities for Neolithic size (N) parameter
den.lin.sm.logit_N<-density(sm.lin.logit$adj.values[,1]) #with logit transformation
den.lin.sm.log_N<-density(sm.lin.log$adj.values[,1]) #with log transformation

plot(den.lin.sm.logit_N, xlab=expression(paste(italic(N[N])," posteriors")),main="") 
lines(den.lin.sm.log_N,col="grey")
lines(density(sm.lin.logit$unadj.values[,1]), lty=2) #rejection only
lines(density(sm.par.sim[,1]),lty=3) #prior
legend("topright",c("prior","rejection only","loclinear with logit","loclinear with log"), lty=c(3,2,1,1), col=c("black","black","black","grey"), bty="n")


##Prior, and posterior probabilities for Palaeolithic size (UP) parameter
den.lin.sm.logit_UP<-density(sm.lin.logit$adj.values[,2]) #with logit transformation
den.lin.sm.log_UP<-density(sm.lin.log$adj.values[,2]) #with log transformation

plot(den.lin.sm.logit_UP, xlab=expression(paste(italic(N[UP])," posteriors")),main="") 
lines(den.lin.sm.log_UP,col="grey")  
lines(density(sm.lin.log$unadj.values[,2]), lty=2) #rejection only
lines(density(sm.par.sim[,2]),lty=3) #prior
legend("topright",c("prior","rejection only","loclinear with logit","loclinear with log"), lty=c(3,2,1,1), col=c("black","black","black","grey"), bty="n")

par(mfrow=c(1,1))
title(main="Split")

#dev.off()



