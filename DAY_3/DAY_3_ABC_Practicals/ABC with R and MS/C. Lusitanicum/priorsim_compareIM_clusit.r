#Example use of ms to infer parameters with ABC in a model of two populations
#that have diverged without migration 

# define the number of simulations
nsim <- 100000
skipms <- F

# sample theta1, theta2, theta3 (for comparison with IM)
theta1 <- runif(nsim,max=10)
theta2 <- runif(nsim,max=10)
theta3 <- runif(nsim,max=10)

# sample the time of divergence (in generations)
t1 <- runif(nsim,max=10)

theta <- theta1
T1 <- t1/theta1
N2 <- theta2/theta1
N3 <- theta3/theta1

# make a matrix of parameter vals
params <- cbind(theta,N2,T1,T1,N3)

# write this to a file for ms to read in using tbs
write(t(params),file="params.txt",ncol=5)

# perform the ms command
if(!skipms)shell(paste("ms 41",formatC(nsim,format="d")," -t tbs -I 2 15 26 -n 2 tbs  -ej tbs 2 1 -en tbs 1 tbs < params.txt > file1"))

# compute summary stats
shell("sample_stat_mab < file1 > res2")
shell("sample_stat_mab 0 14 < file1 > res3")
shell("sample_stat_mab 15 40 < file1 > res4")

#read summary stats into R

res2 <- read.table("res2")
res3 <- read.table("res3")
res4 <- read.table("res4")

# put the summ stats into a matrix
sstable <- data.matrix(cbind(res2[c(2,4,6,8)],res3[c(2,4,6,8)],res4[c(2,4,6,8)]))

#get the ABC function

source("make_pd2.r")

#define the target summary statistics
target <- scan("clusit_ss.txt")

# using the regression adjustment method obtain approximate samples
# from the posterior distribution

res.theta1 <- makepd4(target,theta1,sstable,0.01,rej=F,transf="logit",bb=c(0,10))
res.theta2 <- makepd4(target,theta2,sstable,0.01,rej=F,transf="logit",bb=c(0,10))
res.theta3 <- makepd4(target,theta3,sstable,0.01,rej=F,transf="logit",bb=c(0,10))
res.t1 <- makepd4(target,t1,sstable,0.01,rej=F,transf="logit",bb=c(0,10))

par(ask=T)
# plot prior for theta in population 1
plot(density(theta1),lty=2,lwd=3,,xlim=c(0,10),ylim=c(0,0.5),xlab="theta1",main="Fish ABC analysis with ms")
#plot posterior for theta in population 1
lines(density(res.theta1$x),col="blue",lwd=2)
dev.print(file="fish_theta1_abc.pdf",dev=pdf)


# plot prior for theta in population 1
plot(density(theta2),lty=2,lwd=3,,xlim=c(0,10),ylim=c(0,1),xlab="theta2",main="Fish ABC analysis with ms")
#plot posterior for theta in population 1
lines(density(res.theta2$x),col="blue",lwd=2)
dev.print(file="fish_theta2_abc.pdf",dev=pdf)


# plot prior for theta in population 1
plot(density(theta3),lty=2,lwd=3,,xlim=c(0,10),ylim=c(0,0.5),xlab="theta3",main="Fish ABC analysis with ms")
#plot posterior for theta in population 1
lines(density(res.theta3$x),col="blue",lwd=2)
dev.print(file="fish_theta3_abc.pdf",dev=pdf)


# plot prior for theta in population 1
plot(density(t1),lty=2,lwd=3,,xlim=c(0,10),ylim=c(0,2),xlab="t1",main="Fish ABC analysis with ms")
#plot posterior for theta in population 1
lines(density(res.t1$x),col="blue",lwd=2)
dev.print(file="fish_t1_abc.pdf",dev=pdf)

par(ask=F)




