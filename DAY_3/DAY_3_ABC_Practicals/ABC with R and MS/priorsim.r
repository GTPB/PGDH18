#Example use of ms to infer parameters with ABC in a model of two populations
#that have diverged without migration 

# define the number of simulations
nsim <- 100000

# sample the population sizes of pops 1, 2, and 3 (3 is ancestral)
n1 <- 10^runif(nsim,max=4)
n2 <- 10^runif(nsim,max=4)
n3 <- 10^runif(nsim,max=4)

# sample the time of divergence (in generations)
t1 <- 10^runif(nsim,max=5)

# length of sequence

length=5000

# set the reference population (for scaling everything in ms)

n0 <- n1

# sample mutation rates (per site)
mu <- 10^runif(nsim,min=-9,max=-7)

# compute scaled parameters: theta, scaled time, and two size ratios (relative to pop 1)

theta <- 4*n0*mu*length

T1 <- t1/(4*n0)
N2 <- n2/n0 
N3 <- n3/n0

# make a matrix of parameter vals
params <- cbind(theta,N2,T1,T1,N3)
lparams <- log10(params)

# write this to a file for ms to read in using tbs
write(t(params),file="params.txt",ncol=5)

# perform the ms command
shell(paste("ms 40",formatC(nsim,format="d")," -t tbs -I 2 20 20 -n 2 tbs  -ej tbs 2 1 -en tbs 1 tbs < params.txt > file1"))

# compute summary stats
shell("sample_stat_mab < file1 > res2")
shell("sample_stat_mab 0 19 < file1 > res3")
shell("sample_stat_mab 20 39 < file1 > res4")

#read summary stats into R

res2 <- read.table("res2")
res3 <- read.table("res3")
res4 <- read.table("res4")

# put the summ stats into a matrix
sstable <- data.matrix(cbind(res2[c(2,4,6,8)],res3[c(2,4,6,8)],res4[c(2,4,6,8)]))

#get the ABC function

source("make_pd2.r")

#define the target summary statistics
target <- scan("target_summstats.txt")

# using the regression adjustment method obtain approximate samples
# from the posterior distribution

res.theta <- makepd4(target,lparams[,1],sstable,0.01,rej=F)
res.n2 <- makepd4(target,lparams[,2],sstable,0.01,rej=F)
res.t1 <- makepd4(target,lparams[,3],sstable,0.01,rej=F)
res.n3 <- makepd4(target,lparams[,5],sstable,0.01,rej=F)

# plot prior for theta in population 1
plot(density(lparams[,1]),lty=2,lwd=3,ylim=c(0,0.5),xlim=c(-5,3))
#plot posterior for theta in population 1
lines(density(res.theta$x),col="blue",lwd=2)



# Exercises
# 1) Explore the effect of different priors. For example try to give a flat uniform on 
#     parameters used by IM, so that we can compare with IM. 
# 2) Explore different demographic models (you may want to simulate test data from these).
#
