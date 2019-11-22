#Example use of ms to infer parameters with ABC in a model of two populations
#that have diverged without migration 

#NOTE!!! The example given here uses the parameterization that is commonly used in 
#        IM. In IM the times are divided by mutation rate, whereas in ms the times 
#        are divide by 4N0 where N0 is the 'reference' population. 

# define the number of simulations
nsim <- 100000

# sample theta1, theta2, theta3 (for comparison with IM)
# Note that IM assumes diploids so theta = 4Nmu
# runif is a uniform random number generator and we are asking for nsim random numbers between 0 and 10.
theta1 <- runif(nsim,max=10)
theta2 <- runif(nsim,max=10)
theta3 <- runif(nsim,max=10)

# sample the time of divergence (in generations divided by mutation rate - for comparison with IM)
t1 <- runif(nsim,max=10)

#ms requires that all the parameters be scaled in a certain way. We take the theta of pop 1 (theta1) to be 4N0mu
theta <- theta1
T1 <- t1/theta1      # i.e. (T*mu)/4N0mu = T/4N0
N2 <- theta2/theta1
N3 <- theta3/theta1

# make a matrix of parameter vals. cbind() creates matrices by gluing vectors column-wise [rbind() does it row-wise]
params <- cbind(theta,N2,T1,T1,N3)

# write this to a file for ms to read in using tbs
write(t(params),file="params.txt",ncol=5)

# perform the ms command
system(paste("./ms 40",formatC(nsim,format="d")," -t tbs -I 2 20 20 -n 2 tbs  -ej tbs 2 1 -en tbs 1 tbs < params.txt > file1"))

#explanation:
#'ms' runs the ms executable. The remaining space-separated items are arguments to the ms command. 
#The order is very important.
#40                          is the total number of gene-copies to simulate. In this example there are two populations and 
#                              we want 20 copies to be sampled from each population. So the total number of copies is 40.
#formatC(nsim,format="d")    is a complicated way of writing 'nsim' [we need stop R writing a decimal point].
#-t tbs                      take the first column from the file params.txt (theta for deme 1)
#-I 2 20 20                  there are 2 demes each with 20 gene copies sampled
#-n 2 tbs                    take the second column from the file params.txt and assign to theta for deme 2 
#-ej tbs 2 1                 at time given by the third column from the file params.txt lineages from deme 2 move to deme 1
#-en tbs 1 tbs               at time given by the fourth column from the file params.txt deme 1 changes in size to size given
#                              by the fifth column from the file params.txt.
#NOTE!!! Ignore warning messages that 'cmd' execution failed.

# compute summary stats
system("./sample_stat_mab < file1 > res2")       # summary stats for all demes together
system("./sample_stat_mab 0 19 < file1 > res3")  # summary stats for deme 1
system("./sample_stat_mab 20 39 < file1 > res4") # summary stats for deme 2
#NOTE!!! Ignore warning messages that 'cmd' execution failed.

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

res.theta1 <- makepd4(target,theta1,sstable,0.01,rej=F,transf="log")
res.theta2 <- makepd4(target,theta2,sstable,0.01,rej=F,transf="log")
res.theta3 <- makepd4(target,theta3,sstable,0.01,rej=F,transf="log")
res.t1 <- makepd4(target,t1,sstable,0.01,rej=F,transf="log")

#Note warning messages about collinearity. This reflects some feature of the 
#summary statistics of points closest to target data. 
#Worth investigating if you get such a message, but not necessarily a problem. 
#Ignore in this case.

library(locfit)
par(ask=T)
#plot posterior density for theta for population 1
plot(locfit(~res.theta1$x,xlim=c(0,10)),type="l",col="blue",xlab="Theta Deme 1")
#plot approximation of prior by running locfit on sims from prior (should be uniform on 0,10)
lines(locfit(~theta1,xlim=c(0,10)),lty=2)

#plot posterior density for theta for population 1
plot(locfit(~res.theta2$x,xlim=c(0,10)),type="l",col="blue",xlab="Theta Deme 2")
#plot approximation of prior by running locfit on sims from prior (should be uniform on 0,10)
lines(locfit(~theta2,xlim=c(0,10)),lty=2)

#plot posterior density for theta for ancestral population 
plot(locfit(~res.theta3$x,xlim=c(0,10)),type="l",col="blue",xlab="Theta Ancestral Population")
#plot approximation of prior by running locfit on sims from prior (should be uniform on 0,10)
lines(locfit(~theta3,xlim=c(0,10)),lty=2)

#plot posterior density for t/mu for ancestral population 
plot(locfit(~res.t1$x,xlim=c(0,10)),type="l",col="blue",xlab="T/mu")
#plot approximation of prior by running locfit on sims from prior (should be uniform on 0,10)
lines(locfit(~t1,xlim=c(0,10)),lty=2)
par(ask=F)

