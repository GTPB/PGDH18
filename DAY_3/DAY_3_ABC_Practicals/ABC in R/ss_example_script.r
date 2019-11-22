# This exercise approximates the posterior distribution for theta using 
# Approximate Bayesian computation (ABC). 

# The data are microsatellites.

# The summary statistics we use are 

# * heterozygosity
# * variance of allele length
# * number of alleles.

# In order to determine how good the ABC algorithm is 
# we use the result of a full likelihood method
# 
# This "real" posterior distribution, was obtained from an MCMC run using MSVAR
# This output is in the file called "out1"


m1 <- t(matrix(scan("out1"),nrow=10))
plot(density(m1[,4]))



# 2. Source the functions necessary to do the ABC


source("urn_example_r.txt")
source("make_pd2.r")

# 3. Input the microsatellites alleles lengths


target <- scan("targetvals.txt")


# 3. Get values of the summary statistics:
# het
# var
# no. alleles  
# from the observed data

ss.target <- cbind(het.measure(target),var(target),length(table(target)))

# 4. Simulate values of theta from a prior.

theta.sim <- 10^(runif(2000,-2,2))

# 5. Run the microsatellite genealogical simulations

ss.sim <- simulate.hets(100,theta.sim,50)

# 6. Estimate the posterior distribution using 
# ABC with SIMPLE REJECTION

theta.ss.pdf <- reject.func(log10(theta.sim),ss.sim,ss.target,0.5)

# 7. estimate the posterior distribution using 
# ABC with REGRESSION

result <- makepd4(ss.target,log10(theta.sim),ss.sim,0.1,rejmethod=F)


# 8. Plot it and compare with the "true" posterior distribution.


lines(density(theta.ss.pdf),col="red")
lines(density(result$x),col="green")



#######################################################################
#
# EXERCICE:
# a. Compare the effect of different tolerances in rejection
# b. Compare pure rejection with regression.
# c. Simulate other data sets with different (or same) value of theta.
#
#######################################################################