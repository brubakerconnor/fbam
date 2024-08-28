# Code to apply the BMARD method to a simulation defined from a MA(4) process
# the periodogram matrix can be loaded or computed new simulations , uncomment the appropiate section

argsbash<- commandArgs()[6]

library(Rcpp)
library(TSA)
library(doParallel)

## change the directory accordingly
Rcpp::sourceCpp('~/BMARD_V112020.cpp')

# time series lenght is 500 points
# last value of the periodogram filtered, (frequency = .5)
# epsilon and epsilonBW are specified based on practical reasons are the MH step proposal
# for other parameters see specification in cpp code 
multirun <- function(S){  
  trunc<- sample(20:30,1);
  SpectralBDP(Tsize = 500, omegas =seq(.002,.498,by=.002), Pgram = S , Nsamp=Nsamp, L=trunc, epsilon=rep( runif(1, .05,.2)  ,trunc)  , epsilon_BW=runif(1,.05,.2)  , SamplingRate = 1, 
              MaxFreq = 0.5, Sup = Nsamp, Lambda_init = sample(1:20,1) , Lambda_prior = "flat", 
              Lambda_max = 50L, Lambda_prop_dist = "up1down1", alpha_const = F, 
              alpha_init = 1, a_alpha_L = .1, b_alpha_L = .1, q = 1) 
}

#number of mcmc iterations
Nsamp=100000
#number of chains
chains<-6

BDP=list()
size=500

# for new simulations use the following code and generate new MA(4) processes  
S<-arima.sim(n = size, list(ma=c(-.3,-.6,-.3,.6)  ), sd = 1, n.start = 10000)
#standarized serie
S<-(S-mean(S))/sd(S)

# for use the simulated data load the matrix by chaing accordingly the directory where you saved the matrix
database<-readRDS("~/pgrammatrixMA4.rds")

B<-list()
for(k in 1:chains){   
  #uncomment the appropiate option for new simualtion (first) or the computed periodogram (second)
  
  #   B[[k]]<- periodogram(S, plot=F)$spec[-250]            
  #  B[[k]]<-database[as.numeric(argsbash)  ,] 
}
#run the method in parallel, the cores is thinking the cpus are going to be partitioned in 3 clusters using 33% of the cpus each with gives the best performance
BDP<-mclapply(B, multirun, mc.cores = 3)


# save the results in a database for posterior use, change the directory accordingly 
saveRDS(BDP, file= paste( "~/BayesiansimulMA4_V"  ,argsbash,  ".rds", sep = "" )   ) 
