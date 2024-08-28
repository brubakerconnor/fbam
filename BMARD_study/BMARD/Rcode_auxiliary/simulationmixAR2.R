# Code to apply the BMARD method to a simulation defined from a mixture of AR(2) processes
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
  SpectralBDP(Tsize = 500, omegas =seq(.002,.498,by=.002), Pgram = S , Nsamp=Nsamp, L=trunc, epsilon=rep( runif(1, .01,.05 ) ,trunc)  , epsilon_BW=runif(1,.01,.05)  , SamplingRate = 1, 
              MaxFreq = 0.5, Sup = Nsamp, Lambda_init = sample(1:20,1) , Lambda_prior = "flat", 
              Lambda_max = 50L, Lambda_prop_dist = "up1down1", alpha_const = F, 
              alpha_init = 1, a_alpha_L = .1, b_alpha_L = .1, q = 1) 
}

#number of mcmc iterations
Nsamp=100000
#number of chains
chains<-6

#time series length
size=500

#weights
phi<-c(.1,.6,.3)
#phase parameters at 8, 30, and 60 Hz for the particular time series length
psi<- (1/500)*c(4,15,30)
#log modulus parameter scale of the kernel
L<-rep(.03, 3)  

#parameters in the time domain
par1<- c(exp(-L[1])*2*cos(2*pi*psi[1]), -exp(-2*L[1])   )
par2<-c(exp(-L[2])*2*cos(2*pi*psi[2]), -exp(-2*L[2])   )
par3<-c(exp(-L[3])*2*cos(2*pi*psi[3]), -exp(-2*L[3])   )
#individual time series simulations
S1<-arima.sim(n = size, list(ar = par1 ), sd = 1, n.start = 10000)
S2<-arima.sim(n = size, list(ar = par2 ), sd = 1, n.start = 10000)
S3<-arima.sim(n = size, list(ar = par3 ), sd = 1, n.start = 10000)
#mixture of AR(2) processes
S<-(phi[1]^.5)*S1+(phi[2]^.5)*S2+(phi[3]^.5)*S3

#standarized series
S<-(S-mean(S))/sd(S)

# for use the simulated data load the matrix by chaing accordingly the directory where you saved the matrix
database<-readRDS("~/pgrammatrixAR2.rds")

B<-list()
for(k in 1:chains){   
  #uncomment the appropiate option for new simualtion (first) or the computed periodogram (second)
  
  #   B[[k]]<- periodogram(S, plot=F)$spec[-250]            
  #  B[[k]]<-database[as.numeric(argsbash)  ,] 
}

BDP<-mclapply(B, multirun, mc.cores = 3)

#run the method in parallel, the cores is thinking the cpus are going to be partitioned in 3 clusters using 33% of the cpus each with gives the best performance
saveRDS(BDP, file= paste( "~/BayesiansimmixAR2_V"  ,argsbash,  ".rds", sep = "" )   ) 
