#BMARD EXAMPLE

#### set the following parameters ####

#number of mcmc iterations
Nsamp=1000
#number of chains
chains<-3
#time series length
size=1000
#sample frequency 
Sample_frequency=1000  # measured in Hertz (Hz), time points per second.


### loodaing and installing requiered packages

## next code is taken from the next article 
## https://vbaliga.github.io/verify-that-r-packages-are-installed-and-loaded/
## If a package is installed, it will be loaded. If any 
## are not, the missing package(s) will be installed 
## from CRAN and then loaded.

## First specify the packages of interest
packages = c("Rcpp", "TSA",
             "doParallel", "mclust")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)


## sourcing the necessary code
 codepath<-dirname(rstudioapi::getActiveDocumentContext()$path )
Rcpp::sourceCpp(paste(codepath, '/CPPcode/BMARD_V112020.cpp',sep="" ) ) 
source(paste(codepath, '/Rcode/ExtractionBMARDmaincomponentsmodes.R',sep="" ))


# time series lenght is 1000 points in this test 500 in the paper
# last value of the periodogram filtered, (frequency = .5)
# epsilon and epsilonBW are specified based on practical reasons are the MH step proposal
# for other parameters see specification in cpp code 

# function to run chains in parallel
#  Lambda_init = sample(1:20,1) set at random an intial number of components
multirun <- function(S){  
  trunc<- sample(20:30,1); # random initial level of trunctation
  SpectralBDP(Tsize = 500, omegas =omegas, Pgram = S , Nsamp=Nsamp, L=trunc, epsilon=rep( runif(1, .01,.05 ) ,trunc)  , epsilon_BW=runif(1,.01,.05)  , SamplingRate = 1, 
              MaxFreq = 0.5, Sup = Nsamp, Lambda_init = sample(1:20,1) , Lambda_prior = "flat", 
              Lambda_max = 50L, Lambda_prop_dist = "up1down1", alpha_const = F, 
              alpha_init = 1, a_alpha_L = .1, b_alpha_L = .1, q = 1) 
}


## same definiton as in Section 3
#weights
phi<-c(.1,.6,.3)
#phase parameters at 8, 30, and 60 Hz for the particular time series length
psi<- (1/Sample_frequency)*c(4,30,60)
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
#frequency values
omegas =seq(1/size , .5-1/size,by=1/size)

### true AR2 mixture SDF  
W1<-ARMAspec(model=list(ar = par1 ), plot = F , freq = omegas )
W2<-ARMAspec(model=list(ar = par2 ), plot = F , freq = omegas )
W3<-ARMAspec(model=list(ar = par3 ), plot = F , freq = omegas )
specmix<-phi[1]*W1$spec+phi[2]*W2$spec+phi[3]*W3$spec
specmix<-specmix/simpson(specmix,1/size ) 

## BMARD fitting
B<-list()
for(k in 1:chains){   
B[[k]]<- periodogram(S, plot=F)$spec[-250]            
}
BDP<-mclapply(B, multirun, mc.cores = 3)


#extraction of main components

listmodes<- modes_extraction_gaussianmix(data=BDP,threshold_weight=.005, chains=3,Nsamp=Nsamp, freq_rate=size, quant=.9,percsamp=.6 ,qlow=.025, qsup=.975  )

# main global curves
plot(omegas, B[[1]] ,type = "l",col=4, xlab = "Frequencies", ylab = "Power",main = "BMARD applied to a mixture of AR2(2) processes" ,   xlim = c(0,.1))
lines(omegas, listmodes$mediancurveglobal,type = "l" )
lines(omegas, listmodes$qsupcurveglobal,type = "l", col=2 )
lines(omegas, listmodes$qlowcurveglobal,type = "l",col=3 )
lines(omegas, specmix ,type = "l",col=5, lwd=3 )

# components 
abline(v=listmodes$globalcenter[1,], col=6)
legend(.06, 150, legend = c("BMARD median","BMARD 5% quantile", "BMARD 95% quantile", "True Spectrum", "MAP locations", "perdiodogram"), col=c(1,3,2,5,6,4) , bty="n", cex=.5 , lwd=c(1,1,1,3,1,1) )


##### application to a AR(12) process ####

S<-arima.sim(n = size, list(ar = c(0,0,0,0.9,0,0,0,.7,0,0,0,-.63)  ), sd = 1, n.start = 10000)
#standarized serie
S<-(S-mean(S))/sd(S)

## AR(12) model 

W<-ARMAspec(model=list(ar = c(0,0,0,0.9,0,0,0,.7,0,0,0,-.63) ), plot = F,freq = omegas)
Wint<- simpson(W$spec,1/size )
Wstand<-W$spec/Wint

B<-list()
for(k in 1:chains){   
B[[k]]<- periodogram(S, plot=F)$spec[-250]            
}

#run the method in parallel, the cores is thinking the cpus are going to be partitioned in 3 clusters using 33% of the cpus each with gives the best performance
BDP<-mclapply(B, multirun, mc.cores = 3)


#extraction of main components

listmodes<- modes_extraction_gaussianmix(data=BDP,threshold_weight=.005, chains=3,Nsamp=Nsamp, freq_rate=size, quant=.9,percsamp=.6 ,qlow=.025, qsup=.975  )

# main global curves
plot(omegas, B[[1]] ,type = "l",col=4, xlab = "Frequencies", ylab = "Power",main = "BMARD applied to a AR(12) process"    )
lines(omegas, listmodes$mediancurveglobal,type = "l" )
lines(omegas, listmodes$qsupcurveglobal,type = "l", col=2 )
lines(omegas, listmodes$qlowcurveglobal,type = "l",col=3 )
lines(omegas, Wstand ,type = "l",col=5, lwd=2 )

# components 
abline(v=listmodes$globalcenter[1,],col=6)
legend(.06, 100, legend = c("BMARD median","BMARD 5% quantile", "BMARD 95% quantile", "True Spectrum", "MAP locations", "perdiodogram"), col=c(1,3,2,5,6,4) , bty="n", cex=.5 , lwd=c(1,1,1,3,1,1) )




#### application to a MA(4) process ####

S<-arima.sim(n = size, list(ma=c(-.3,-.6,-.3,.6)  ), sd = 1, n.start = 10000)
#standarized serie
S<-(S-mean(S))/sd(S)

## MA(4) process
W<-ARMAspec(model=list(ma = c(-.3,-.6,-.3,.6 ) ), plot = F,freq = omegas)
Wint<- simpson(W$spec,1/size )
WMA4st<-W$spec/Wint

B<-list()
for(k in 1:chains){   
  B[[k]]<- periodogram(S, plot=F)$spec[-250]            
}

#run the method in parallel, the cores is thinking the cpus are going to be partitioned in 3 clusters using 33% of the cpus each with gives the best performance
BDP<-mclapply(B, multirun, mc.cores = 3)


#extraction of main components

listmodes<- modes_extraction_gaussianmix(data=BDP,threshold_weight=.005, chains=3,Nsamp=Nsamp, freq_rate=size, quant=.9,percsamp=.6 ,qlow=.025, qsup=.975  )

# main global curves
plot(omegas, B[[1]] ,type = "l",col=4, xlab = "Frequencies", ylab = "Power",main = "BMARD applied to a AR(12) process"    )
lines(omegas, listmodes$mediancurveglobal,type = "l" )
lines(omegas, listmodes$qsupcurveglobal,type = "l", col=2 )
lines(omegas, listmodes$qlowcurveglobal,type = "l",col=3 )
lines(omegas, WMA4st ,type = "l",col=5, lwd=2 )

# components 
abline(v=listmodes$globalcenter[1,], col=6)
legend(.01, 10, legend = c("BMARD median","BMARD 5% quantile", "BMARD 95% quantile", "True Spectrum", "MAP locations", "perdiodogram"), col=c(1,3,2,5,6,4) , bty="n", cex=.5 , lwd=c(1,1,1,3,1,1) )

