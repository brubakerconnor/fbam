### individual run of the Method of Choudhuri (2014) for a HPC cluster

# capture the variable from a for loop in bash
argsbash<- commandArgs()[6]

library(Rcpp)
#library(R.matlab)
#library(TSA)
library(doParallel)

# change the directories accordingly
Rcpp::sourceCpp( '~/CPPcode/NDPrcppBen_pgram.cpp' ) 



multirun <- function(S){  
  trunc<- sample(20:30,1);
  SpectralBDP(Pgram=S, omegas=omegas, Nsamp=Nsamp, L=trunc, epsilon=rep( runif(1, .01,.05)  ,trunc)    , SamplingRate = 1, 
              MaxFreq = 0.5, Sup = (Nsamp-1), Lambda_init = sample(1:20,1) , Lambda_prior = "flat", 
              Lambda_max = 300L, Lambda_prop_dist = "up1down1", alpha_const = TRUE, 
              alpha_init = 1, a_alpha_L = 1, b_alpha_L = 1) 
  }

omegas<-seq(.002, .498,by=.002)


# change the directory accordingly
# change the database name as necessary 
database<-readRDS("~/SimulationsDATABASES/pgrammatrixAR12.rds")

Nsamp=100000

chains<-6

BDP=list()

S<-database[    as.numeric( argsbash)     , ]   
  
B<-list()
for(k in 1:chains){   
  B[[k]]<- S
}

BDP<-mclapply(B, multirun, mc.cores = 3)

saveRDS(BDP, file= paste( "ChoudhurimixAR12_V"  , as.numeric(argsbash)  ,  ".rds", sep = "" )   ) 




