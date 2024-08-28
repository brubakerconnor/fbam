### functions to extract information form DP MCMC chains

#### auxiliary function to detect na values derived from possible errors ####
sum_nna<-function (a){sum( !is.na(a) ) }

#####  auxiliary function for numerical computation of integrals: Simpson method ####
simpson <- function(fun,h) {
  # numerical integral using Simpson's rule
  # assume a < b and n is an even positive integer
  n=length(fun)
  if (n == 2) {
    s <- fun[1] + 4*fun[2] +fun[3]
  } else {
    s <- fun[1] + fun[n] + 2*sum(fun[seq(3,n-2,by=2)]) + 4 *sum(fun[seq(2,n-1, by=2)] )
  }
  s <- s*h/3
  return(s)
}

#### Integrated absolute error function ####
IAE<-function(x,y,h){
  return( simpson(  abs(x-y),h ) )  
}

#### auxiliary function for numerical computation of integrals: trapezoidal integral ####
trapezoid_integral=function(v, h ){h*(mean(v[1],v[length(v)])+sum(v[2:(length(v)-1)]) )}

#### Bartlett–Priestley kernel ####
BPresker<-function(x, h  ){
  y=x  
  for(i in 1:length(x)){
    if(abs(x[i]) <h ){
      y[i]<-  (1- (x[i]/h)^2 )
    }else{
      y[i]<-   0
    }
  }
  return(y)  
}

##### Bartlett–Priestley kernel smoothing ####
kernBP<-function(x, Y, X, h){ 
  DD<-c()
  for(i in 1:length(x)){
    DD[i]<-    sum(   BPresker(x[i]-X   ,h)*Y  )          /sum(  BPresker(x[i]-X   ,h) ) 
  }
  return(DD)
}

#### Bartlett–Priestley kernel smoothing MSE ####
JBP=function(h, dataX, dataY){
  fhati<-rep(0,length(dataX))
  for(i in 1:length(dataX)){
    fhati[i]<- kernBP(x=dataX[i], X=dataX[-i], Y=dataY[-i], h = h)
  }
  return(mean( (fhati - dataY)^2 ))
}

#### function to extract the mode of components parameters ####
#this library is for clustering the parameter based on a Gaussian mixture
#library(mclust) # loaded in the wrapper code BMARD_Replication.R

# The next function extract the MCMC iterations that achieve high loglikelihood values,
# quant define the quantile level for the loglikelihood values .95 is used in our paper
# persamp define the percetage of the last MCMC iterations from the chains used in the computations before filtering by likelihood 
# qlow and qsup are the limits for the pointwise curves estimated for each chain and globally
# threshold_weight=.001 filter the components with weigth below the threshold
# data is assumed to be a list result from the function SpectralBDP()
# chains is the numbers of chains specified to run the method in parallel with mcapply()
# Nsamp is the number of MCMC iterations
# freq_rate is the number of points taken from the epoch to apply the method

modes_extraction_gaussianmix<-function(data, chains,Nsamp,threshold_weight=.001 ,  freq_rate, quant,percsamp,qlow, qsup  ){
  #variables definiton
  sampconsider<-round( Nsamp*percsamp,0)
  psimatrix<-matrix( rep(NA,10*chains), nrow = chains )
  weightmatrix<-matrix( rep(NA,10*chains), nrow = chains )
  BWmatrix<-matrix( rep(NA,10*chains), nrow = chains )
  meanmatrix<-matrix( rep(NA,(freq_rate/2-1)*chains), nrow = chains )
  medianmatrix<-matrix( rep(NA,(freq_rate/2-1)*chains), nrow = chains )
  qsupmatrix<-matrix( rep(NA,(freq_rate/2-1)*chains), nrow = chains )
  qlowmatrix<-matrix( rep(NA,(freq_rate/2-1)*chains), nrow = chains )
  alllambda<-c()
  extweight<-c()
  psivecs<-c()
  BEWvecs<-c()
  
  if(chains>1){
    
    for(h in 1:chains){
      alllambdatemp<-c()
      extweighttemp<-c()
      psivecstemp<-c()
      BEWvecstemp<-c()
      # extraction of iteration with high log-likelihood values
      MLElist<- which( data[[h]]$Whittle[(   sampconsider  +1 )  :Nsamp] > quantile(data[[h]]$Whittle[(   sampconsider  +1 )  :Nsamp],quant ))
      indexes_chain<- (sampconsider+ MLElist)     
      for(k in 1:length(MLElist) ){
        extweighttemp<-c( extweighttemp, data[[h]]$weights[indexes_chain[k],1:data[[h]]$Lambda[( indexes_chain[k])] ]  )
        psivecstemp<-c(psivecstemp, data[[h]]$psi[indexes_chain[k],1:data[[h]]$Lambda[( indexes_chain[k])]] )
        BEWvecstemp<-c(BEWvecstemp, data[[h]]$BW[indexes_chain[k],1:data[[h]]$Lambda[( indexes_chain[k])]])
      }
      alllambdatemp<-data[[h]]$Lambda[( indexes_chain)]
      alllambda<-c(alllambda, alllambdatemp)  
      extweight<-c(extweight,extweighttemp)
      psivecs<-c(psivecs,psivecstemp)
      BEWvecs<-c(BEWvecs,BEWvecstemp)
      curvesmatrixtemp<-data[[h]]$CurveEsts[indexes_chain,]
      #putting ell the data together
      if(h==1){ curvesmatrix<-curvesmatrixtemp}else{
        curvesmatrix<-rbind(curvesmatrix, curvesmatrixtemp )
      }
      #possible candidates for the mean values
      priormeans =as.numeric(names (  table (alllambdatemp) ))
      df<-data.frame(psi= psivecstemp,BW= BEWvecstemp, weight= extweighttemp )
      iclvec<-c()
      listcenters=list()
      #compute cluster based on the possible cluster values taken form the mcmc and decide base on
      #maximizing the integrated conditional likelihood
      
      for( g in  1: length( priormeans)  ){
        
        k2 <-suppressWarnings( Mclust(df,G=priormeans[g]) )
        iclvec[g]<- k2$icl #to maximize
        listcenters[[g]]<-k2$parameters$mean 
      }
      
      indexlambda <-which.max(iclvec)
      
      # computing pointwise mean, median, and quanitle curves
      mediancurve=apply( curvesmatrixtemp,2,median  )
      meancurve=apply( curvesmatrixtemp,2,mean  )
      quantilcurve1<-apply( curvesmatrixtemp,2,quantile,qsup  )
      quantilcurve2<-apply( curvesmatrixtemp,2,quantile,qlow  )
      
      # matrix filling of individual chain results
      singcomp<-which( listcenters[[indexlambda]][3,  ]>threshold_weight )
      psimatrix[h, 1:length(singcomp) ]<-as.vector( listcenters[[indexlambda]][1,singcomp]   )
      BWmatrix[h, 1:length(singcomp)]<-as.vector( listcenters[[indexlambda]][2,singcomp]   )
      weightmatrix[h, 1:length(singcomp)]<-as.vector( listcenters[[indexlambda]][3,singcomp]  ) 
      meanmatrix[h, ]<-meancurve
      medianmatrix[h,]<-mediancurve
      qsupmatrix[h,]<-quantilcurve1
      qlowmatrix[h,]<-quantilcurve2
    }#end for h
    
    #same the method of clustering and decision but with all the chains together,   
    # also are computer as global components , and global curves   
    priormeans =as.numeric(names (  table (alllambda) ))
    df<-data.frame(psi= psivecs,BW= BEWvecs, weight= extweight )
    iclvec<-c()
    listcenters=list()
    for( g in  1: length( priormeans)  ){
      k2 <-suppressWarnings( Mclust(df,G=priormeans[g]) )
      iclvec[g]<- k2$icl 
      listcenters[[g]]<-k2$parameters$mean 
    }
    if(length( priormeans)==1){indexlambda=1}else{
      indexlambda <-which.max(iclvec)
    }
    mediancurveglobal=apply( curvesmatrix,2,median  )
    meancurveglobal=apply( curvesmatrix,2,mean  )
    quantilcurve1global<-apply( curvesmatrix,2,quantile,qsup  )
    quantilcurve2global<-apply( curvesmatrix,2,quantile,qlow  )
    return(list(
      psimatrix=psimatrix,
      BWmatrix=BWmatrix,
      weightmatrix=weightmatrix,
      meanmatrix=meanmatrix,
      medianmatrix=medianmatrix,
      qsupmatrix=qsupmatrix,
      qlowmatrix=qlowmatrix, 
      globalcenter=listcenters[[indexlambda]],
      mediancurveglobal=mediancurveglobal,
      meancurveglobal=meancurveglobal,
      qsupcurveglobal=quantilcurve1global,
      qlowcurveglobal=quantilcurve2global
    ) )
  } # end if chains 
  else{
    
    alllambdatemp<-c()
    extweighttemp<-c()
    psivecstemp<-c()
    BEWvecstemp<-c()
    
    MLElist<- which( data$Whittle[(   sampconsider  +1 )  :Nsamp] > quantile(data$Whittle[(   sampconsider  +1 )  :Nsamp],quant ))
    indexes_chain<- (sampconsider+ MLElist)     
    for(k in 1:length(MLElist) ){
      extweighttemp<-c( extweighttemp, data$weights[indexes_chain[k],1:data$Lambda[( indexes_chain[k])] ]  )
      psivecstemp<-c(psivecstemp, data$psi[indexes_chain[k],1:data$Lambda[( indexes_chain[k])]] )
      BEWvecstemp<-c(BEWvecstemp, data$BW[indexes_chain[k],1:data$Lambda[( indexes_chain[k])]])
    }
    alllambdatemp<-data$Lambda[( indexes_chain)]
    alllambda<- alllambdatemp
    extweight<-extweighttemp
    psivecs<-psivecstemp
    BEWvecs<-BEWvecstemp
    curvesmatrixtemp<-data$CurveEsts[indexes_chain,]
    
    priormeans =as.numeric(names (  table (alllambdatemp) ))
    df<-data.frame(psi= psivecstemp,BW= BEWvecstemp, weight= extweighttemp )
    iclvec<-c()
    listcenters=list()
    for( g in  1: length( priormeans)  ){
      k2 <-suppressWarnings( Mclust(df,G=priormeans[g]))
      iclvec[g]<- k2$icl 
      listcenters[[g]]<-k2$parameters$mean 
    }
    if(length( priormeans)==1){indexlambda=1}else{
      indexlambda<- which.max(iclvec)
    }
    mediancurve=apply( curvesmatrixtemp,2,median  )
    meancurve=apply( curvesmatrixtemp,2,mean  )
    quantilcurve1<-apply( curvesmatrixtemp,2,quantile,qsup  )
    quantilcurve2<-apply( curvesmatrixtemp,2,quantile,qlow  )
    
    # matrix filling
    singcomp<-which( listcenters[[indexlambda]][3,  ]>threshold_weight )
    psimatrix[1, 1:length(singcomp) ]<-as.vector( listcenters[[indexlambda]][1,singcomp]   )
    BWmatrix[1, 1:length(singcomp)]<-as.vector( listcenters[[indexlambda]][2,singcomp]   )
    weightmatrix[1, 1:length(singcomp)]<-as.vector( listcenters[[indexlambda]][3,singcomp]  ) 
    meanmatrix[1, ]<-meancurve
    medianmatrix[1,]<-mediancurve
    qsupmatrix[1,]<-quantilcurve1
    qlowmatrix[1,]<-quantilcurve2 
    return(list(
      psimatrix=psimatrix,
      BWmatrix=BWmatrix,
      weightmatrix=weightmatrix,
      meanmatrix=meanmatrix,
      medianmatrix=medianmatrix,
      qsupmatrix=qsupmatrix,
      qlowmatrix=qlowmatrix
    ) )
  } 
  
  return(list(
    psimatrix=psimatrix,
    BWmatrix=BWmatrix,
    weightmatrix=weightmatrix,
    meanmatrix=meanmatrix,
    medianmatrix=medianmatrix,
    qsupmatrix=qsupmatrix,
    qlowmatrix=qlowmatrix, 
    globalcenter=listcenters[[indexlambda]],
    mediancurveglobal=mediancurveglobal,
    meancurveglobal=meancurveglobal,
    qsupcurveglobal=quantilcurve1global,
    qlowcurveglobal=quantilcurve2global
  ) )
} #end function


#auxiliar library to simulate AR(2) process
#library(TSA)

#### function to simualte AR2 mixture from components ####
# psi is a vector of phase parameter of the characteristic roots of the AR(2) 
# phi is a vector of weights of the mixture
# L is a vector of the log of the modulus of the characteristic roots of the AR(2) process
# n is the numer of points to simulate
# sd is the standard deviation of the time series
AR2mixsimulation<-function(psi, phi, L, n, sd ){
  #variables definition  
  phi1<- exp(-L)*2*cos(2*pi*psi)
  phi2<--exp(-2*L) 
  simindividual<-matrix(rep(1:n,length(psi)  ), nrow = length(psi) )
  
  for( i in 1:length(psi) ){
    simindividual[i,]<-sqrt(phi[i])*arima.sim(model=list(ar=c(phi1[i], phi2[i] ) )  , n=n)    
  }
  simmixture<-apply(simindividual, 2, sum)  
  
  return( list( individualmatrix=simindividual, simmixture=simmixture   )  )  
}


  
  
