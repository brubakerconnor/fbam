# Code to extract the main information of the BMARD posterior samples



###### kernel and spline estimation  smoothing ####

## sourcing the necessary code
 codepath<-dirname(rstudioapi::getActiveDocumentContext()$path )
# source the code ExtractionBMARDmaincomponentsmodes.R
#adjust the directories accordingly
source(paste(codepath,"/Rcode_auxiliary/ExtractionBMARDmaincomponentsmodes.R",sep="" )   )
# optimize the Leave One Out Crossvalidation and find the optimal bandwidth for both method

#adjust the directories accordingly
MA4pgrammatrix<-readRDS(paste(codepath,"/SimulationsDATABASES/pgrammatrixMA4.rds",sep="" ))
AR12pgrammatrix<-readRDS(paste(codepath,"/SimulationsDATABASES/pgrammatrixAR12.rds",sep="" ))
pgram_matrix<-readRDS(paste(codepath,"/SimulationsDATABASES/pgrammatrixAR2.rds",sep="" ))

for(i in 1:1000){
  #spline ma4
  Spline<-smooth.spline(seq(.002,.498,by=.002), MA4pgrammatrix[i,], cv=T )
  splineMA4matrix[i,]<-Spline$y/simpson(Spline$y,.002)
  
  #spline AR2
  Spline<-smooth.spline(seq(.002,.498,by=.002), pgram_matrix[i,], cv=T )
  splineAR2matrix[i,]<-Spline$y/simpson(Spline$y,.002)
  
  #spline AR12
  Spline<-smooth.spline(seq(.002,.498,by=.002), AR12pgrammatrix[i,], cv=T )
  splineAR12matrix[i,]<-Spline$y/simpson(Spline$y,.002)
  
  #kernel MA4
  OP<-optimize(JBP,interval=c(.00001,.2) ,dataX=seq(.002,.498,by=.002), MA4pgrammatrix[i,] )
  KSM<-kernBP(x=seq(.002,.498,by=.002), Y=MA4pgrammatrix[i,], X=seq(.002,.498,by=.002), h=OP$minimum)
  kernelMA4matrix[i,]<-KSM/simpson(KSM,.002)
  
  # #kernel AR12
  OP<-optimize(JBP,interval=c(.00001,.2) ,dataX=seq(.002,.498,by=.002), AR12pgrammatrix[i,] )
  KSM<-kernBP(x=seq(.002,.498,by=.002), Y=AR12pgrammatrix[i,], X=seq(.002,.498,by=.002), h=OP$minimum)
  kernelAR12matrix[i,]<-KSM/simpson(KSM,.002)
  
  # #kernel AR2
  OP<-optimize(JBP,interval=c(.00001,.2) ,dataX=seq(.002,.498,by=.002), pgram_matrix[i,] )
  KSM<-kernBP(x=seq(.002,.498,by=.002), Y=pgram_matrix[i,], X=seq(.002,.498,by=.002), h=OP$minimum)
  kernelAR2matrix[i,]<-KSM/simpson(KSM,.002)
}

#save the databases for posterior analysis
saveRDS(splineMA4matrix,"splineMA4matrix.rds" )
saveRDS(splineAR2matrix,"splineAR2matrix.rds")
saveRDS(splineAR12matrix,"splineAR12matrix.rds")
saveRDS( kernelMA4matrix,"kernelMA4matrix.rds")
saveRDS(kernelAR12matrix,"kernelAR12matrix.rds")
saveRDS(kernelAR2matrix,"kernelAR2matrix")


#### true spectrums from the simulations setting ####

## MA(4) process
W<-ARMAspec(model=list(ma = c(-.3,-.6,-.3,.6 ) ), plot = F,freq = seq(0.002, 0.498, .002))
Wint<- simpson(W$spec,W$freq[1] )
WMA4st<-W$spec/Wint


###  AR2 mixture 

#weights
phi<-c(.1,.6,.3)
#phase parameters at 8, 30, and 60 Hz for the particular time series length
psi<- (1/1000)*c(4,30,60)
#log modulus parameter scale of the kernel
L<-rep(.03, 3)  
#parameters in the time domain
par1<- c(exp(-L[1])*2*cos(2*pi*psi[1]), -exp(-2*L[1])   )
par2<-c(exp(-L[2])*2*cos(2*pi*psi[2]), -exp(-2*L[2])   )
par3<-c(exp(-L[3])*2*cos(2*pi*psi[3]), -exp(-2*L[3])   )
freq = seq(.0005,.5-.0005, by=0.0005) 
W1<-ARMAspec(model=list(ar = par1 ), plot = T , freq = freq )
W2<-ARMAspec(model=list(ar = par2 ), plot = T , freq = freq )
W3<-ARMAspec(model=list(ar = par3 ), plot = T , freq = freq )

specmix<-phi[1]*W1$spec+phi[2]*W2$spec+phi[3]*W3$spec



## AR(12) model 

W<-ARMAspec(model=list(ar = c(0,0,0,0.9,0,0,0,.7,0,0,0,-.63) ), plot = F,freq = seq(0.002, 0.498, .002))
Wint<- simpson(W$spec,W$freq[1] )
Wstand<-W$spec/Wint






#### log curves #### 

## the next code assumes was run the extraction modes sources individually and obtain 
## all the sumaries, the code read the data bases sumaries and get the median curves and find the median IAE value 
## plot the ,atrix of curves and highlight the one with median IAE

#####  locating the median IAE ##### 

median_matrixBDPAR2<-matrix(rep(1:249, 1000), nrow = 1000 )
median_matrixBDPAR12<-matrix(rep(1:249, 1000), nrow = 1000 )
median_matrixBDPMA4<-matrix(rep(1:249, 1000), nrow = 1000 )
median_matrixCAR2<-matrix(rep(1:249, 1000), nrow = 1000 )
median_matrixCAR12<-matrix(rep(1:249, 998), nrow = 998 )
median_matrixCMA4<-matrix(rep(1:249, 1000), nrow = 1000 )
IAEAR2<-c()
IAEAR12<-c()
IAEMA4<-c()
ar2modes<- list()
mamodes<- list()
ar12modes<- list()
for(i in 1:1000){
  base1=readRDS(paste( "~/Desktop/AR2mixcomponents/compAR2mix_",i,".rds",sep=""  ) ) 
  ar2modes[[i]]<-base1
    median_matrixBDPAR2[i,]<- 2* base1$mediancurveglobal 
  IAEAR2[i]<- IAE(median_matrixBDPAR2[i,], specmix,W1$freq[1] )
}
for(i in 1:1000){
  base1=readRDS(paste( "~/Desktop/AR12components/compAR12_",i,".rds",sep=""  ) ) 
  ar12modes[[i]]<-base1
  median_matrixBDPAR12[i,]<- 2* base1$mediancurveglobal 
  IAEAR12[i]<- IAE(median_matrixBDPAR12[i,], Wstand,W1$freq[1] )
}
for(i in 1:1000){
  base1=readRDS(paste( "~/Desktop/MA4components/compMA4_",i,".rds",sep=""  ) ) 
  mamodes[[i]]<-base1
  median_matrixBDPMA4[i,]<- 2* base1$mediancurveglobal 
  IAEMA4[i]<- IAE(median_matrixBDPMA4[i,], WMA4st,W1$freq[1] )
}
for(i in 1:1000){
  base1=readRDS(paste( "~/Desktop/AR2mixchoud/meancurvechoudAR2_",i,".rds",sep=""  ) ) 
  median_matrixCAR2[i,]<- base1$globalmedian
  IAEAR2[i]<- IAE(median_matrixCAR2[i,], specmix,W1$freq[1] )
}

vecind<-c(1:1000)[-c(102,149)]  
for(i in 1:998){
  base1=readRDS(paste( "~/Desktop/AR12choud/meancurvechoudAR12_",vecind[i],".rds",sep=""  ) ) 
  median_matrixCAR12[i,]<- base1$globalmedian
  IAEAR12[i]<- IAE(median_matrixCAR12[i,], Wstand,W1$freq[1] )
}
for(i in 1:1000){
  base1=readRDS(paste( "~/Desktop/MA4choud/meancurvechoudMA4_",i,".rds",sep=""  ) ) 
  median_matrixCMA4[i,]<- base1$globalmedian
  IAEMA4[i]<- IAE(median_matrixCMA4[i,], WMA4st,W1$freq[1] )
}


# the index of the median value was found manually

#AR2 bmard  

matplot(seq(2,70,by=2), t(log(median_matrixBDPAR2[,1:35])), type="l", xlab = "Hz", ylab = "Standardized power", main="BMARD: AR(2) Mixture" )
lines(seq(2,70,by=2),log(specmix[1:35]), col=2, lwd=3)  
lines(seq(2,70,by=2),log(median_matrixBDPAR2[812,1:35]), type="l", lwd=3)
abline(v=base1$globalcenter[1,]*1000,col=3)
legend(x=25, y=4.8, legend=c("True log SSDF","log Median IAE BMARD",expression("MAP"~ bar(phi)~"Median BMARD" ) ), col=c(2,1,3), bty="n",lwd=c(3,3,1),cex=.8)

#AR12 bmard  

matplot(seq(2,498,by=2), t(log(median_matrixBDPAR12)), type="l", xlab = "Hz", ylab = "Standardized power", main="BMARD: AR(12) process" ,ylim=c(-5,6))
lines(seq(2,498,by=2),log(median_matrixBDPAR12[339,]), type="l", lwd=3)
lines(seq(2,498,by=2),log(Wstand), col=2, lwd=3)
abline(v=base1$globalcenter[1,]*1000,col=3)
legend(x=0, y=6, legend=c("True log SSDF","log Median IAE BMARD",expression("MAP"~ bar(phi)~"Median BMARD" ) ), col=c(2,1,3), bty="n",lwd=c(3,3,1),cex=.8)

#MA4 bmard  

matplot(seq(2,498,by=2), t(log(median_matrixBDPMA4 )), type="l", xlab = "Hz", ylab = "Standardized power", main="BMARD: MA(4) process", ylim=c(-3,3) )
lines(seq(2,498,by=2),log(WMA4st), col=2, lwd=3)
lines(seq(2,498,by=2),log(median_matrixBDPMA4[993,]), type="l", lwd=3)
abline(v=base1$globalcenter[1,]*1000,col=3)
legend(x=0, y=3, legend=c("True log SSDF","log Median IAE BMARD",expression("MAP"~ bar(phi)~"Median BMARD" ) ), col=c(2,1,3), bty="n",lwd=c(3,3,1),cex=.8)

#choud AR2  
matplot(seq(2,70,by=2), t(log(median_matrixCAR2[,1:35])), type="l", xlab = "Hz", ylab = "Standardized power", main="Bernstein: AR(2) Mixture" )
lines(seq(2,70,by=2),log(specmix[1:35]), col=2, lwd=3)
lines(seq(2,70,by=2),log(median_matrixCAR2[376,1:35]), type="l", lwd=3)
legend(x=25, y=4, legend=c("True log SSDF","log Median IAE Bernstein" ), col=c(2,1), bty="n",lwd=c(3,3),cex=.8)

#choud AR12  

matplot(seq(2,498,by=2), t(log(median_matrixCAR12)), type="l", xlab = "Hz", ylab = "Standardized power", main="Bernstein: AR(12) process",ylim=c(-5,6) )
lines(seq(2,498,by=2),log(median_matrixCAR12[879,]), type="l", lwd=3)
lines(seq(2,498,by=2),log(Wstand), col=2, lwd=3)

legend(x=0, y=6, legend=c("True log SSDF","log Median IAE Bernstein" ), col=c(2,1), bty="n",lwd=c(3,3),cex=.8)

#choud MA4  

matplot(seq(2,498,by=2), t(log(median_matrixCMA4)), type="l", xlab = "Hz", ylab = "Standardized power", main="Bernstein: MA(4) process",ylim=c(-3,3) )
lines(seq(2,498,by=2),log(WMA4st), col=2, lwd=3)
lines(seq(2,498,by=2),log(median_matrixCMA4[154,]), type="l", lwd=3)
legend(x=0, y=3, legend=c("True log SSDF","log Median IAE Bernstein" ), col=c(2,1), bty="n",lwd=c(3,3),cex=.8)

#kernel AR2  

# first found the IAE for the kernel estimates
IAEAR2<-c()
for(i in 1:1000){
  IAEAR2[i]<- IAE(kernelAR2matrix[i,], specmix,W1$freq[1] )
}
IAEAR2[1001]=-1
which(IAEAR2==median(IAEAR2))  
matplot(seq(2,70,by=2), t(log(kernelAR2matrix[,1:35])), type="l", xlab = "Hz", ylab = "Standardized power", main="Kernel: AR(2) Mixture" )
lines(seq(2,70,by=2),log(kernelAR2matrix[78,1:35]), type="l", lwd=3)
lines(seq(2,70,by=2),log(specmix[1:35]), col=2, lwd=3)
legend(x=35, y=5, legend=c("True log SSDF","log Median IAE Kernel" ), col=c(2,1), bty="n",lwd=c(3,3),cex=.8)

#kernel AR12  

IAEAR2<-c()
for(i in 1:1000){
  IAEAR2[i]<- IAE(kernelAR12matrix[i,], Wstand,W1$freq[1] )
}
IAEAR2[1001]=-1
which(IAEAR2==median(IAEAR2)) 
matplot(seq(2,498,by=2), t(log(kernelAR12matrix)), type="l", xlab = "Hz", ylab = "Standardized power", main="Kernel: AR(12) process", ylim=c(-8,6) )
lines(seq(2,498,by=2),log(Wstand), col=2, lwd=3)
lines(seq(2,498,by=2),log(kernelAR12matrix[967,]), type="l", lwd=3)
legend(x=0, y=6, legend=c("True log SSDF","log Median IAE Kernel" ), col=c(2,1), bty="n",lwd=c(3,3),cex=.8)

#kernel MA4  
IAEAR2<-c()
for(i in 1:1000){
  IAEAR2[i]<- IAE(kernelAR12matrix[i,], Wstand,W1$freq[1] )
}
IAEAR2[1001]=-1
which(IAEAR2==median(IAEAR2))   
matplot(seq(2,498,by=2), t(log(kernelMA4matrix)), type="l", xlab = "Hz", ylab = "Standardized power", main="Kernel: MA(4) process", ylim=c(-4,4) )
lines(seq(2,498,by=2),log(WMA4st), col=2, lwd=3)
lines(seq(2,498,by=2),log(kernelMA4matrix[967,]), type="l", lwd=3)
legend(x=0, y=4, legend=c("True log SSDF","log Median IAE Kernel" ), col=c(2,1), bty="n",lwd=c(3,3),cex=.8)


#splines AR2  

# first found the IAE for the spline estimates
IAEAR2<-c()
for(i in 1:1000){
  IAEAR2[i]<- IAE(splineAR2matrix[i,], specmix,W1$freq[1] )
}
IAEAR2[1001]=-1
which(IAEAR2==median(IAEAR2)) 
matplot(seq(2,70,by=2), t(log(splineAR2matrix[,1:35])), type="l", xlab = "Hz", ylab = "Standardized power", main="Cubic Spline: AR(2) Mixture" )
lines(seq(2,70,by=2),log(specmix[1:35]), col=2, lwd=3)
lines(seq(2,70,by=2),log(splineAR2matrix[659,1:35]), type="l", lwd=3)
legend(x=25, y=6, legend=c("True log SSDF","log Median IAE C. Spline" ), col=c(2,1), bty="n",lwd=c(3,3),cex=.8)

#splines AR12  
IAEAR2<-c()
for(i in 1:1000){
  IAEAR2[i]<- IAE(splineAR12matrix[i,], Wstand,.002 )
}
IAEAR2[1001]=-1
which(IAEAR2==median(IAEAR2)) 
matplot(seq(2,498,by=2), t(log(splineAR12matrix)), type="l", xlab = "Hz", ylab = "Standardized power", main="Cubic Spline: AR(12) process",ylim=c(-6,6)  )
lines(seq(2,498,by=2),log(splineAR12matrix[171,]), type="l", lwd=3)
lines(seq(2,498,by=2),log(Wstand), col=2, lwd=3)
legend(x=0, y=6, legend=c("True log SSDF","log Median IAE C. Spline" ), col=c(2,1), bty="n",lwd=c(3,3),cex=.8)

#splines AR12  
IAEAR2<-c()
for(i in 1:1000){
  IAEAR2[i]<- IAE(splineMA4matrix[i,], WMA4st,.002 )
}
IAEAR2[1001]=-1
which(IAEAR2==median(IAEAR2))  
matplot(seq(2,498,by=2), t(log(splineMA4matrix)), type="l", xlab = "Hz", ylab = "Standardized power", main="Cubic Spline: MA(4) process",ylim=c(-7,4)  )
lines(seq(2,498,by=2),log(WMA4st), col=2, lwd=3)
lines(seq(2,498,by=2),log(splineMA4matrix[522,]), type="l", lwd=3)
legend(x=0, y=4, legend=c("True log SSDF","log Median IAE C. Spline" ), col=c(2,1), bty="n",lwd=c(3,3),cex=.8)



####  AR2 mixture Disparity computation ####
# computation of the percentage of correct number of chains finding the true # of peaks
# conditional on true peaks; a how far are the modes using abs()
# last step individual per chain, and also in summamry as average or sum, 
AVG_phasedisparityv2<-matrix(rep(NA,3*1000),nrow = 1000)
AVG_BWdisparityv2<-matrix(rep(NA,3*1000),nrow = 1000)
AVG_weightdisparityv2<-matrix(rep(NA,3*1000),nrow = 1000)
correct_Ncomponents<-c()
maxcomp<-c()
chains=6
for(k in 1:1000){
  lens<-which( apply(ar2modes[[k]]$psimatrix , 1 ,sum_nna_limit, lim=.1  )==3 )
  correct_Ncomponents<-c(correct_Ncomponents,  length( lens )   )  
  modes_matrix<-c()
  modes_BWmatrix<-c()
  modes_Weightmatrix<-c()
  for(h in 1:chains ){
    if (  sum_nna_limit( ar2modes[[k]]$psimatrix[h,] ,  lim=.1  )  ==3  ){
      
      tempindex<-  which( ar2modes[[k]]$psimatrix[h,] <.1)
      tempvec<- ar2modes[[k]]$psimatrix[h,tempindex] 
      tempBW<- ar2modes[[k]]$BWmatrix[h,tempindex]
      tempWei= ar2modes[[k]]$weightmatrix[h,tempindex]
      modesorder<-tempvec[order(tempvec)  ] 
      BWorder<-tempBW[order(tempvec)  ] 
      weiorder<-tempWei[order(tempvec)  ] 
      if(h==1){
        modes_matrix<-modesorder*1000
        modes_BWmatrix<-BWorder
        modes_Weightmatrix<-weiorder
      }else{
        modes_matrix<-rbind(modes_matrix, modesorder*1000 )  
        modes_BWmatrix<- rbind(modes_BWmatrix, BWorder)
        modes_Weightmatrix<-rbind(modes_Weightmatrix,weiorder)
      }
    }
  }
  
  maxcomp[k]<-max(apply(ar2modes[[k]]$psimatrix , 1 ,sum_nna_limit, lim=.1  ) )
  
  if(length(lens)!=1){ 
    phasedisparitymat<-matrix( c( abs(modes_matrix[,1]-8),
                                  abs(modes_matrix[,2]-30),
                                  abs(modes_matrix[,3]-60)  ), ncol  = 3
    )
    BWdisparitymat<-matrix( c( abs(modes_BWmatrix[,1]-.03),
                               abs(modes_BWmatrix[,2]-.03),
                               abs(modes_BWmatrix[,3]-.03)  ), ncol = 3
    )
    Weightdisparitymat<-matrix( c( abs(modes_Weightmatrix[,1]-.1),
                                   abs(modes_Weightmatrix[,2]-.6),
                                   abs(modes_Weightmatrix[,3]-.3)  ), ncol = 3
    )
  }else{
    phasedisparitymat<-matrix( c( abs(modes_matrix[1]-8),
                                  abs(modes_matrix[2]-30),
                                  abs(modes_matrix[3]-60)  ), ncol  = 3
    )
    BWdisparitymat<-matrix( c( abs(modes_BWmatrix[1]-.03),
                               abs(modes_BWmatrix[2]-.03),
                               abs(modes_BWmatrix[3]-.03)  ), ncol = 3
    )
    Weightdisparitymat<-matrix( c( abs(modes_Weightmatrix[1]-.1),
                                   abs(modes_Weightmatrix[2]-.6),
                                   abs(modes_Weightmatrix[3]-.3)  ), ncol = 3
    )
  }
  
  # percentage will be given by the dimension of the average ohase disparity
  AVG_phasedisparityv2[k,]<-apply(phasedisparitymat, 2, mean)
  AVG_BWdisparityv2[k,]<-apply(BWdisparitymat, 2, mean)
  AVG_weightdisparityv2[k,]<-apply(Weightdisparitymat, 2, mean)
  
}# end for k


#report mean and SD changer cc from 1 to 3
cc=3
mean(AVG_phasedisparityv2[,cc] ,na.rm = T )
sd(AVG_phasedisparityv2[,cc],na.rm = T)

mean(AVG_BWdisparityv2[,cc],na.rm = T)
sd(AVG_BWdisparityv2[,cc],na.rm = T)
mean(AVG_weightdisparityv2[,cc],na.rm = T)
sd(AVG_weightdisparityv2[,cc],na.rm = T)




##### AR12 analysis of peaks disparity ####

AVG_phasedisparityar12<-matrix(rep(NA,5*1000),nrow = 1000)
correct_Ncomponents12<-c()
maxcompar12<-c()
for(k in 1:1000){
  
  lens<-which( apply(ar12modes[[k]]$psimatrix , 1 ,sum_nna  )==5 )
  correct_Ncomponents12<-c(correct_Ncomponents12,  length( lens )   ) 
  
  modes_matrix<-c()
  modes_BWmatrix<-c()
  modes_Weightmatrix<-c()
  for(h in 1:chains ){
    if (  sum_nna( ar12modes[[k]]$psimatrix[h,]   )  ==5  ){
      
      tempindex<-  1:5
      tempvec<- ar12modes[[k]]$psimatrix[h,tempindex] 
      tempBW<- ar12modes[[k]]$BWmatrix[h,tempindex]
      tempWei= ar12modes[[k]]$weightmatrix[h,tempindex]
      modesorder<-tempvec[order(tempvec)  ] 
      BWorder<-tempBW[order(tempvec)  ] 
      weiorder<-tempWei[order(tempvec)  ] 
      if(h==1){
        modes_matrix<-modesorder*1000
        modes_BWmatrix<-BWorder
        modes_Weightmatrix<-weiorder
      }else{
        modes_matrix<-rbind(modes_matrix, modesorder*1000 )  
        modes_BWmatrix<- rbind(modes_BWmatrix, BWorder)
        modes_Weightmatrix<-rbind(modes_Weightmatrix,weiorder)
      }
    }
  }  
  
  
  maxcompar12[k]<-max(apply(ar12modes[[k]]$psimatrix , 1 ,sum_nna  ))
  if(length(lens)!=0){

    if(length(lens)!=1){ 
      phasedisparitymat<-matrix( c( abs(modes_matrix[,1]-0),
                                    abs(modes_matrix[,2]-150),
                                    abs(modes_matrix[,3]-250),
                                    abs(modes_matrix[,4]-350),
                                    abs(modes_matrix[,5]-500)), ncol  = 5
      )
    }else{
      phasedisparitymat<-matrix( c( abs(modes_matrix[1]-0),
                                    abs(modes_matrix[2]-150),
                                    abs(modes_matrix[3]-250),
                                    abs(modes_matrix[4]-350),
                                    abs(modes_matrix[5]-500))
                                 , ncol  = 5
      )
    }

        AVG_phasedisparityar12[k,]<-apply(phasedisparitymat, 2, mean)
  }#end if
}# end for k


#report mean and SD
for( cc in 1:5){print( mean(AVG_phasedisparityar12[,cc], na.rm = T  ) )
  print( sd(AVG_phasedisparityar12[,cc], na.rm = T)
  )
  print("")
}


