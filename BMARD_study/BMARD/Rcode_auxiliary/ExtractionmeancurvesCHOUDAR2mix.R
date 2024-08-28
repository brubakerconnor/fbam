
### informaiton extraction from choudhuri, 

# 1 mean curve
# 2 median curve
# 3 quantile curves, 
X<-as.numeric(  commandArgs()[6])

# change the directory accordingly
dbase<-readRDS( paste("~/BMARD_code/Rcode/ChoudhurimixAR2_V",X,".rds", sep="")   )

chains=6
nsamp=100000

processchoud<-function(dbase, chains, nsamp ){
matrixall<-matrix( rep(1:249,chains*5000 ), nrow=chains*5000)
meancurvematrix=matrix( rep(1:249,chains ), nrow=chains)
mediancurvematrix=matrix( rep(1:249,chains ), nrow=chains)
quppermatrix=matrix( rep(1:249,chains ), nrow=chains)
qlowermatrix=matrix( rep(1:249,chains ), nrow=chains)

#take stats global and matrix of individuals

for(k in 1: chains){
  matrixall[ (1+5000*(k-1))  : (5000*k) , ]=dbase[[k]]$CurveEsts[ (nsamp-5000+1):nsamp ,]
  meancurvematrix[k,]= apply(  dbase [[k]]$CurveEsts[ (nsamp-5000+1):nsamp ,],2,mean ) 
  mediancurvematrix[k,]=apply(  dbase [[k]]$CurveEsts[ (nsamp-5000+1):nsamp ,],2,median ) 
  quppermatrix[k,]=apply(  dbase [[k]]$CurveEsts[ (nsamp-5000+1):nsamp ,],2,quantile, .975 ) 
  qlowermatrix[k,]=apply(  dbase [[k]]$CurveEsts[ (nsamp-5000+1):nsamp ,],2,quantile, .025 ) 
}

globalmean<-apply(  matrixall,2,mean ) 
globalmedian=apply(  matrixall,2,median ) 
globalqupper=apply(  matrixall,2,quantile, .975 ) 
globalqlower=apply(  matrixall,2,quantile, .025)
  
return(list(
  matrixall=matrixall,
  meancurvematrix=meancurvematrix,
  mediancurvematrix=mediancurvematrix,
  quppermatrix=quppermatrix,
  qlowermatrix=qlowermatrix,
  globalmean=globalmean,
  globalmedian=globalmedian,
  globalqupper=globalqupper,
  globalqlower=globalqlower
) )  
  
}  

T_T<-processchoud(dbase,6,100000)


# change the directory accordingly
saveRDS(T_T,paste( "~/BMARD_code/Rcode/meancurvechoudAR2_",X,".rds",sep=""  ) )





  
