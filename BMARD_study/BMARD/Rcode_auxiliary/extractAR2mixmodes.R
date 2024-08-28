### computation of all modes#### 
i<-as.numeric(  commandArgs()[6])
#setwd("/ibex/scratch/projects/c2088/")
source( "~/BMARD_code/Rcode/Extractionmaincomponentsmodes.R" )

listmodes<-list()

#for(i in 1:1000){

base<-readRDS(paste( "~/BMARD_code/Rcode/BayesiansimmixAR2100k_V",i,".rds",sep=""  )  )
  
listmodes<- modes_extraction_gaussianmix(data=base,threshold_weight=.005, chains=6,Nsamp=1e5, freq_rate=500, quant=.9,percsamp=.6 ,qlow=.025, qsup=.975  )
  
#}

saveRDS(listmodes,paste( "~/BMARD_code/Rcode/compAR2mix_",i,".rds",sep=""  ) )

#print("I'm done")
