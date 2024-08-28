### computation of all modes#### 
i<-as.numeric(  commandArgs()[6])


# change the directory accordingly
source( "~/BMARD_code/Rcode/Extractionmaincomponentsmodes.R" )

listmodes<-list()

base<-readRDS(paste( "~/BMARD_code/Rcode/BayesiansimulMA40720_V",i,".rds",sep=""  )  )
  
listmodes<- modes_extraction_gaussianmix(data=base,threshold_weight=.005, chains=6,Nsamp=1e5, freq_rate=500, quant=.9,percsamp=.6 ,qlow=.025, qsup=.975  )


# change the directory accordingly
saveRDS(listmodes,paste( "~/BMARD_code/Rcode/compMA4_",i,".rds",sep=""  ) )


