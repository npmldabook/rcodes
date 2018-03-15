### Chapter 7 Ex1 ###
rm(list = ls())
#----------------
library(npmlda)
#----------------
str(NGHS ) #19071*12
NGHS$Black <- (NGHS$RACE==2)*1
NGHS<- NGHS[!is.na(NGHS$SBP) & !is.na(NGHS$BMIPCT) & !is.na(NGHS$HTPCT ),]      
nrow(NGHS) #19320

Ct <-   data.frame(table(NGHS$ID),  1:nrow(Ct))
names(Ct)<- c('ID', 'ni','IDD')
NGHS<- merge(NGHS, Ct, by= 'ID')
nID<- nrow(Ct) #2376


Age.grid <- seq(9, 19, by=0.5) #21

NGHS$HTPCTc<- NGHS$HTPCT-50
NGHS$BMIPCTc<- NGHS$BMIPCT-50

Beta <- with(NGHS, LocalLm.Beta(Age.grid, AGE, X1=Black,  X2=HTPCTc, X3=BMIPCTc, SBP, Bndwdth=3.5, Weight=1/ni))


# Bootstrap CI ##
NGHS2<- NGHS[,c('ID', 'AGE', 'Black','HTPCTc','BMIPCTc','SBP','ni') ]

Bootsample <- function(){ 
  resample.ID <- sample(x= unique(NGHS$ID) ,size= nID ,replace=T) 
  do.call("rbind", lapply(1:nID , function(i) data.frame(subset(NGHS2, ID==resample.ID[i]),ID2=i ) ))
}

#######
AgeBins<- seq(9.3, 18.7, by=.2)  

NN<- length(AgeBins); NN

nBoot<-  1000
Boot.FIT1  <- matrix(NA, nrow= nBoot, ncol= NN)
Boot.FIT2 <- matrix(NA, nrow= nBoot, ncol= NN)
Boot.FIT3 <- matrix(NA, nrow= nBoot, ncol= NN)
Boot.FIT4  <- matrix(NA, nrow= nBoot, ncol= NN)

set.seed(101)

for (j in 1:nBoot)
{
  if ((j-floor(j/20)*20)== 0 )  print(j) 
  Bootdata <- Bootsample()
  
  Beta <-  LocalLm.Beta(AgeBins, Tvec=Bootdata$AGE, X1=Bootdata$Black, X2=Bootdata$HTPCTc, 
                            X3= Bootdata$BMIPCTc, Yvec=Bootdata$SBP, Bndwdth=3.5,  Weight = 1/Bootdata$ni)
  
  Boot.FIT1[j,] <- Beta[,1] 
  Boot.FIT2[j,] <- Beta[,2]
  Boot.FIT3[j,] <- Beta[,3]
  Boot.FIT4[j,] <- Beta[,4]
  
}


## CV choice h=1.6 ##
Beta.CV <- with(NGHS, LocalLm.Beta(AgeBins, AGE, X1=Black,
                                   X2=HTPCTc, X3=BMIPCTc, SBP, Bndwdth=1.6, Weight=1/ni))

## CV function for choosing bandwidth of LS LocalLm fit for NGHS ##
CV.LocalLm.Beta <- function(h, nID)
{
  
    Yest <- unlist(lapply(1:nID,  function(i) { 
    #if ((i-floor(i/100)*100)== 0 ) print(i) 
    DataRmi<- NGHS[NGHS$IDD!=i,] 
    Datai<- NGHS[NGHS$IDD ==i, c('AGE', 'Black', 'HTPCTc','BMIPCTc' )] 
    
    tvec1 <- Datai$AGE
    Betai <-  LocalLm.Beta(tvec1, Tvec=DataRmi$AGE, X1=DataRmi$Black, X2=DataRmi$HTPCTc, X3= DataRmi$BMIPCTc, 
                           Yvec= DataRmi$SBP, Bndwdth=h,  Weight = 1/DataRmi$ni)
    
    DesignX <- as.matrix(data.frame(1, Datai[,2:4] ))
    as.numeric(apply(DesignX* Betai , 1, sum))
   }
  ))
  
  Yvec <- NGHS$SBP
  Weight<- 1/NGHS$ni
  sum( (Weight)*(Yvec- Yest)^2)/nID  
}



CV.LocalLm.Beta(h=3.5, nID=2376) #70.38393
CV.LocalLm.Beta(h=1.6, nID=2376) #70.2622











































