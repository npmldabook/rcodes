### Chapter 8 Ex1 ###
rm(list = ls())
#----------------
library(npmlda)
#----------------
str(NGHS ) #19071*12
NGHS$Black <- (NGHS$RACE==2)*1
NGHS<- NGHS[!is.na(NGHS$SBP) & !is.na(NGHS$BMIPCT) & !is.na(NGHS$HTPCT ),]      
nrow(NGHS) #19320

Ct <-   data.frame(table(NGHS$ID))
names(Ct)<- c('ID', 'ni')
Ct$IDD<- 1:nrow(Ct)
NGHS<- merge(NGHS, Ct, by= 'ID')
nID<- nrow(Ct) #2376


NGHS$HTPCTc<- NGHS$HTPCT-50
NGHS$BMIPCTc<- NGHS$BMIPCT-50

# 100 agebins #
NGHS <- data.frame(NGHS, agebin=numeric(nrow(NGHS))) 

for (i in 9:18){
  print(i)
  for (j in 1:10){
    NGHS$agebin[( NGHS$AGE*10 >= (i*10+ (j-1)) &  NGHS$AGE*10 < (i*10+ j))]  <-   (i*10+ (j-1)) 
    print(c(i+ (j-1)/10))
  } }
summary(NGHS$agebin)

##1. Getting the raw coefficent at each agebin##
AgeCat <-   seq(90, 189, by= 1) 
Beta.raw <- matrix(NA, nrow=100, ncol=4)

for (j in 1:100)
{
  print(c(j, "Agebin", AgeCat[j]/10))
  Datai <- NGHS[(NGHS$agebin== AgeCat[j]) ,]
  Beta.raw[j,]<-  lm(SBP~  Black+HTPCTc +BMIPCTc , data=Datai )$coef 
}

##2. Getting the raw coefficent at each agebin##

Var1 <-  Beta.raw[,1] 

hh<- seq(0.2, 5,by=0.1)

Nh<- length(hh) 
CVh<- numeric(Nh)
N<-100

 # leaving one time-point CV#
for (j in 1:Nh)
{
  #print(j)
  CVh[j] <- CVlm(AgeCat/10, Var1 , bw=hh[j], ID=1:N, Wt=rep(1/N,N))   
}
plot(hh, CVh, type='o')
hh[CVh==min(CVh)]  #1.3
# similarly get CV hh for other coefficients #


Beta.f0 <- LocalLm(AgeCat/10, AgeCat/10, Beta.raw[,1], bw =1.3 )  
Beta.f1 <- LocalLm(AgeCat/10, AgeCat/10, Beta.raw[,2], bw =4 )    
Beta.f2 <- LocalLm(AgeCat/10, AgeCat/10, Beta.raw[,3], bw =2.4 )  
Beta.f3 <- LocalLm(AgeCat/10, AgeCat/10, Beta.raw[,4], bw =1.5 )  

#Bootstrap sample is generated similar as previous sections##

N <- length(unique(NGHS$ID)) #2376
NGHS2<- NGHS[,c('ID', 'agebin', 'Black','HTPCTc','BMIPCTc','SBP') ]

Bootsample <- function(){ 
  resample.ID <- sample(x= unique(NGHS$ID) ,size= N ,replace=T) 
  do.call("rbind", lapply(1:N , function(i) data.frame(subset(NGHS2, ID==resample.ID[i]),ID2=i ) ))
}



