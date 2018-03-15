### Chapter 13 Ex1: transformation model ###
rm(list = ls())
#----------------
library(npmlda)
library(survival)

#str(NGHS ) #19071*12

NGHS <- data.frame(NGHS, agebin=numeric(nrow(NGHS))) 

for (i in 9:18){
  print(i)
  for (j in 1:10){
    NGHS$agebin[( NGHS$AGE*10 >= (i*10+ (j-1)) &  NGHS$AGE*10 < (i*10+ j))]  <-   (i*10+ (j-1)) 
    print(c(i+ (j-1)/10))
  } }

#summary(NGHS$agebin)

#####################################################################################
#---- Time-varying transformation (proportional odds) model with 2 covariate -- #
TVtrans.fit<- function( Agebins,  Y=SBP, X1= (RACE==2)*1, X2=HTPCT)
{
  
  nbins<- length(Agebins)
  OR<- matrix(NA, nrow=nbins, ncol=2)
  
  for (i in 1:nbins)
  { 
    print(c(i, "Agebin", Agebins[i]))
    Datai <- data.frame(Y,  X1, X2, censor=1)[agebin== Agebins[i], ]
    
    FIT1<- survreg(Surv( Y, censor) ~  X1  +X2, data=Datai ,   dist="loglogistic") 
    
    
    Tvec  <-  Datai$Y; 
    nT<- length(Tvec)
    
    Tmatrix <- Tvec%*%t(rep(1, nT))  
    Indicator <- (Tmatrix - t(Tmatrix))>=0   #Indicator(i,j) is Ti > Tj
    
    # get 2 matrix Zi-Zj
    
    Z1vec    <-  Datai$X1 ;
    Z1matrix <-  Z1vec%*%t(rep(1, nT))
    Zij1     <-  as.vector(Z1matrix - t(Z1matrix))  #Z(ij)(1) = Z(i) - Z(j)
    
    
    Z2vec    <-  Datai$X2;
    Z2matrix <-  Z2vec %*%t(rep(1, nT))
    Zij2     <-  as.vector(Z2matrix - t(Z2matrix))  #Z(ij)(2) = Z(i) - Z(j)
    
    Zij<- cbind(Zij1, Zij2 ) 
    
    ## 2.1. calculate the U(bo) and U'(b0), initial b0
    #initial value:
    b0<- as.vector((-FIT1$coef/FIT1$scale)[2:3])
    
    #Three estimation function
    Diff<-  as.vector(Indicator)*1-  Xi(Zij%*%b0);  #nT*nT*1
    
    U1  <-  sum(Zij1* Diff)
    U2  <-  sum(Zij2* Diff)
    
    Ub <- c(U1, U2) 
    
    # Now write a function for the Newton solve:
    OR[i,]<- Newton2var(Zij, b0, Ub, Indicator, difflmt=1e-12)
  }
  OR
}
##############################################################################

NGHS.sbp <- NGHS[!is.na(NGHS$SBP) & !is.na(NGHS$HTPCT ),] 
dim(NGHS.sbp) #19439 *13
attach(NGHS.sbp)

Agebins <- seq(90, 189, by= 1) #100 bins
OR.SBP <- TVtrans.fit(Agebins,  Y=SBP, X1= (RACE==2)*1, X2=HTPCT)

SBP.beta1.lm <-LocalLm(Agebins/10, Agebins/10, OR.SBP[,1], bw=1.6)
SBP.beta2.lm <-LocalLm(Agebins/10, Agebins/10, OR.SBP[,2], bw=2.4)
detach(NGHS.sbp)

#----------------------------------------------------------------------------#
NGHS.dbp <- NGHS[!is.na(NGHS$DBP) & !is.na(NGHS$HTPCT ),]      
dim(NGHS.dbp)#19404*13
attach(NGHS.dbp)

OR.DBP <- TVtrans.fit(Agebins,  Y=DBP, X1= (RACE==2)*1, X2=HTPCT)

DBP.beta1.lm <-LocalLm(Agebins/10, Agebins/10, OR.DBP[,1], bw=2.4)
DBP.beta2.lm <-LocalLm(Agebins/10, Agebins/10, OR.DBP[,2], bw=2.4)
detach(NGHS.dbp)

# Generate resampling bootstrap by the same approach from earlier chapters #














