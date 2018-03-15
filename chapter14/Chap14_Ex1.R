### Chapter 14 Ex1: transformation model ###
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

# summary(NGHS$agebin)

NGHS.sbp <- NGHS[!is.na(NGHS$SBP) & !is.na(NGHS$HTPCT ),] 
dim(NGHS.sbp) #19439 *13

attach(NGHS.sbp)
Agebins <- seq(90, 189, by= 1) #100 bins

##############################################################################
# Conditional probability estimate based on Time-varying transformation 
#  (proportional odds) model with 2 covariate 
#----------------------------------------------------------------------------

Cond.Prob <- function(Agebins, Y=SBP, X1= (RACE==2)*1, X2=HTPCT, Y0=100, X10=1, X20=50)
{
  nbins<- length(Agebins)
  OR<- matrix(NA, nrow=nbins, ncol=2)
  hy0 <- Prob<-  rep(NA, nbins )
  
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
    OR[i,] <- Newton2var(Zij, b0, Ub, Indicator, difflmt=1e-12)
    
    #----------------
    h0<- 0
    Z12vec  <- cbind(Z1vec   , Z2vec)    
    HZB <- h0 + Z12vec%*% OR[i,]
    
    Ind.Y <-  sum(Datai$Y > Y0)
    Vh <- Ind.Y - sum(1/(1+exp(HZB)))  
    
    hy0[i] <-  Newton1var(Z12vec, h0 , Vh ,HZB, Ind.Y , Diff=1e-8,  ORR= OR[i,], MaxIter=100)
    
    Prob[i] <- 1/( 1+exp(hy0[i]+ X10* OR[i,1] + X20*OR[i,2] ) ) 
  }
  Prob
}

##############################################################################

# Get raw estimates of the conditional probability for given X10=1 (black),X20=50%
Prob.Y100 <- Cond.Prob(Agebins, Y=SBP, X1= (RACE==2)*1, X2=HTPCT, Y0=100, X10=1, X20=50)

lines( Agebins/10, Prob.Y100.lm , col='dark green', lwd=2 )
# Local linear smoothing estimate
Prob.Y100.lm <- LocalLm(Agebins, Agebins, Prob.Y100, bw=16)






