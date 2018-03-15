### Chapter 15 Ex ###
rm(list = ls())
#----------------
library(npmlda)
library(lme4)
library(splines)

#----------------
# read CDC quantiles
BMIq85<- c(19.11937,19.26034,19.40282,19.61895,19.76436,19.98400,20.13118,20.27876,20.50052,20.64838,20.86984,21.01703,21.16371,21.38246,21.52727
            ,21.74263,21.88480,22.02573,22.23458,22.37196,22.57506,22.70835,22.83987,23.03366,23.16045,23.34689,23.46861,23.58823,23.76363,23.87784
            ,24.04503,24.15370,24.26015,24.41564,24.51653,24.66372,24.75913,24.85240,24.98836,25.07643,25.20482,25.28802,25.36937,25.48810,25.56517
            ,25.67786,25.75118,25.82317,25.92887,25.99799)
#----------------
str(NGHS ) #19071*12

NGHS$Black <- (NGHS$RACE==2)*1
NGHS<- NGHS[!is.na(NGHS$BMI),]
dim(NGHS) #19398 obs

NGHS.B <- NGHS[NGHS$Black ==1,]
nID <- length(unique(NGHS.B$ID)) #1213

# fit the a linear mixed model with cubic B-splines
KN1 <- seq(from=9, to=18.9, length=4)[-c(1,4)]
Bs.age <- bs(NGHS.B$AGE, knots=KN1)
fm1 <- lmer(BMI ~ 1+ Bs.age +(1+ Bs.age|ID), data=NGHS.B)

# generate the BLUP predictions on a grid
T2<- seq(from=9, to=18.8, by=0.2)
S.BS <- bs(T2, knots=KN1 )
IDlevel<- rownames(coef(fm1)[[1]])
NGHSped<- data.frame(ID= NA, AGE=rep(T2, nID), BMI=NA)
nT2<- 50
for (i in 1:nID)
{
  KK <- i-1
  Datai <- NGHS.B[NGHS.B$ID==IDlevel[i],]
  mean.pred <- cbind(1,S.BS) %*% t(as.vector(coef(fm1)[[1]][i,]))
  NGHSped[(KK*nT2+1):(KK*nT2+nT2),]$ID <- IDlevel[i]
  NGHSped[(KK*nT2+1):(KK*nT2+nT2),]$BMI <- mean.pred
}
NGHSped$BMIp <- NGHSped$BMI + rnorm(nrow(NGHSped) ,
                                      mean = 0, sd = sigma(fm1) )

###########

NN<- nrow(NGHSped)

S1cat <- seq(9.0, 18.8, by=0.2)
S12cat <- seq(9.0, 15.8, by=0.2)
KK1 <- length(S1cat)
KK2 <- length(S12cat)

Prob.S1 <- numeric(KK1 )
for (i in 1:KK1)
{
  SeqKK1 <- seq(from=i, to =NN-nT2+i, by=nT2)
  Datai <- NGHSped[SeqKK1,]
  Prob.S1[i] <- mean(Datai$BMIp>=BMIq85[i])
  # BMIq85 is the CDC percentile curve
}


## Compute the PA(x,s1,s2): grid*15=3 year##
Prob.S1S2 <- numeric( KK2 )
for (i in 1:KK2)
{ print(i)
  SeqKK1 <- seq(from=i, to =NN-nT2+i, by=nT2)
  Datai <- cbind(NGHSped[SeqKK1,]$BMIp, NGHSped[SeqKK1+15,]$BMIp)
  Prob.S1S2[i]<- mean((Datai[,1]>=BMIq85[i])&
                        (Datai[,2]>=BMIq85[i+15]))
}

RTP <- Prob.S1S2/Prob.S1[S1cat <= 15.8 ]
RTPR<- RTP/(Prob.S1[S1cat   >= 12.0 ]) 

summary(RTP)
summary(RTPR)



