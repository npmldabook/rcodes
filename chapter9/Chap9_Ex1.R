### Chapter 9 Ex1 ###
rm(list = ls())
#----------------
library(npmlda)
library(splines)

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

# center variables#
NGHS$HTPCTc<- NGHS$HTPCT-50
NGHS$BMIPCTc<- NGHS$BMIPCT-50

# set up knots
nknots <- c(4,0,0,0)+2 

KN0 <- seq(from=9, to=19, length=nknots[1])[-c(1,nknots[1])]
KN1 <- seq(from=9, to=19, length=nknots[2])[-c(1,nknots[2])]
KN2 <- seq(from=9, to=19, length=nknots[3])[-c(1,nknots[3])]
KN3 <- seq(from=9, to=19, length=nknots[4])[-c(1,nknots[4])]

# generate cubic B-spline regression basis
 Bs.age <- bs(NGHS$AGE, knots=KN0, intercept=T)
 Bs.race <- bs(NGHS$AGE, knots=KN1, intercept=T)* NGHS$Black
 Bs.HT <- bs(NGHS$AGE, knots=KN2, intercept=T)* NGHS$HTPCTc
 Bs.BMI <- bs(NGHS$AGE, knots=KN3, intercept=T)* NGHS$BMIPCTc

 fit.WLS <- lm(SBP ~ 0 + Bs.age + Bs.race + Bs.HT + Bs.BMI,
                 weights=1/ni, data=NGHS)
 summary(fit.WLS)

##obtain the coefficient curves##
 tgrid <- seq(from=9, to=19, by=0.1)
 Bs0 <- bs(tgrid, knots=KN0, intercept=T)
 Bs1 <- bs(tgrid, knots=KN1, intercept=T)
 Bs2 <- bs(tgrid, knots=KN2, intercept=T)
 Bs3 <- bs(tgrid, knots=KN3, intercept=T)

  nn <- cumsum(c(ncol(Bs0), ncol(Bs1), ncol(Bs2), ncol(Bs3)))
 Beta0 <- Bs0 %*% fit.WLS$coef[1:nn[1]]
 Beta1 <- Bs1 %*% fit.WLS$coef[(nn[1]+1):nn[2]]
 Beta2 <- Bs2 %*% fit.WLS$coef[(nn[2]+1):nn[3]]
 Beta3 <- Bs3 %*% fit.WLS$coef[(nn[3]+1):nn[4]]

 # not just plot beta0-beta3, no random curve like SII paper since did not assume all coefficient are random 
 par(mfrow=c(2,2))
  plot(tgrid, Beta0, ylab=expression(paste(beta, "0(t)")),  main= " Baseline ", type="l", ylim=c(95,110))
  plot(tgrid, Beta1, ylab=expression(paste(beta, "1(t)")),  main= " Race ", type="l", ylim=c(-2,3))
  plot(tgrid, Beta2, ylab=expression(paste(beta, "2(t)")),  main= " Height percentile", type="l",ylim=c(0,0.1))
  plot(tgrid, Beta3, ylab=expression(paste(beta, "3(t)")),  main= " BMI percentile", type="l",ylim=c(0, 0.15))
 
 # similar bootstrap step to get confidence interval as chapter 8#
 
 
 
 
 
 
 
 