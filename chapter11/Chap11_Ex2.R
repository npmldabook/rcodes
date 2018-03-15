### Chapter 11 Ex2 ###
rm(list = ls())
#----------------
library(npmlda)
library(splines)
library(lme4)

## starting here for the chapter 11 #3
str(NGHS ) #19071*12
NGHS$Black <- (NGHS$RACE==2)*1
NGHS<- NGHS[!is.na(NGHS$SBP) & !is.na(NGHS$HTPCT ),]      
nrow(NGHS) #19439


KN1 <- seq(from=9, to=19, length=5)[-c(1,5)]
Bs.age <- bs(NGHS$AGE, knots=KN1)
fm.Ht <- lmer(HTPCT ~ 1+ Bs.age +(1+ Bs.age|ID), data=NGHS)
NGHS$Htfitted <- fitted(fm.Ht ) #BLUP estimate for height

#--------fig 11.3:BLUP for height -------- #

T.range<- range(NGHS$AGE)  
Tgrid2 <- seq(from=T.range[1], to=T.range[2], length=20)  
BS <-  bs(Tgrid2, knots=KN1 , degree=3,  intercept=F)  
mean.hat <-  cbind(1,BS) %*% fixef(fm.Ht)
IDlevel<- rownames(coef( fm.Ht )[[1]])

postscript("fig11-3-new.ps", horizontal=T)

 par(mfrow=c(2,3), mar=c(5,4.1,2.5,1), cex.lab=1.5, cex.axis=1.5, cex=0.8, cex.main=1.5, font.main=1)

 k<-1
 for (i in c(3,14, 16,9,20,15))
 {   
   Datai <- NGHS[NGHS$ID==IDlevel[i],]
   plot(Datai$AGE, Datai$HTPCT, ylim=c(0,100), pch=1, main=paste("Girl",k),
       xlab='Age', ylab="Height(%)", xlim=c(9,19) )
   lines( Datai$AGE,   Datai$Htfitted, col=1,lty=1)
   lines(Tgrid2 ,mean.hat, lty=2)
  
   mtext(paste(LETTERS[k],".",sep =""), cex=1.3, side=3, line=0.7, font=2, at=6.3)
   k<-k+1
 }

dev.off()
#-------------------------------

KN2 <- seq(from=9, to=19, length=4)[-c(1,4)]
Bs.Age <- bs(NGHS$AGE, knots=KN2, intercept=T)
Bs.Race <- bs(NGHS$AGE, knots=NULL,intercept=T)*NGHS$Black
Bs.Ht <- bs(NGHS$AGE, knots=NULL,intercept=T)*(NGHS$Htfitted-50)
fm.SBP <- lmer(SBP~ 0+ Bs.Age+ Bs.Race+ Bs.Ht+(0+ Bs.Age|ID),data=NGHS)

Tgrid<- seq(from=9, to=19, by=0.5)
nn   <- cumsum(c(ncol(Bs.Age ),ncol(Bs.Race),ncol(Bs.Ht)))
Bhat0b <- bs(Tgrid, knots=KN2,  intercept=T) %*% fixef(fm.SBP)[1:nn[1]]
Bhat1b <- bs(Tgrid, knots=NULL, intercept=T) %*% fixef(fm.SBP)[(nn[1]+1):nn[2]]
Bhat2b <- bs(Tgrid, knots=NULL, intercept=T) %*% fixef(fm.SBP)[(nn[2]+1):nn[3]]

SD.b0 <- SD.b1 <- SD.b2 <- numeric(length(Tgrid) )
nn1<- nn[1]
nn2<- nn[2]
nn3<- nn[3]

cov1<-  vcov(fm.SBP)[1:nn1,1:nn1] 
cov2<-  vcov(fm.SBP)[(nn1+1):nn2, (nn1+1):nn2] 
cov3<-  vcov(fm.SBP)[(nn2+1):nn3, (nn2+1):nn3] 

for (i in 1:length(Tgrid))
{
  print(i)
  Bs.agei <- bs(Tgrid[i], knots=KN2, degree=3,  intercept=T, Boundary.knots=c(9, 19)) 
  Bs2i  <-  bs(Tgrid[i], knots=NULL, degree=3,  intercept=T,Boundary.knots=c(9, 19))  #race 
  Bs3i  <-  bs(Tgrid[i], knots=NULL,  degree=3,  intercept=T,Boundary.knots=c(9, 19)) #HT 
  
  SD.b0[i] <- (t(as.numeric(Bs.agei)) %*%cov1 %*% as.numeric(Bs.agei))[1,1]
  SD.b1[i] <- (t(as.numeric(Bs2i))    %*%cov2 %*% as.numeric( Bs2i ))[1,1]
  SD.b2[i] <- (t(as.numeric(Bs3i))    %*%cov3 %*% as.numeric( Bs3i ))[1,1]
}

Bhat0b.LCL<- Bhat0b -  1.96*sqrt(SD.b0)
Bhat0b.UCL<- Bhat0b +  1.96*sqrt(SD.b0)

Bhat1b.LCL<- Bhat1b -  1.96*sqrt(SD.b1)
Bhat1b.UCL<- Bhat1b +  1.96*sqrt(SD.b1)

Bhat2b.LCL<- Bhat2b*20 -  1.96*sqrt(SD.b2)*20
Bhat2b.UCL<- Bhat2b*20 +  1.96*sqrt(SD.b2)*20

#------------------------------------------------#
postscript("fig11-4.ps", horizontal=T)

par(mfrow=c(1,3), mar=c(30.5,4.5,2.5,1),  cex.lab=1.5, cex.axis=1.5, cex=0.8, cex.main=1.5 , font.main=1)

plot(Tgrid, Bhat0b, ylab=expression(paste(beta, "0(t)")), xlab="Age", main= " Baseline ", type="l", ylim=c(95, 110))
lines(Tgrid, Bhat0b.LCL, lty=2)
lines(Tgrid, Bhat0b.UCL, lty=2)
mtext("A", cex=1.3, side=3, line=0.7, font=2, at= 6.3)

plot(Tgrid, Bhat1b, ylab=expression(paste(beta, "1(t)")),xlab="Age",main="Race effect", type="l", ylim=c(-0.6, 3.6) )
lines(Tgrid, Bhat1b.LCL, lty=2)
lines(Tgrid, Bhat1b.UCL, lty=2)
mtext("B.", cex=1.3, side=3, line=0.7, font=2, at=6.3)

plot(Tgrid, Bhat2b*20, ylab=expression(paste(beta, "2(t)")),xlab="Age", main="Height effect per 20% increase", type="l", ylim=c(-0.4, 3.4))
lines(Tgrid, Bhat2b.LCL, lty=2)
lines(Tgrid, Bhat2b.UCL, lty=2)
mtext("C.", cex=1.3, side=3, line=0.7, font=2, at=6.3)

dev.off()


## predict SBP :BLUP   ####
NGHS$SBPfitted <- fitted(fm.SBP ) #BLUP estimate for SBP

#------------------------------------------------#
main115<- c('Caucasian girl with 11% height', 'Caucasian girl with 79% height', 'African American girl with 25% height')

postscript("fig11-5.ps", horizontal=T)

par(mfrow=c(1,3), mar=c(30.5,4.5,2.5,1),  cex.lab=1.5, cex.axis=1.5, cex=0.8, cex.main=1.25 , font.main=1)
k<-1
for (i in c(124,210, 558))
{
  Datai <- NGHS[NGHS$ID==IDlevel[i],]
  plot(Datai$AGE, Datai$SBP, ylim=c(85,130), pch=1, main=main115[k],  xlim=c(9,19), xlab="Age", ylab="SBP" )
  
  Bs1 <-   bs(Datai$AGE, knots=KN2,  intercept=T)
  Bs2  <-  bs(Datai$AGE, knots=NULL, intercept=T)  
  Bs3  <-  bs(Datai$AGE, knots=NULL, intercept=T) 
 
  Bhat0 <-   Bs1  %*% (fixef(fm.SBP)[1:nn[1]])
  Bhat1 <-   Bs2  %*% (fixef(fm.SBP)[(nn[1]+1):nn[2]])
  Bhat2 <-   Bs3  %*% (fixef(fm.SBP)[(nn[2]+1):nn[3]])
  
  HTpcti <-  Datai[Datai$ID==IDlevel[i],]$HTped
  
  mean.hat.SBP <-   Bhat0 + Bhat1* Datai$Black[1] + Bhat2*(Datai$Htfitted-50) 
  
  lines( Datai$AGE, mean.hat.SBP, col=1)
  lines( Datai$AGE, Datai$SBPfitted  , col=1, lty=2)
  
  mtext(paste(LETTERS[k],".",sep =""), cex=1.3, side=3, line=0.7, font=2, at=6.3)
  k<-k+1
  
  }

dev.off()



















