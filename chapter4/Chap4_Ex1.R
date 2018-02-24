### Chapter 4 Ex1 ###
# remove (almost) everything in the working environment.
rm(list = ls())
#----------------

library(npmlda)
library(splines)

HSCT<- HSCT[!is.na(HSCT$Granu),]
str(HSCT)
summary(HSCT$Granu)

#------
Ct <- data.frame(table(HSCT$ID))
names(Ct)<- c("ID", "ni")
HSCT<- merge(HSCT, Ct, by= "ID")

## Code for Figure 4.1 analysis ##
 attach(HSCT)
 Granu.log <- log10(Granu)

 # 2 knots ,cubit 
 ( KN2 <- quantile(Days, c(.33, .66)))
 bs.Days<- bs(Days, knots=KN2, degree=3)

 # Obtain coefficients for the spline basis, subject-uniform weights
 Spline.fit <- lm(Granu.log ~ bs.Days, weights=1/ni)
# Obtain fitted estimates for a given x
 New.Days<- bs(-7:35, knots=KN2, degree=3)
 Spline.Est <- cbind(1, New.Days) %*% coef(Spline.fit)

#---measurement uniform
 Spline.fit0 <- lm(Granu.log ~ bs.Days)
 Spline.Est0 <- cbind(1, New.Days) %*% coef(Spline.fit0)
 
##1 knot, linear 

 (KN1 <- median(HSCT$Days)) #8
 bs.DAY.L1 <- bs(HSCT$Days, knots= KN1 , degree=1) 
 Fit.1L   <- lm(Granu.log ~ bs.DAY.L1)
 Fit.1L.W <- lm(Granu.log ~ bs.DAY.L1, weights=1/ni)

 New.Days.L1 <- bs(-7:35, knots=KN1, degree=1)
 Spline.fit.L1 <-   cbind(1, New.Days.L1 ) %*% coef(Fit.1L )
 Spline.fit.L1W <-  cbind(1, New.Days.L1 ) %*% coef(Fit.1L.W)
 
 ##1 knot, Cubic 
 (KN1 <- median(HSCT$Days)) #8
 bs.DAY.C1 <- bs(HSCT$Days, knots= KN1 , degree=3) 
 Fit.1C   <- lm(Granu.log ~ bs.DAY.C1)
 Fit.1C.W <- lm(Granu.log ~ bs.DAY.C1, weights=1/ni)
 
 New.Days.C1 <- bs(-7:35, knots=KN1, degree=3)
 Spline.fit.C1 <-   cbind(1, New.Days.C1 ) %*% coef(Fit.1C )
 Spline.fit.C1W <-  cbind(1, New.Days.C1 ) %*% coef(Fit.1C.W)
 
 ##2 knots, linear  
 ( KN2 <- quantile(Days, c(.33, .66)))
 bs.DAY.L2 <- bs(HSCT$Days, knots= KN2 , degree=1) 
 Fit.2L   <- lm(Granu.log ~ bs.DAY.L2)
 Fit.2L.W <- lm(Granu.log ~ bs.DAY.L2, weights=1/ni)
 
 New.Days.L2 <- bs(-7:35, knots=KN2, degree=1)
 Spline.fit.L2 <-   cbind(1, New.Days.L2 ) %*% coef(Fit.2L )
 Spline.fit.L2W <-  cbind(1, New.Days.L2 ) %*% coef(Fit.2L.W)
 

 postscript("fig4.1.granu.ps", horizontal=T)
  par(mar=c(4.5, 5.5, 3,1),cex=2,mfrow=c(2,2), cex.lab=1.6, cex.axis=1.4, cex.main=1.6)
 
  ## Subplot 1:   1 knot, linear 
  plot(Days , Granu.log,   main='Granulocyte (1 knot, linear)',  col='gray50', pch=1, cex=0.8,
      xlab='Days post-transplantation', ylab='',  xlim=c(-8,35) , ylim=c(-3,1.5), axes=F)
 
  axis(2,at= -3:1,  labels= c(1, 10,100,1000,10000), cex.axis=1.6 , las=2 )
  axis(1, at=seq(-7, 35, by=7),cex.axis=1.6); box() 
  lines(-7:35, Spline.fit.L1  ,col='gray40', lwd=2)
  lines(-7:35, Spline.fit.L1W ,col=1, lwd=2, lty=2)
   mtext("A.", cex=1.5, side=3, line=0.7, font=2, at=-15)
 
 ## Subplot 2:   1 knot, Cubic 
 plot(Days , Granu.log,   main='Granulocyte (1 knot, cubic)',  col='gray50', pch=1, cex=0.8,
      xlab='Days post-transplantation', ylab='',  xlim=c(-8, 35) , ylim=c(-3,1.5), axes=F)
   axis(2,at= -3:1,  labels= c(1, 10,100,1000,10000), cex.axis=1.6 , las=2 )
   axis(1, at=seq(-7, 35, by=7),cex.axis=1.6); box() 
   lines(-7:35, Spline.fit.C1  ,col='gray40', lwd=2)
   lines(-7:35, Spline.fit.C1W ,col=1, lwd=2, lty=2)
   mtext("B.", cex=1.5, side=3, line=0.7, font=2, at=-15)
 
  ## Subplot 3:   2 knot, linear 
  plot(Days , Granu.log,   main='Granulocyte (2 knots, linear)',  col='gray50', pch=1, cex=0.8,
      xlab='Days post-transplantation', ylab='',  xlim=c(-8, 35) , ylim=c(-3,1.5), axes=F)
  axis(2,at= -3:1,  labels= c(1, 10,100,1000,10000), cex.axis=1.6 , las=2 )
  axis(1, at=seq(-7, 35, by=7),cex.axis=1.6); box() 
  lines(-7:35, Spline.fit.L2  ,col='gray40', lwd=2)
  lines(-7:35, Spline.fit.L2W ,col=1, lwd=2, lty=2)
   mtext("C.", cex=1.5, side=3, line=0.7, font=2, at=-15)
 
 ## Subplot 4:   2 knot, cubic 
  plot(Days , Granu.log,   main='Granulocyte (2 knots, cubic)',  col='gray50', pch=1, cex=0.8,
      xlab='Days post-transplantation', ylab='',  xlim=c(-8, 35) , ylim=c(-3,1.5), axes=F)
   axis(2,at= -3:1,  labels= c(1, 10,100,1000,10000), cex.axis=1.6 , las=2 )
   axis(1, at=seq(-7, 35, by=7),cex.axis=1.6) ; box()
   lines(-7:35,  Spline.Est0 , col='gray40', lwd=2)
   lines(-7:35,  Spline.Est  , col=1,lwd=2, lty=2) # subject-uniform
    mtext("D.", cex=1.5, side=3, line=0.7, font=2, at=-15)
 dev.off()
 
 detach(HSCT)