### Chapter 4, Ex2 ###

## remove (almost) everything in the working environment.
rm(list = ls())

#library
library (npmlda) 
library (splines) 


### Fig 4.2 analysis ##
Ct <- data.frame(table(BMACS$ID))
Ct <- data.frame( Ct, 1:nrow(Ct) )
names(Ct)<- c("ID", "ni", "IDD")
BMACS<- merge(BMACS, Ct, by= "ID")
str(BMACS)

attach(BMACS)
newX <- seq(min(Time), max(Time), by=0.1)

fit5  <- spline.fit(newX, Time, CD4, nKnots=5, Degree=3)
fit5W <- spline.fit(newX, Time, CD4, nKnots=5, Degree=3, Wt=1/ni)

fit1  <- spline.fit(newX, Time, CD4, nKnots=1, Degree=3)
fit1W <- spline.fit(newX, Time, CD4, nKnots=1, Degree=3, Wt=1/ni)

#--------------
postscript("fig4.2.ps", horizontal=T)

par(mar=c(4.5, 4.5, 3, 1),cex=2, mfrow=c(2, 2), pty='m', cex.lab=1.2,cex.axis=1.1,  cex.main=1.2)

plot(CD4 ~ Time, data = BMACS,  xlab = "Time since infection (years)", ylab = "CD4 percentage",  ylim=c(0,65), cex=0.7, col='gray50', main="Spline fit: 1 knot") 
lines(newX, fit1, col=1,  lwd=1.5)
lines(newX, fit1W, col=1, lwd=1.5, lty=2)
mtext("A.", cex=1.5, side=3, line=0.7, font=2, at=-0.9)

#
plot(CD4 ~ Time, data = BMACS,  xlab = "Time since infection (years)", ylab = "CD4 percentage",  ylim=c(0,65), cex=0.7, col='gray50', main="Spline fit: 5 knots") 
lines(newX, fit5  , col=1,  lwd=1.5)
lines(newX, fit5W  , col=1, lwd=1.5, lty=2)
mtext("B.", cex=1.5, side=3, line=0.7, font=2, at=-0.9)

dev.off()
#--------------

fit10 <- spline.fit(newX, Time, CD4, nKnots=10, Degree=3)
fit15 <- spline.fit(newX, Time, CD4, nKnots=15, Degree=3)
fit20 <- spline.fit(newX, Time, CD4, nKnots=20, Degree=3)

# Cross-validation ##

NK <- 20  
CVh <- numeric(NK )

for (j in 1:NK )
{
  print(j)
  CVh[j]  <-   CVspline( Time,  CD4 , IDD,  nKnots=j, Degree=3, Wt=rep(1/nrow(BMACS), nrow(BMACS)))
}

(1:NK)[which.min(CVh)] #5 knots is selected for uniform weight fit 

postscript("fig4.3_CD4.ps", horizontal=T)

par(mar=c(4.5, 4.5, 3, 1),cex=2, mfrow=c(2, 2), pty='m', cex.lab=1.2,cex.axis=1.2,  cex.main=1.2)
 plot(CD4 ~ Time, data = BMACS,  xlab = "Time since infection (years)", ylab = "CD4 percentage",  ylim=c(0,65), cex=1, col='gray50', main="", type='n') 
 lines(newX, fit1,  lty=1, lwd=1, col='gray20')
 lines(newX, fit5, col=1, lty=2, lwd=1.5)
 lines(newX, fit10,col=1, lty=3, lwd=1.5)
 lines(newX, fit15,col=1, lty=4, lwd=1.5)
 lines(newX, fit20,col=1, lty=5, lwd=1.5)
 legend('topright', lty=1:5, col=c('gray30','black','black','black','black'), legend=c("nk=1","nk=5","nk=10", "nk=15" ,"nk=20"), bty='n', cex=1.1)
#---
plot(1:NK,CVh, pch=16, col='gray40', xlab='Number of knots (nk)', ylab="cross-validation score", axes=F)
lines(1:NK,CVh, lty=2)
axis(2)
axis(1, at=c(1,5,10,15,20))
box()

dev.off()

## Bootstrap CI: method1: use Chap3.5 or use loop below ##

Time.int<- seq(0.1,5.9,  by=0.1)
NN<- length(Time.int) ; NN

nBoot<- 1000
BootCD4FIT <- matrix(NA, nrow= nBoot, ncol= NN)

set.seed(101)
nID <- length(unique(BMACS$ID))

for (j in 1:nBoot)
{
  if ((j-floor(j/50)*50)== 0 ) print(j) 
  Index.ID <-  sample(nID, nID, replace=T) 
  Bootdata <- NULL 
  
  for ( i in 1:nID)
  {  
    new <- BMACS[BMACS$IDD== Index.ID[i],] 
    Bootdata<- rbind( Bootdata, new)
  } 
  BootCD4FIT[j,] <-  spline.fit(Time.int, Bootdata$Time, Bootdata$CD4, nKnots=5, Degree=3)    
}

# by estimate +/- SD
CD4.SD <-  sqrt(apply( BootCD4FIT,  2, var )) 

EstUpp.CI <-   fit5   + qnorm(.975) * CD4.SD 
EstLow.CI <-   fit5   + qnorm(.025) * CD4.SD 


EstUpp.Simul.CI <-  fit5   + qnorm(1-.025/60) * CD4.SD 
EstLow.Simul.CI <-  fit5   + qnorm(.025/60)* CD4.SD 

# percentile CI 
EstUpp.pCI <-    apply( BootCD4FIT,  2, quantile, prob=.025 )
EstLow.pCI <-   apply( BootCD4FIT,  2, quantile, prob=.975 )

EstUpp.Simul.pCI <-  apply( BootCD4FIT,  2, quantile, prob=.025/60 )
EstLow.Simul.pCI <-  apply( BootCD4FIT,  2, quantile, prob=(1-.025/60) )


postscript("fig4.4_CI.ps", horizontal=T)

par(mar=c(20.5, 4.5, 3, 1),cex=2, mfrow=c(1, 2), pty='m', cex.lab=1.4, cex.axis=1.2,  cex.main=1.2)

plot(CD4 ~ Time, data = BMACS,  xlab = "Time since infection (years)", ylab = "CD4 percentage",   cex=0.3, col='gray70', main="") 
 polygon(c(Time.int[1], Time.int, rev(Time.int)),  c(EstLow.CI[1], EstUpp.CI, rev(EstLow.CI)),  col="gray60", border=NA)
 points(CD4 ~ Time, data = BMACS,  xlab = "", ylab = "",   cex=0.3, col='gray10',pch=16,  main="") 
 lines(Time.int,fit5, lwd=2.5, col=1)
 mtext("A.", cex=1.5, side=3, line=0.7, font=2, at=-0.9)
 grid()

#---
plot(CD4 ~ Time, data = BMACS,  xlab = "Time since infection (years)", ylab = "CD4 percentage",   cex=0.3, col='gray70', main="") 
 polygon(c(Time.int[1], Time.int, rev(Time.int)),  c(EstLow.Simul.CI[1], EstUpp.Simul.CI, rev(EstLow.Simul.CI)),  col="gray60", border=NA)
 points(CD4 ~ Time, data = BMACS,  xlab = "", ylab = "",   cex=0.3, col='gray10',pch=16,  main="")  
 lines(Time.int,fit5, lwd=2.5, col=1)
 mtext("B.", cex=1.5, side=3, line=0.7, font=2, at=-0.9)
 grid()

dev.off()































