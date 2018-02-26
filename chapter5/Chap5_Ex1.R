### Chapter 5 Ex1 ###
# remove (almost) everything in the working environment.
rm(list = ls())
#----------------

library(npmlda)
#----------------

HSCT$LYM.log   <- log10(HSCT$LYM)
HSCT<- HSCT[!is.na(HSCT$LYM.log ),]

Ct <- data.frame(table(HSCT$ID))
names(Ct)<- c("ID", "ni")
HSCT<- merge(HSCT, Ct, by= "ID")

attach(HSCT)

plot(Days, LYM.log, xlab="Days post-transplantation",ylab="")

smfit<-smooth.spline(Days, LYM.log,spar=0.7,cv=NA)
smfit.w<-smooth.spline(Days, LYM.log, spar=0.7, cv=NA, w=1/ni)

lines(predict(smfit,-8:35), col ="gray40", lwd=1.5)
lines(predict(smfit.w, -8:35), lwd=1.5, lty=2)

## Code for Figure 5.1 analysis ##

LYM.spl.1 <- smooth.spline(Days, LYM.log, all.knots = T, spar = 0.2, cv=NA)
LYM.spl.1W <- smooth.spline(Days, LYM.log, all.knots = T, spar = 0.2, cv=NA, w=1/ni)

LYM.spl.2 <- smooth.spline(Days, LYM.log, all.knots = T, spar = 0.5, cv=NA)
LYM.spl.2W <- smooth.spline(Days, LYM.log, all.knots = T, spar = 0.5, cv=NA, w=1/ni)

LYM.spl.3 <- smooth.spline(Days, LYM.log, all.knots = T, spar = 0.7, cv=NA)
LYM.spl.3W <- smooth.spline(Days, LYM.log, all.knots = T, spar = 0.7, cv=NA, w=1/ni)

LYM.spl.4 <- smooth.spline(Days, LYM.log, all.knots = T, spar = 1.5, cv=NA)
LYM.spl.4W <-smooth.spline(Days, LYM.log, all.knots = T, spar = 1.5, cv=NA, w=1/ni)


postscript("fig5.1.LYM.ps", horizontal=T)

par(mar=c(4.5, 5.5, 3,1),cex=2,mfrow=c(2,2), cex.lab=1.4, cex.axis=1.4, cex.main=1.6)

plot(Days, LYM.log,   main=expression(paste("spar=0.2,  ", lambda, "= 2 x", 10^-6)),  col='gray50', pch=1, cex=0.8,
     xlab='Days post-transplantation', ylab='',  xlim=c(-8,35) , ylim=c(-3,0.5), axes=F)

axis(2, at= -3:0,  labels= c(1, 10,100,1000), cex.axis=1.6 , las=2 ) 
axis(1, at=seq(-7, 35, by=7), cex.axis=1.6, las=1 ,mgp=c(3,1,0)) 
box()

mtext("A.", cex=1.5, side=3, line=0.7, font=2, at=-15)

lines(predict(LYM.spl.1, (-8:35)), col ='gray40', lwd=1.5)
lines(predict(LYM.spl.1W, (-8:35)), lwd=1.5, lty=2)

#----------------------#
plot(Days, LYM.log,   main=expression(paste("spar=0.5,  ", lambda, "= 3.1 x", 10^-4)),  col='gray50', pch=1, cex=0.8,
     xlab='Days post-transplantation', ylab='',  xlim=c(-8,35) , ylim=c(-3,0.5), axes=F)

axis(2, at= -3:0,  labels= c(1, 10,100,1000), cex.axis=1.6 , las=2 ) 
axis(1, at=seq(-7, 35, by=7), cex.axis=1.6, las=1,mgp=c(3,1,0)) 
box()

mtext("B.", cex=1.5, side=3, line=0.7, font=2, at=-15)

lines(predict(LYM.spl.2, (-8:35)), col ='gray40', lwd=1.5)
lines(predict(LYM.spl.2W, (-8:35)), lwd=1.5, lty=2)

#----------------------#
plot(Days, LYM.log,   main=expression(paste("spar=0.7,  ", lambda, "= 8.7 x", 10^-3)),  col='gray50', pch=1, cex=0.8,
     xlab='Days post-transplantation', ylab='',  xlim=c(-8,35) , ylim=c(-3,0.5), axes=F)

axis(2, at= -3:0,  labels= c(1, 10,100,1000), cex.axis=1.6 , las=2 ) 
axis(1, at=seq(-7, 35, by=7), cex.axis=1.6, las=1,mgp=c(3,1,0)) 
box()

mtext("C.", cex=1.5, side=3, line=0.7, font=2, at=-15)

lines(predict(LYM.spl.3, (-8:35)), col ='gray40', lwd=1.5)
lines(predict(LYM.spl.3W, (-8:35)), lwd=1.5, lty=2)


#----------------------#
plot(Days, LYM.log,   main=expression(paste("spar=1.5,  ", lambda, "= 5216")),  col='gray50', pch=1, cex=0.8,
     xlab='Days post-transplantation', ylab='',  xlim=c(-8,35) , ylim=c(-3,0.5), axes=F)

axis(2, at= -3:0,  labels= c(1, 10,100,1000), cex.axis=1.6 , las=2 ) 
axis(1, at=seq(-7, 35, by=7), cex.axis=1.6, las=1,mgp=c(3,1,0)) 
box()

mtext("D.", cex=1.5, side=3, line=0.7, font=2, at=-15)

lines(predict(LYM.spl.4, (-8:35)), col ='gray40', lwd=1.5)
lines(predict(LYM.spl.4W, (-8:35)), lwd=1.5, lty=2)

dev.off()

## Code for Figure 5.2 analysis ##
# leave one-subject out Cross-validation Function # 
CV.smfit <- function( Xvec, Yvec,  ID, Spar, Wt )
{
  NN <-length( Yvec)
  Yest <- numeric(NN)
  nID <- length(unique(ID))
  for (i in 1:nID)
  {
    # print(i)
    Xsub <- Xvec[ ID != i ]
    Ysub <- Yvec[ ID != i ]
    X.IDi <- Xvec[ ID == i ]
    Wtsub<- Wt[ ID != i ]
    Smfit <-  smooth.spline( Xsub , Ysub , all.knots = T, spar = Spar , cv=NA , w=Wtsub) 
    Yest[ID == i] <-  predict(Smfit , X.IDi )$y
  }
  sum(Wt*(Yvec- Yest)^2 )  
}
#--------------------------------------------------#


Par.range <- seq(0.1, 0.9, by =0.001)
NK <-  length(Par.range)

CVh <- CVh.W <- numeric(NK )

for (j in 1:NK )
{
  #print(j)
  CVh[j]  <-    CV.smfit( Days, LYM.log, ID,  Spar =Par.range[j] ,Wt=rep(1/length(LYM.log),length(LYM.log) ))
  CVh.W[j] <-   CV.smfit( Days, LYM.log, ID,  Spar =Par.range[j], Wt=1/(ni*length(unique(ID))) )
}

Par.range[which.min(CVh)]#0.555
Par.range[which.min(CVh.W)] # 0.568

(LYM.spl.CV   <- smooth.spline( Days, LYM.log, all.knots = T, spar =0.555, cv=NA))
(LYM.spl.CV.W <- smooth.spline( Days, LYM.log, all.knots = T, spar =0.568, cv=NA, w= 1/ni))

#---------
postscript("fig5.2.LYM-CV.ps", horizontal=T)

par(mar=c(20.5, 4.5, 3, 1),cex=2, mfrow=c(1, 2), pty='m', cex.lab=1.4, cex.axis=1.2,  cex.main=1.2)
## plot 1: CV scores

plot( Par.range, CVh, type='l', pch=16,   xlab='spar', ylab='CV score', lwd=1.5, axes=T, ylim=c(0.20, 0.40))
lines( Par.range ,CVh.W , col=1, type='l',lty=2, lwd=1.5)
box()
abline(v=0.56, col='gray20')
abline(v=0.57, col='gray20', lty=2)

legend('topleft', lty=c(1,2), legend=c('measurement uniform weight','subject uniform weight'), bty='n',lwd=1.5, cex=1.1)
mtext("A.", cex=1.5, side=3, line=0.7, font=2, at=-0.05)


## plot : LYM with CV fitted curve
plot(Days, LYM.log,   main="" ,  col='gray50', pch=1, cex=0.8,
     xlab='Days post-transplantation', ylab='Lymphocyte count',  xlim=c(-8,35) , ylim=c(-3.4,0.5), axes=F, )

axis(2, at= -3:0,  labels= c(1, 10,100,1000), las=2, mgp=c(3,0.6,0) ) 
axis(1, at=seq(-7, 35, by=7),  las=1,mgp=c(3,1,0)) 
box()

mtext("B.", cex=1.5, side=3, line=0.7, font=2, at=-17)

lines(predict(LYM.spl.CV, (-8:35)), col ='gray40', lwd=1.5)
lines(predict(LYM.spl.CV.W, (-8:35)), lwd=1.5, lty=2)

legend('bottomright', lty=c(1,2), legend=c(expression(paste("spar=0.56, ", lambda, "=0.00078")),
                                          expression(paste("spar=0.57, ", lambda, "=0.00096"))), bty='n',
       lwd=1.5, cex=1.1)

dev.off()

#-------------
detach(HSCT)





