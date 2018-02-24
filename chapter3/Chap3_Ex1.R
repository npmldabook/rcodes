### Chapter 3 Ex1 ###

## remove (almost) everything in the working environment.
rm(list = ls())
#----------------

library(npmlda)

str(HSCT)
HSCT[HSCT$ID==1,]

summary(HSCT$Granu)
summary(HSCT$LYM)
summary(HSCT$MON)

HSCT$Granu.log <- log10(HSCT$Granu)
HSCT$LYM.log   <- log10(HSCT$LYM)
HSCT$MON.log   <- log10(HSCT$MON)

# Code for Figure 3.1 analysis ##
Ct <- data.frame(table(HSCT$ID))
names(Ct)<- c("ID", "ni")
HSCT<- merge(HSCT, Ct, by= "ID")

Fit.Granu  <- with(HSCT[!is.na(HSCT$Granu.log),],
                  kernel.fit(sort(unique(Days)),Days, Granu.log, bw=4, Kernel="Ep", Wt=1/ni ))      
                      
Fit2.Granu  <- with(HSCT[!is.na(HSCT$Granu.log),],
                  kernel.fit(sort(unique(Days)),Days, Granu.log, bw=4, Kernel="Ep", Wt=1 ))      

#-- LYM fit  --
Fit.LYM  <- with(HSCT[!is.na(HSCT$LYM.log),],
                kernel.fit(sort(unique(Days)),Days, LYM.log, bw=3, Kernel="Nm", Wt=1/ni ))      

Fit2.LYM  <- with(HSCT[!is.na(HSCT$LYM.log),],
                 kernel.fit(sort(unique(Days)),Days, LYM.log, bw=3, Kernel="Nm", Wt=1 ))      
#-- MON fit --

Fit.MON  <- with(HSCT[!is.na(HSCT$MON.log),],
                kernel.fit(sort(unique(Days)),Days, MON.log, bw=5, Kernel="Bw", Wt=1/ni ))      

Fit2.MON  <- with(HSCT[!is.na(HSCT$MON.log),],
                kernel.fit(sort(unique(Days)),Days, MON.log, bw=5, Kernel="Bw", Wt=1 ))      


postscript("fig3.1.ps", horizontal=T)

par(mar=c(4.5, 5.5, 3,1),cex=2,mfrow=c(1,3), pty='s', cex.lab=1.6, cex.axis=1.4, cex.main=1.6)

#-----------PMN-----------
plot(HSCT$Days , HSCT$Granu.log,   main='Granulocyte (Epanechnikov)',  col='gray50', pch=1, cex=0.8,
     xlab='Days post-transplantation', ylab='',  xlim=c(-8,28) , ylim=c(-3,1.5), axes=F)

axis(2,at= -3:1,  labels= c(1, 10,100,1000,10000), cex.axis=1.6 , las=2 )
axis(1, at=seq(-7, 35, by=7),cex.axis=1.6) ;box()

GDays<- sort(unique(HSCT$Days[!is.na(HSCT$Granu.log)]))
lines(GDays , Fit2.Granu , col='gray40',lwd=2)
lines(GDays , Fit.Granu  , col=1, lwd=2, lty=2) 


#-----------LYM  -------------
plot(HSCT$Days , HSCT$LYM.log ,   main='Lymphocyte (Normal)',  col='gray50', pch=1, cex=0.8,
     xlab='Days post-transplantation', ylab='',  xlim=c(-8,28) , ylim=c(-3,0.5), axes=F)

axis(2,at= -3:0,  labels= c(1, 10,100,1000), cex.axis=1.6 , las=2 )
axis(1, at=seq(-7, 28, by=7),cex.axis=1.6) ;box()

LDays<- sort(unique(HSCT$Days[!is.na(HSCT$LYM.log)]))
lines(LDays , Fit2.LYM , col='gray40', lwd=2)
lines(LDays , Fit.LYM , col=1, lwd=2, lty=2)


#-----------MON -------------
plot(HSCT$Days , HSCT$MON.log,   main='Monocyte (Biweight)',  col='gray50', pch=1, cex=0.8,
     xlab='Days post-transplantation', ylab='',  xlim=c(-8,28) , ylim=c(-3,0.5), axes=F)

axis(2,at= -3:0,  labels= c(1, 10,100,1000), cex.axis=1.6 , las=2 )
axis(1, at=seq(-7, 28, by=7),cex.axis=1.6) ;box()

MDays <- sort(unique(HSCT$Days[!is.na(HSCT$MON.log)]))
lines(MDays , Fit2.MON, col='gray40', lwd=2)
lines(MDays , Fit.MON , col=1, lwd=2, lty=2)

dev.off()












