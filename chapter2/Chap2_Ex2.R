rm(list = ls())
#----------------

library(npmlda)
search()
ls()

#data(BDIdata) #optional
str(BDIdata)
BDIdata[BDIdata$ID==1,]

library(nlme)
# recode time in months
 
 BDIdata$Tijm <- BDIdata$time*12/365.25
 BDIsub <- subset(BDIdata, med.time >=0 & med.time < 200)
 dim(BDIsub)# 1465    6

# Model (2.50)
 BDI.Model3 <- lme(BDI ~ Tijm, data=BDIsub, random=~Tijm|ID)
 summary(BDI.Model3)
 
 # Model (2.51)
 BDIsub$Sim <- BDIsub$med.time*12/365.25
 BDIsub$Rijm <- with(BDIsub, med*(Tijm -Sim))
 BDI.model4 <- lme(BDI ~ Tijm+ med + Rijm , data=BDIsub,
                     random=~Tijm+ med + Rijm|ID)
 summary(BDI.model4)
#--------------------------
 
## Generate figure 2.2 ##
# patient-1 has 18 observations
 
 ID1<-3 
 table(BDIsub$ID) 

 BDI.pt1<-  BDIsub[BDIsub$ID==ID1,]


 new.obs<- cbind(1, 3.4,  0, 0) 
 xx<-  as.matrix(cbind(1, BDI.pt1[, c('Tijm', 'med', 'Rijm')]))
 xx<- rbind(xx, new.obs)
 xx<- xx[order(xx[,2]),] 
 
 #18 *4 
 fit.y1 <- (xx[1:11,]) %*% (BDI.model4$coef$fixed)
 fit.y2 <- (xx[12:19,])%*% (BDI.model4$coef$fixed)
 
 # 
 ranfit.y1 <- (xx[1:11,]) %*% t(coef(BDI.model4)[ID1,] )
 ranfit.y2 <- (xx[12:19,])%*% t(coef(BDI.model4)[ID1,] )
 
 
 #Pt-2 is pt #11 
 
 ID2<- 4
 BDI.pt2<-  BDIsub[BDIsub$ID==ID2,]
 
 new.obs<- rbind(c(1,  2,  0, 0), c(1,  2.01,  1,  2.01-2.00411) ) 
 xx2<-  as.matrix(cbind(1, BDI.pt2[, c('Tijm', 'med', 'Rijm')]))
 xx2<- rbind(xx2, new.obs)
 xx2<- xx2[order(xx2[,2]),] 
 
 #18 *4 
 P2.fit.y1 <- (xx2[1:10,]) %*% (BDI.model4$coef$fixed)
 P2.fit.y2 <- (xx2[11:27,])%*% (BDI.model4$coef$fixed)
 
 # 
 P2.ranfit.y1 <- (xx2[1:10,]) %*% t(coef(BDI.model4)[ID2, ] )
 P2.ranfit.y2 <- (xx2[11:27,])%*% t(coef(BDI.model4)[ID2, ] )
 
 
#----------------------------------------------------
 
 #CD4fit1$coef$fixed == 'fixed effects' #
 
 ##  note the following code is not most efficient. can later write a shorter code ##
 postscript("fig2.2.ps", horizontal=T)
 
 par(mar=c(4.5, 4.5, 3, 1),cex=2,mfrow=c(2,2), cex.lab=1.6, cex.axis=1.3, cex.main=1.3)
 
 ## random model (2.50) ; Pt1
 
 plot( BDI~ Tijm, data=BDI.pt1, pch=1, cex=0.8, xlab="",ylim=c(5 , 30),  ylab="BDI score",
       main="Model (2.50): Subject 1")
 abline(BDI.Model3$coef$fixed, lty=1, col=1, lwd=1)
 abline(v= BDI.pt1$Sim[1], lty=2, lwd=1.5, col=1)
 abline(coef(BDI.Model3)[ID1,1] , coef(BDI.Model3)[ID1,2]  , lty=2, lwd=1.5, col=4)
 
 
 ## random model 3 ; Pt2
 
 plot( BDI~ Tijm, data=BDI.pt2, pch=1, cex=0.8, xlab="",ylim=c(5 , 30),  ylab="", main="Model (2.50): Subject 2")
 abline(BDI.Model3$coef$fixed, lty=1, col=1, lwd=1)
 abline(v= BDI.pt2$Med.Time[1], lty=2, lwd=1.5, col=1)
 abline(coef(BDI.Model3)[ID2,1] , coef(BDI.Model3)[ID2,2]  ,lwd=1.5,  lty=2, col=4)
 
 
 ##  random model (2.51) ; Pt1
 
 plot( BDI~ Tijm, data=BDI.pt1, pch=1, cex=0.8, ylim=c(5, 30),  xlab="Time (months)", 
        ylab="BDI score", main="Model (2.51): Subject 1")
 abline(v= BDI.pt1$Sim[1], lty=2, col=1)
 lines(xx[1:11,2], fit.y1, lty=1, lwd=1.5, col=1)
 lines(xx[12:19,2], fit.y2, lty=1, lwd=1.5, col=1)
 
 lines(xx[1:11,2], ranfit.y1, lty=2, lwd=1.5, col=4)
 lines(xx[12:19,2],ranfit.y2, lty=2,lwd=1.5,  col=4)
 
 
 ## random model (2.51): pt2
 
 plot( BDI~ Tijm, data=BDI.pt2, pch=1, cex=0.8, ylim=c(5 , 30),
       xlab="Time (months)",  ylab="BDI score", main="Model (2.51): Subject 2")
 abline(v= BDI.pt2$Sim[1], lty=2, col=1)
 lines(xx2[1:10,2], P2.fit.y1, lty=1,lwd=1.5,  col=1)
 lines(xx2[11:27,2], P2.fit.y2, lty=1,lwd=1.5,  col=1)
 
 lines(xx2[1:10,2], P2.ranfit.y1, lty=2,lwd=1.5,  col=4)
 lines(xx2[11:27,2],P2.ranfit.y2, lty=2, lwd=1.5, col=4)
 
 dev.off()
 
 
 
 
 
 
 
 
 
 