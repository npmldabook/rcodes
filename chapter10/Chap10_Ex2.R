## Code for Ch 10, Sec 10.7.1##

rm(list = ls())
#----------------

library(npmlda)
library(nlme)
library(MASS)
------------------

#data(BDIdata)  #optional
str(BDIdata)

BDI.firstob <- do.call("rbind", as.list(by(BDIdata, BDIdata$ID, head,  n=1)))

dim(BDI.firstob)#557

table(BDI.firstob$med.time)  
#       11,   45,              47 or      454 patients used med
#   before, at baseline, or during 6 m, or after 6m

# recode time in months
 BDIdata$Tijm <- BDIdata$time*12/365.25
 BDIsub <- subset(BDIdata, med.time >=0 & med.time < 200)
 dim(BDIsub)# 1465    6

# LME Model (2.51), output in Table 10.3
 BDIsub$Sim <- BDIsub$med.time*12/365.25
 BDIsub$Rijm <- with(BDIsub, med*(Tijm -Sim))
 BDI.model4 <- lme(BDI ~ Tijm+ med + Rijm , data=BDIsub,
                     random=~Tijm+ med + Rijm|ID)
 summary(BDI.model4)

 ##  Firgure 10.1  ##
 
 BDI.pt1 <-  BDIsub[ BDIsub$ID== 19,]
 BDI.pt2 <-  BDIsub[ BDIsub$ID==20,]
 BDI.pt3 <-  BDIsub[ BDIsub$ID==4,]
 BDI.pt4 <-  BDIsub[ BDIsub$ID==16,]
 

 postscript("fig10.1ABC.ps", horizontal=T)
 
 par(mar=c(4.5, 4.5, 3, 1.5),cex=2,mfrow=c(2,2), cex.lab=1.6, cex.axis=1.5, cex.main=1.5, font.main=1)
 
 plot( BDI~ time, data=BDI.pt1, pch=1, cex=0.9, xlab="Days",ylim=c(0 ,57),  ylab="BDI score", main="Subject 1: med started at day 13")
 lines(BDI~ time, lty=2 , data=BDI.pt1)
 abline(v= BDI.pt1$med.time[1], lty=1, lwd=1.5, col='gray30')
 mtext("A.", cex=1.3, side=3, line=0.9, font=2, at= -30)
 
 plot( BDI~ time, data=BDI.pt2, pch=1, cex=0.9, xlab="Days",ylim=c(0 ,57),  ylab="BDI score", main="Subject 2: med startedat day 10")
 lines(BDI~ time, lty=2 , data=BDI.pt2)
 abline(v= BDI.pt2$med.time[1], lty=1, lwd=1.5, col='gray30')
 mtext("B.", cex=1.3, side=3, line=0.9, font=2, at= -30)
 
 plot( BDI~ time, data=BDI.pt3, pch=1, cex=0.9, xlab="Days",ylim=c(0 ,57),  ylab="BDI score", main="Subject 3: med started at day 61")
 lines(BDI~ time, lty=2 , data=BDI.pt3)
 abline(v= BDI.pt3$med.time[1], lty=1, lwd=1.5, col='gray30')
 mtext("C.", cex=1.3, side=3, line=0.9, font=2, at= -30)
 
 plot( BDI~ time, data=BDI.pt4, pch=1, cex=0.9, xlab="Days",ylim=c(0 ,57),  ylab="BDI score", main="Subject 4: med started at day 127")
 lines(BDI~ time, lty=2 , data=BDI.pt4)
 abline(v= BDI.pt4$med.time[1], lty=1, lwd=1.5, col='gray30')
 mtext("D.", cex=1.3, side=3, line=0.9, font=2, at= -30)
 
 dev.off()
 
 #-------- VCME ---- #

 BDIsub$Tijm <- BDIsub$time*12/365.25
 BDIsub$Sim <- BDIsub$med.time*12/365.25
 BDIsub$TijSim <- with(BDIsub,Tijm*Sim)
 BDIsub$Med <- with(BDIsub, ifelse(time-med.time>=0, 1,0))
 BDIsub$Rijm <- with(BDIsub, Med*(Tijm -Sim))
 # Fit the VCME model
 VCME.fit <- lme(BDI ~ 1+Sim + Tijm + TijSim + Med + Rijm,
                   random=~Tijm + Med + Rijm|ID, data=BDIsub)
 summary(VCME.fit)

 # Fixed effects: BDI ~ 1 + Sim + Tijm + TijSim + Med + Rijm 
 #                  Value Std.Error   DF   t-value p-value
 # (Intercept) 25.640517 1.4418923 1369 17.782546  0.0000
 # Sim         -1.377896 0.5901003   90 -2.335021  0.0218
 # Tijm        -0.237399 0.8153425 1369 -0.291165  0.7710
 # TijSim       0.069672 0.1728152 1369  0.403161  0.6869
 # Med         -4.489470 1.0506389 1369 -4.273086  0.0000
 # Rijm        -2.052332 0.7675756 1369 -2.673785  0.0076
 
 
 ## Code for Ch 10, Sec 10.7.2##
 
  #Sec 10.7.2 : Naive Model  #
 dim(BDIdata ) # No. of observations
 length(unique(BDIdata$ID)) # No. of patients
 
 # Recode the covariates and time variables in months
 BDIdata$Tijm <- BDIdata$time*12/365.25
 BDIdata$med.time[BDIdata$med.time <0] <- 0
 BDIdata$Med <- with(BDIdata, ifelse(time-med.time>=0, 1,0))
 BDIdata$Rijm <- with(BDIdata, Med*(time-med.time)*12/365.25)
 Naive.LME <- lme(BDI ~ Tijm + Med + Rijm, data=BDIdata,
                  random=~ Tijm + Med + Rijm|ID )
 summary(Naive.LME)
 
 # Fixed effects: BDI ~ Tijm + Med + Rijm 
 #                 Value Std.Error   DF   t-value p-value
 # (Intercept) 14.453470 0.3116090 6557  46.38336  0.0000
 # Tijm        -1.887150 0.0669985 6557 -28.16707  0.0000
 # Med          3.579455 0.8252207 6557   4.33757  0.0000
 # Rijm         0.035481 0.2269312 6557   0.15635  0.8758
 
 
 
 
 