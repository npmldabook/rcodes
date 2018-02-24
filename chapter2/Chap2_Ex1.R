### Chapter 2###

library (npmlda)  
library(nlme)
library(lme4)
#-----------------

str(BMACS)
head(BMACS)

# random intercept
CD4fit1 <- lme(CD4~ Time , random=~1 |ID, data= BMACS)
summary(CD4fit1 )

# random intercept and slop
CD4fit2 <- lme(CD4~ Time , random=~Time |ID, data= BMACS)
summary(CD4fit2 )

# compare models
anova(CD4fit1, CD4fit2)


# intercept and slope

CD4fit2b <- lmer(CD4~ Time + (Time|ID), data= BMACS)
summary(CD4fit2b )

# add covariates
BMACS$preCD4c <- BMACS$preCD4 - mean(BMACS$preCD4)
CD4fit3<- lme(CD4~ Time + preCD4c + Smoke + age,
                 random=~Time|ID, data=BMACS)
summary(CD4fit3)

#############################

CD.pt1<-BMACS[BMACS$ID==  2114,]
CD.pt2<-BMACS[BMACS$ID==  9044,]

postscript("fig2.1.ps", horizontal=T)

par(mar=c(4.5, 4.5, 3, 1),cex=2,mfrow=c(2,3), cex.lab=1.6, cex.axis=1.3, cex.main=1.4)

## mixed effect model (2.48) : random intercept

# first line of plot
plot( CD4~ Time, data=BMACS, pch=1, ylim=c(0 , 63), cex=0.7,col='gray30',  xlab="", ylab="CD4 percentage (%)", main="Model 1: All subjects")
abline(CD4fit1$coef$fixed, lty=1, col=1, lwd=2)


plot( CD4~ Time, data=CD.pt1, pch=1, ylim=c(0 , 63),  xlab="", ylab="",  main="Model(2.48): Subject 1")
abline(CD4fit1$coef$fixed, lty=1,  lwd=1.5,col=1)
abline(coef(CD4fit1 )[29,1] , coef(CD4fit1 )[29,2]   , lwd=1.5,lty=2, col=1)

plot( CD4~ Time, data=CD.pt2, pch=1, ylim=c(0 , 63), xlab="",  ylab="",  main="Model(2.48): Subject 2")
abline(CD4fit1$coef$fixed, lty=1,  lwd=1.5,col=1)
abline(coef(CD4fit1 )[259,1] , coef(CD4fit1 )[259,2]  , lty=2, lwd=1.5, col=1)


## mixed effect model (2.49): random  intercept and slope

plot( CD4~ Time, data=BMACS, pch=1, ylim=c(0 , 63), cex=0.7, col='gray30',  ylab="CD4 percentage (%)",  main="Model 2: All subjects")
abline(CD4fit2$coef$fixed, lty=1, col=1, lwd=2)

plot( CD4~ Time, data=CD.pt1, pch=1, ylim=c(0 , 63), ylab="",  main="Model(2.49): Subject 1")
abline(CD4fit2$coef$fixed, lty=1, lwd=1.5, col=1)
abline(coef(CD4fit2 )[29,1] , coef(CD4fit2 )[29,2]   , lty=2,  lwd=1.5,col=1 )


plot( CD4~ Time, data=CD.pt2, pch=1, ylim=c(0 , 63), ylab="", main="Model(2.49): Subject 2")
abline(CD4fit2$coef$fixed, lty=1, lwd=1.5, col=1)
abline(coef(CD4fit2 )[259,1] , coef(CD4fit2 )[259,2]   , lty=2,  lwd=1.5, col=1)

dev.off()





