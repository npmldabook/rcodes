### Chapter 11 Ex1 ###
rm(list = ls())
#----------------
library(npmlda)
library(splines)
library(lme4)

## starting here for the chapter 11 #3

#data(BMACS)
T.range<- range(BMACS$Time)
Nk <- 4
KN <- seq(from=T.range[1], to=T.range[2], length=Nk)[-c(1,Nk)]
bs.time <- bs(BMACS$Time, knots=KN, degree=2, intercept=T)
fmCD4 <- lmer(CD4 ~ 0+bs.time+(0+bs.time|ID), data=BMACS)
summary(fmCD4)

####
N <- 50
Tgrid <- seq(from=T.range[1], to=T.range[2], length=N)
BS <- bs(Tgrid, knots=KN, degree=2, intercept=T)

# Obtain the estimate mean curve #
mean.hat <- BS %*% fixef(fmCD4)

# obtain the estimate of covariance for random effects and
# outcome

GAMMA0 <- matrix(NA, 5, 5)
VC <- as.data.frame(VarCorr(fmCD4))
diag(GAMMA0) <- VC$vcov[1:5]
GAMMA0[lower.tri(GAMMA0)] <- VC$vcov[6:15]
GAMMA <- t(GAMMA0)
GAMMA[lower.tri(GAMMA)] <- GAMMA0[lower.tri(GAMMA0)]
COV.MAT <- matrix(NA, N, N)
Index <- expand.grid(1:N, 1:N)
f <- function(i){sum( GAMMA * (BS[Index[i,1],] %o% BS[Index[i,2],]))}

COV.MAT <- matrix(do.call('rbind', lapply(1:(N^2),f)), nrow=N, byrow=F) #50*50

# obtain the eigenvalue and eigenvectors
Eigenfun <- eigen(COV.MAT)

(prop<-Eigenfun$values/sum(Eigenfun$values))
round(prop*100, 1)
#[1] 85.3 11.3  2.0  1.4#


Nsub <- 283
IDlevel<- rownames(coef(fmCD4)[[1]])
BLUP.est <- matrix(NA, nrow= Nsub, ncol= N)
Proj1 <- Proj2 <- numeric(Nsub)

for (i in 1:Nsub)
{
  Datai <- BMACS[BMACS$ID==IDlevel[i],]
  BLUP.est[i,] <- BS %*% t(as.vector(coef(fmCD4)[[1]][i,]))
  Proj1[i] <- BLUP.est[i,] %*% Eigenfun$vectors[,1]
  Proj2[i] <- BLUP.est[i,] %*% Eigenfun$vectors[,2]
}

#---------figure 11.1 ---------------------#

postscript("fig11-1.ps", horizontal=T)

par(mar=c(4.5, 4.5, 3, 1),cex=1, mfrow=c(2, 2), pty='m',  cex.axis=1.2, cex.main=1.3, lwd=1)

## overall fit
plot(BMACS$CD4 , BMACS$Time,  
     xlab = "Time since infection (years)", ylab = "CD4 percentage", type = "n", xlim=c(0,6), 
     ylim=c(0,66), cex.lab=1.1)
for (i in unique(BMACS$ID)) { lines(CD4 ~ Time, data = BMACS[BMACS$ID  == i, ], lty = 1, col='gray40')}
lines( Tgrid, mean.hat, col=1, lwd=2)

mtext("A.", cex=1.5, side=3, line=0.7, font=2, at=-0.9)

#---covariance structure

par(mar=c(0.5, 2.5, 1.5, 0),cex=0.8, pty='m', cex.axis=1)

mtext("B.", cex=1.5, side=3, line=-1, font=2)

Z.surf <-COV.MAT

### interpolate a set of given colors to create new color palette ###
col.pal <- colorRampPalette(c("gray90", "gray40"))
colors <- col.pal(100)

z.facet.center <- (Z.surf[-1, -1] + Z.surf[-1, -ncol(Z.surf)] + Z.surf[-nrow(Z.surf), -1]
                   + Z.surf[-nrow(Z.surf), -ncol(Z.surf)])/4

z.facet.range <- cut(z.facet.center, 100)

Surf <- persp(Tgrid, Tgrid, COV.MAT, theta=-60, phi=20, r=2, shade=0.5, axes=TRUE,scale=TRUE, box=TRUE, 
              nticks=5, ticktype="simple", zlim = c(0, 210) , col=colors[z.facet.range], xlab="Time", 
              ylab="Time", border=NA, zlab="Covariance", main="" , mgp= c(3, 1.2, 0),
              expand = 0.7, cex.lab=1.1, cex.axis=1.2)

lines (trans3d(x=Tgrid, y=Tgrid, z=diag(Z.surf), pm=Surf),lwd=2, col='gray80')


#-- eigen 1 ---#
par(mar=c(4.5, 4.5, 3, 1),cex=1, pty='m', cex.lab=1.1, cex.axis=1, cex.main=1.3)

plot(Tgrid, Eigenfun$vectors[,1], main='', xlab='Time', ylab="", type='l')
mtext("C.", cex=1.5, side=3, line=0.7, font=2, at=-0.9)

##--eigen 2 ---#

plot(Tgrid, Eigenfun$vectors[,2], main='', xlab='Time', ylab="",  type='l')
mtext("D.", cex=1.5, side=3, line=0.7, font=2, at=-0.9)

dev.off()

#---------figure 11.2 ---------------------#

cbind(No=1:Nsub, Proj1,Proj2)[Proj1== min(Proj1),]   #no. 171 
cbind(No=1:Nsub, Proj1,Proj2)[Proj1 == max(Proj1),]  #no. 78
cbind(No=1:Nsub, Proj1,Proj2)[Proj2== min(Proj2),]   #no.  14
cbind(No=1:Nsub, Proj1,Proj2)[Proj2 == max(Proj2),]  #no.  59


postscript("fig11-2.ps", horizontal=T)

par(mar=c(4.5, 4.5, 3, 1),cex=1, mfrow=c(2, 2), pty='m', cex.lab=1.2, cex.axis=1.2, cex.main=1.3, lwd=1)


## check using CD4 data or use HSCT data ### 
IDD<-  171   
Subject1   <- BMACS[BMACS$ID==IDlevel[IDD],]

plot(Subject1$Time, Subject1$CD4, xlim=c(0,6), ylim=c(0, 66), main= "" , 
     xlab = "Time since infection (years)", ylab = "CD4 percentage")

for (i in unique(BMACS$ID)) { lines(CD4 ~ Time, data = BMACS[BMACS$ID  == i, ], lty = 2, col='gray60')}

points(Subject1$Time, Subject1$CD4, xlim=c(0,6), ylim=c(0, 70),lwd=1.5 )

lines( Tgrid, mean.hat, col=1, lwd=1.5 )
lines( Tgrid, BLUP.est[IDD,], col=1, lty=2, lwd=1.5 )
mtext("A.", cex=1.5, side=3, line=0.7, font=2, at=-0.9)

## plot 1 ##
IDD<-  78 #171   
Subject1   <- BMACS[BMACS$ID==IDlevel[IDD],]

plot(Subject1$Time, Subject1$CD4, xlim=c(0,6), ylim=c(0, 66), main= "" , 
     xlab = "Time since infection (years)", ylab = "CD4 percentage")

for (i in unique(BMACS$ID)) { lines(CD4 ~ Time, data = BMACS[BMACS$ID  == i, ], lty = 2, col='gray60')}

points(Subject1$Time, Subject1$CD4, xlim=c(0,6), ylim=c(0, 70),lwd=1.5 )

lines( Tgrid, mean.hat, col=1, lwd=1.5 )
lines( Tgrid[BLUP.est[IDD,]>0], (BLUP.est[IDD,])[BLUP.est[IDD,]>0], col=1, lty=2, lwd=1.5 )
mtext("B.", cex=1.5, side=3, line=0.7, font=2, at=-0.9)

#  increasing trend  tp go up: max P2
IDD<-  59  
Subject1   <- BMACS[BMACS$ID==IDlevel[IDD],]

plot(Subject1$Time, Subject1$CD4, xlim=c(0,6), ylim=c(0, 66), main= "" , 
     xlab = "Time since infection (years)", ylab = "CD4 percentage")

for (i in unique(BMACS$ID)) { lines(CD4 ~ Time, data = BMACS[BMACS$ID  == i, ], lty = 2, col='gray60')}

points(Subject1$Time, Subject1$CD4, xlim=c(0,6), ylim=c(0, 70),lwd=1.5 )

lines( Tgrid, mean.hat, col=1, lwd=1.5 )
lines( Tgrid, BLUP.est[IDD,], col=1, lty=2, lwd=1.5 )
mtext("C.", cex=1.5, side=3, line=0.7, font=2, at=-0.9)


### sharpest trend to go down
IDD<-  14  #min P2 
Subject1   <- BMACS[BMACS$ID==IDlevel[IDD],]

plot(Subject1$Time, Subject1$CD4, xlim=c(0,6), ylim=c(0, 66), main= "" , 
     xlab = "Time since infection (years)", ylab = "CD4 percentage")

for (i in unique(BMACS$ID)) { lines(CD4 ~ Time, data = BMACS[BMACS$ID  == i, ], lty = 2, col='gray60')}

points(Subject1$Time, Subject1$CD4, xlim=c(0,6), ylim=c(0, 70),lwd=1.5 )

lines( Tgrid, mean.hat, col=1, lwd=1.5 )
lines( Tgrid[BLUP.est[IDD,]>0], (BLUP.est[IDD,])[BLUP.est[IDD,]>0], col=1, lty=2, lwd=1.5 )
mtext("D.", cex=1.5, side=3, line=0.7, font=2, at=-0.9)

dev.off()




















