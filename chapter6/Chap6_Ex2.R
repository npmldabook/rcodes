### Chapter 6, Ex2 : Simulation  ###
## remove (almost) everything in the working environment.
rm(list = ls())

# -library-
library (npmlda) 
library(MASS)

#--We consider a design that is similar to the nature of the BMACS CD4 study #
Sigma.i <- matrix(NA, 31, 31)
diag(Sigma.i)<-0.0625 # sd=0.25 for each obs

for (i in 1:31)
{  for (j in 1:31)
 {
   Sigma.i[i,j] <- 0.0625*exp(-abs(i-j))
 }
}

#sum( Sigma.i-t( Sigma.i)) #0  # symmetric#

set.seed(119)
Npt<-  400
Nobs<- 31  
Ntotal<- Npt*Nobs # 12400

ii<- rep(Nobs, Npt); 
id1<-rep(1:Npt,ii) ; 

Time<- 0:30
Time.Pt <- rep(Time, Npt)
# by pt data, initial same time grid, then randomly drop 60% #

Nsim  <-  1000  
Sm.SimFIT1  <-Sm.SimFIT2  <-Sm.SimFIT3  <- matrix(NA, nrow= Nsim, ncol= 31)
Ker.SimFIT1 <- Ker.SimFIT2 <- Ker.SimFIT3  <- matrix(NA, nrow= Nsim, ncol= 31)

for (j in 1: Nsim)
{
  if ((j-floor(j/20)*20)== 0 ) print(j) 
  
  # ---For each simulated sample ---
  #--- Variable X1, time invariant---
  
  X1 <- rbinom(Npt, 1, 0.5)  # X1: binary baseline
  X1.Pt <- rep(X1 , ii)       # by pt data, repeat baseline 
  
  X2  <- rnorm( Npt, mean = 0, sd = 4)    # X2: Normal (0, sd=4)
  X2.Pt <- rep(X2 , ii)                    # by pt data, repeat baseline 
  
  #--------Coefficient curve ----- # 
  #---- baseline curve beta0 ---  
  
  beta0t <- 3.5 + 6.5*sin(Time.Pt*pi/60)   
  
  #--- parameter curve b1 for X1 ---
  
  beta1t <-  -0.2 - 1.6*cos((Time.Pt-30)*pi/60)
  
  #--- parameter curve b2 for X2---
  
  beta2t <-  0.25 - 0.0074*((30-Time.Pt)/10)^3
  
  #--- generate  random errors: mvrnorm in MASS library
  
  Error <-  as.vector(t(mvrnorm(Npt, mu=rep(0, 31), Sigma= Sigma.i)))
  YY  <-  ( beta0t +  beta1t*  X1.Pt  + beta2t* X2.Pt )+ Error 
  
  # finish data generation #
  ID <-  rbinom(Ntotal, 1, 0.4)  # binary indicator to include or not
  SimData <- (data.frame(ID =id1,  Time= Time.Pt,  X1 =X1.Pt, X2  = X2.Pt,  Y=YY))[ID==1,]
  
  ## summary(as.vector(table(SimData$ID))) # mean is 12.5, which is 40% OK ##
  
  ##-----------Start with data fitting --------------------##
  # Formula (6.4-6.5)
  Xi <-  as.matrix(cbind(1,X1,X2))  #400* 3
  EXX<-  Reduce('+', lapply(1:Npt, function(i) Xi[i,] %*% t(Xi[i,]))) # 3*3
  
  EstInv  <- solve(EXX /Npt)  
  
  #Obtain pseudo sample 
  eX <- data.frame(ID=1:Npt, eX1 = apply(t(t(Xi)*EstInv[1,]), 1, sum),
                   eX2 = apply(t(t(Xi)*EstInv[2,]), 1, sum),
                   eX3 = apply(t(t(Xi)*EstInv[3,]), 1, sum) )
  
  
  SimData <- merge(SimData , eX , by='ID')
  SimData$PseudoY <- with(SimData , cbind(eX1,eX2,eX3)*Y )  # this is a matrix, not a vector in 1 col
  
  ## get Ni weights ##
  Ct <-   data.frame(table(SimData$ID))
  names(Ct)<- c('ID', 'Ni')
  SimData<- merge(SimData, Ct, by= 'ID') # Good this is a fast way to get the counts#
  
  #-----------smooth.spline ----------
  
  FitL1   <-  smooth.spline(SimData$Time, SimData$PseudoY[,1] ,  spar= 0.675, cv=NA, w= 1/SimData$Ni, all.knots = T)
  FitL2   <- smooth.spline(SimData$Time, SimData$PseudoY[,2]  ,  spar= 0.824, cv=NA, w= 1/SimData$Ni, all.knots = T)
  FitL3   <- smooth.spline(SimData$Time, SimData$PseudoY[,3]  ,  spar= 0.769, cv=NA, w= 1/SimData$Ni, all.knots = T)
  
  ### for kernel estimate ##
  
  Ker.SimFIT1[j,]  <- kernel.fit(Time, SimData$Time, SimData$PseudoY[,1], bw= 2 ,Kernel="Nm", Wt= 1/SimData$Ni)
  Ker.SimFIT2[j,]  <- kernel.fit(Time, SimData$Time, SimData$PseudoY[,2], bw= 3 ,Kernel="Nm", Wt= 1/SimData$Ni)
  Ker.SimFIT3[j,]  <- kernel.fit(Time, SimData$Time, SimData$PseudoY[,3], bw= 3 ,Kernel="Nm", Wt= 1/SimData$Ni)
  
  ####
  
  Sm.SimFIT1[j,] <- predict(FitL1  , Time)$y
  Sm.SimFIT2[j,] <- predict(FitL2  , Time)$y
  Sm.SimFIT3[j,] <- predict(FitL3  , Time)$y
  
}

## Save data if needed ##
# write.table(Sm.SimFIT1,  file = "Sm.SimFIT1.txt", append = FALSE)
# write.table(Sm.SimFIT2,  file = "Sm.SimFIT2.txt", append = FALSE)
# write.table(Sm.SimFIT3,  file = "Sm.SimFIT3.txt", append = FALSE)
# 
# write.table(Ker.SimFIT1,  file = "Ker.SimFIT1.txt", append = FALSE)
# write.table(Ker.SimFIT2,  file = "Ker.SimFIT2.txt", append = FALSE)
# write.table(Ker.SimFIT3,  file = "Ker.SimFIT3.txt", append = FALSE)


Sm.UpperCI1 <-  apply( Sm.SimFIT1,  2, quantile,.975 )
Sm.LowerCI1 <-  apply( Sm.SimFIT1,  2, quantile,.025 )
Sm.mean1    <-  apply( Sm.SimFIT1,  2, mean)

Sm.UpperCI2 <-  apply( Sm.SimFIT2,  2, quantile,.975 )
Sm.LowerCI2 <-  apply( Sm.SimFIT2,  2, quantile,.025 )
Sm.mean2    <-  apply( Sm.SimFIT2,  2, mean)


Sm.UpperCI3 <-  apply( Sm.SimFIT3,  2, quantile,.975 )
Sm.LowerCI3 <-  apply( Sm.SimFIT3,  2, quantile,.025 )
Sm.mean3    <-  apply( Sm.SimFIT3,  2, mean)


Ker.UpperCI1 <-  apply( Ker.SimFIT1,  2, quantile,.975 )
Ker.LowerCI1 <-  apply( Ker.SimFIT1,  2, quantile,.025 )
Ker.mean1    <-  apply( Ker.SimFIT1,  2, mean)

Ker.UpperCI2 <-  apply( Ker.SimFIT2,  2, quantile,.975 )
Ker.LowerCI2 <-  apply( Ker.SimFIT2,  2, quantile,.025 )
Ker.mean2    <-  apply( Ker.SimFIT2,  2, mean)

Ker.UpperCI3 <-  apply( Ker.SimFIT3,  2, quantile,.975 )
Ker.LowerCI3 <-  apply( Ker.SimFIT3,  2, quantile,.025 )
Ker.mean3    <-  apply( Ker.SimFIT3,  2, mean)

#--- fig 6.3 --- #
Time<- 0:30

# baseline curve beta0:  #------
beta0t.True <- 3.5 + 6.5*sin(Time*pi/60)   # baseline beta0

#--- parameter curve b1 for X1 
beta1t.True <-  -0.2 - 1.6*cos((Time-30)*pi/60)

#--- parameter curve b2  for X2
beta2t.True <-  0.25 - 0.0074*((30-Time)/10)^3


postscript("fig6.3.ps", horizontal=T)

par(mar=c(4.5, 5.5, 3,1),cex=2,mfrow=c(1,3), pty='s', cex.lab=1.6, cex.axis=1.4, cex.main=1.6)
#----------
plot(Time, beta0t.True , lwd=1, type='l', col='gray40', xlab="Time", ylab= expression(beta[0](t)), main="", axes=F, ylim=c(3,10))
axis(2, las=2, at= seq(2,12,2))
axis(1)
box()

lines(Time , Ker.mean1,   col= 1,  lty=3, lwd=2)
lines(Time , Sm.mean1,   col= 1, lwd=1, lty=2)
mtext("A.", cex=1.3, side=3, line=0.7, font=2, at=-6.2)


#---------
plot(Time, beta1t.True , type='l',  xlab="Time",col='gray40', ylab= expression(beta[1](t)), main="", axes=F, ylim=c(-2.2,0.2))
axis(2, las=2,mgp=c(3,0.5,0), tck=-0.015)
axis(1)
box()

lines(Time , Ker.mean2,   col= 1,  lty=3, lwd=2)
lines(Time , Sm.mean2,   col= 1, lwd=1, lty=2)
mtext("B.", cex=1.3, side=3, line=0.7, font=2, at=-6.2)

#---------

plot(Time, beta2t.True , type='l', xlab="Time",col='gray40', ylab= expression(beta[2](t)), main="", axes=F, ylim=c(0,0.3))
axis(2, las=2,mgp=c(3,0.5,0), tck=-0.015)
axis(1)
box()

lines(Time , Ker.mean3,   col= 1,  lty=3, lwd=2)
lines(Time , Sm.mean3,   col= 1, lwd=1, lty=2)
mtext("C.", cex=1.3, side=3, line=0.7, font=2, at=-6.2)

dev.off()


#############






















