### Chapter 6, Ex1 ###
## remove (almost) everything in the working environment.
rm(list = ls())

#library
library (npmlda) 
data(BMACS)

### Fig 6.1 analysis ##
## center Age ##
Sample.Age <- as.numeric(unlist( with(BMACS, by(age ,ID , mean)))) 
mean.Age   <- mean(Sample.Age) #[1] 34.18154
BMACS$ageC <- BMACS$age- mean.Age

## center preCD4 ##
Sample.preCD4 <- as.numeric(unlist( with(BMACS, by(preCD4,ID , mean)))) 
mean.preCD4 <- mean(Sample.preCD4 ) #[1] 42.91935
BMACS$preCD4C  <- BMACS$preCD4- mean.preCD4

Ct <- data.frame(table(BMACS$ID))
Ct$IDD<- 1:nrow(Ct)
names(Ct)<- c("ID", "ni", "IDD")
BMACS<- merge(BMACS, Ct, by= "ID")

# Obtain the sample size, N=283 for BMACS
 N <- length(unique(BMACS$ID))
# Obtain the baseline covariate:
# first observation per subject ID
 Xi <- do.call("rbind", as.list(by(BMACS[,c("preCD4C", "Smoke","ageC")], BMACS$ID, head, n=1)))
# Formula (6.4-6.5)
 Xi <- as.matrix(cbind(1,Xi))
 EXX <- Reduce("+", lapply(1:N, function(i) Xi[i,] %*% t(Xi[i,])))
 EstInv <- solve(EXX/N)
 
# Obtain the four pseudo longitudinal samples
 eX <- data.frame(IDD=1:N,
                 eX1 = apply(t(t(Xi)*EstInv[1,]), 1, sum),
                 eX2 = apply(t(t(Xi)*EstInv[2,]), 1, sum),
                 eX3 = apply(t(t(Xi)*EstInv[3,]), 1, sum),
                 eX4 = apply(t(t(Xi)*EstInv[4,]), 1, sum))
 BMACS <- merge(BMACS, eX, by="IDD")
 BMACS$PseudoY <- with(BMACS, cbind(eX1, eX2, eX3, eX4)*CD4)

# apply the kernel estimate as in Chapter 3.

 FitL1.Ni <- with(BMACS,  kernel.fit(sort(unique(Time)), Time, PseudoY[,1], bw=1.5, Kernel="Nm", Wt=1/ni))
 FitL2.Ni <- with(BMACS,  kernel.fit(sort(unique(Time)), Time, PseudoY[,2], bw=1.5, Kernel="Nm", Wt=1/ni))
 FitL3.Ni <- with(BMACS,  kernel.fit(sort(unique(Time)), Time, PseudoY[,3], bw=1.5, Kernel="Nm", Wt=1/ni))
 FitL4.Ni <- with(BMACS,  kernel.fit(sort(unique(Time)), Time, PseudoY[,4], bw=1.5, Kernel="Nm", Wt=1/ni))
 
 Time.int<-  sort(unique(BMACS$Time))
 
 
 # run similar bootstrap as in Chapters 3-5 to get the confidence intervals#
 # bootstrap sampler #
 
 Bootsample <- function(){ 
   resample.ID <- sample(x= unique(BMACS$ID) ,size= N ,replace=T) 
   do.call("rbind", lapply(1:N , function(i) data.frame(subset(BMACS, ID==resample.ID[i]),ID2=i ) ))
 }

 
 #------Figure 6.1--------#
 postscript("fig6.1.ps", horizontal=T)
 
 par(mar=c(4.5, 4.5, 3, 1),cex=2, mfrow=c(2, 2), pty='m', cex.lab=1.2, cex.axis=1.2, cex.main=1.3)

 plot(Time.int, FitL1.Ni ,  lwd=2, xlab = "Years", ylab = "baseline CD4",   cex=0.7, col=1, main="baseline CD4", ylim=c(10,40), type='l') 
  lines(Time.int,UpperCI1 , lwd=1, lty=2, col=1)
  lines(Time.int,LowerCI1 , lwd=1, lty=2,col=1)
   mtext("A.", cex=1.5, side=3, line=0.7, font=2, at=-0.9)
 
  plot(Time.int, FitL2.Ni ,  lwd=2, xlab = "Years", ylab = "Coefficient",   cex=0.7, col=1, main="Pre-infection CD4",  type='l', ylim=c(-1,1)) 
   lines(Time.int,UpperCI2 , lwd=1, lty=2,col=1)
   lines(Time.int,LowerCI2 , lwd=1, lty=2,col=1)
   mtext("B.", cex=1.5, side=3, line=0.7, font=2, at=-0.9)
 
  plot(Time.int, FitL3.Ni ,  lwd=2,  xlab = "Years", ylab = "Coefficient",   cex=0.7, col=1, main="Smoking",  type='l',ylim=c(-20,20) ) 
   lines(Time.int,UpperCI3 , lwd=1,lty=2, col=1)
   lines(Time.int,LowerCI3 , lwd=1, lty=2,col=1)
   mtext("C.", cex=1.5, side=3, line=0.7, font=2, at=-0.9)
 
  plot(Time.int, FitL4.Ni ,   xlab = "Years", ylab = "Coefficient",  lwd=2, cex=0.7, col=1, main="Age",  type='l',ylim=c(-1.5,1.5) ,axes=F) 
   lines(Time.int,UpperCI4 , lwd=1, lty=2,col=1)
   lines(Time.int,LowerCI4 , lwd=1, lty=2, col=1)
    mtext("D.", cex=1.5, side=3, line=0.7, font=2, at=-0.9)
   axis(2, at=seq(-1.5,1.5, by=0.5),  labels= seq(-1.5,1.5, by=0.5) )
   axis(1); box()
  dev.off()
 
 ### Fig 6.2 analysis: using smooth.spline as in Chapter 5. ##
  
  FitL1.Ni.sm  <- smooth.spline(BMACS$Time, BMACS$PseudoY[,1] , spar=0.95,  all.knots = T,  cv=F, w= 1/BMACS$ni) 
  #lamda 0.8398814
  
  FitL2.Ni.sm  <- smooth.spline(BMACS$Time, BMACS$PseudoY[,2] , spar=1.05, all.knots = T,  cv=F, w= 1/BMACS$ni) 
  #lamda  4.43292
 
  FitL3.Ni.sm  <- smooth.spline(BMACS$Time, BMACS$PseudoY[,3] , spar=0.97, all.knots = T,  cv=F, w= 1/BMACS$ni)
  #lamda  1.171419
  
  FitL4.Ni.sm  <- smooth.spline(BMACS$Time, BMACS$PseudoY[,4] , spar=1.06, all.knots = T,  cv=F, w= 1/BMACS$ni) #
  #lambda= 5.235247

  Pred.FitL1.Ni <- predict(FitL1.Ni.sm, Time.int)$y 
  Pred.FitL2.Ni <- predict(FitL2.Ni.sm, Time.int)$y 
  Pred.FitL3.Ni <- predict(FitL3.Ni.sm, Time.int)$y 
  Pred.FitL4.Ni <- predict(FitL4.Ni.sm, Time.int)$y 
  


