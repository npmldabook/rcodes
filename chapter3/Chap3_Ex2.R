### Chapter 3, Ex2 ###

## remove (almost) everything in the working environment.
rm(list = ls())

#library
library (npmlda) 
# need to use LocalLm() and CVlm() from npmlda#

### Fig 3.2 analysis ##
fit.linear.1 <- loess(CD4 ~ Time , span =0.1 , degree=1, data=BMACS)
fit.linear.5  <- loess(CD4 ~ Time , span =0.5 , degree=1, data=BMACS)

fit.Qr.1 <- loess(CD4 ~ Time , span =0.1 , degree=2,  data=BMACS)
fit.Qr.5 <- loess(CD4 ~ Time , span =0.5,  degree=2,  data=BMACS)

Time.int<- seq(0.1,5.9,  by=0.1)

postscript("fig3.2.ps", horizontal=T)

par(mar=c(4.5, 4.5, 3, 1),cex=2, mfrow=c(2, 2), pty='m', cex.lab=1.2, cex.axis=1.2, cex.main=1.2)

plot(CD4 ~ Time, data=BMACS,  xlab = "Time since infection (years)", ylab = "CD4 percentage",  ylim=c(0,65), cex=0.7, col='gray50', main="Local linear: span=0.1") 
lines(Time.int, predict(fit.linear.1 , data.frame(Time  = Time.int)),col=4, lwd=2)

plot(CD4 ~ Time, data=BMACS,  xlab = "Time since infection (years)", ylab = "CD4 percentage",  ylim=c(0,65), cex=0.7, col='gray50', main="Local linear: span=0.5") 
lines(Time.int, predict(fit.linear.5 , data.frame(Time  = Time.int)),col=4, lwd=2)


plot(CD4 ~ Time, data=BMACS,  xlab = "Time since infection (years)", ylab = "CD4 percentage",  ylim=c(0,65), cex=0.7, col='gray50', main="Local quadratic: span=0.1") 
lines(Time.int, predict(fit.Qr.1 , data.frame(Time  = Time.int)),col=4, lwd=2)

plot(CD4 ~ Time, data=BMACS,  xlab = "Time since infection (years)", ylab = "CD4 percentage",  ylim=c(0,65), cex=0.7, col='gray50', main="Local quadratic: span=0.5") 
lines(Time.int, predict(fit.Qr.5 , data.frame(Time  = Time.int)),col=4, lwd=2)

dev.off()

### Fig 3.3 analysis ##
# Obtain CV bandwidth
data(BMACS)

Ct <- data.frame(table(BMACS$ID))
Ct <- data.frame( Ct, 1:nrow(Ct) )
names(Ct)<- c("ID", "ni", "IDD")
BMACS<- merge(BMACS, Ct, by= "ID")

#str(BMACS)
#head(BMACS,30)

hh<- seq(0.5, 3, by =.1)
Nh<- length(hh)  
CVh<- numeric(Nh)

for (j in 1:Nh)
{
  #print(j)
  CVh[j] <-  with(BMACS, CVlm( Time,  CD4 , bw= hh[j],  IDD, Wt=1/ni))
}

plot(hh, CVh,type='l' )
hh[CVh==min(CVh)]  # h(CV) = 0.9

# Obtain bootstrap CI using local linear fit with CV.bandwidth 
Time.int<- seq(0.1,5.9,  by=0.1)
LocalFit.Y <- with(BMACS, LocalLm(Time.int, Time, CD4, bw=0.9, Wt=1/ni))

#--------
IDlist <- unique(BMACS$ID)
nID    <- length(IDlist)

# Bootstrap function
Bootsample <- function(){
  resample.ID <- sample(IDlist ,nID ,replace=T)
  do.call("rbind", lapply(1:nID,
                          function(i) subset(BMACS, ID==resample.ID[i])))}

# Obtain fitted value at a time grid
LocalLm.Fit<- function(Data, Time.int){
  with(Data, LocalLm(Time.int, Time,CD4,bw=0.9, Wt=1/ni))}

# Compute the 95% CI based on B=1000 bootstrap replicates
Boot.Fit <- replicate(1000, LocalLm.Fit(Bootsample(), Time.int))
UpperCI <- apply(Boot.Fit, 1, quantile,.975)
LowerCI <- apply(Boot.Fit, 1, quantile,.025)


# plot with local linear fit and 95% CI ##
plot(CD4 ~ Time, data = BMACS,
      xlab = "Time since infection(years)",
      ylab = "CD4 percentage", cex=0.3, col="gray70", main="")
polygon(c(Time.int[1], Time.int, rev(Time.int)),
         c(LowerCI[1], UpperCI, rev(LowerCI)), col="gray60", border=NA)
lines(Time.int, LocalFit.Y, lwd=2.5, col=1)







