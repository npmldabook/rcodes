### Chapter 5 Ex2 ###
# remove (almost) everything in the working environment.

rm(list = ls())
#----------------
library(npmlda)
#----------------
str(NGHS ) #19071*12
head(NGHS) 

NGHS$Black <- (NGHS$RACE==2)*1

NGHS<- NGHS[!is.na(NGHS$BMI),]

#--get 100 agebins --- #

NGHS <- data.frame(NGHS, agebin=numeric(nrow(NGHS))) 

for (i in 9:18){
  print(i)
  for (j in 1:10){
    NGHS$agebin[( NGHS$AGE*10 >= (i*10+ (j-1)) &  NGHS$AGE*10 < (i*10+ j))]  <-   (i*10+ (j-1)) 
    print(c(i+ (j-1)/10))
  } }
summary(NGHS$agebin) 

NGHS.W <- NGHS[NGHS$Black ==0,]
NGHS.B <- NGHS[NGHS$Black ==1,]

NGHS.sub50w <- NGHS.W[1:413,]
length(table(  NGHS.sub50w$ID))  #50 white girls.

NGHS.sub50B <- NGHS.B[1:411,]
length(table(  NGHS.sub50B$ID))  #50 black girls.

attach(NGHS)
LYM.spl.W <- smooth.spline(agebin[Black ==0 ], BMI[Black ==0 ], all.knots = T, cv=T)#spar=1.023583
LYM.spl.B <- smooth.spline(agebin[Black ==1], BMI[Black ==1 ], all.knots = T, cv=T) #spar= 1.035914 

Weight.spl.W <- smooth.spline(NGHS.W$agebin, NGHS.W$BMI, all.knots = T, cv=T, w=1/NGHS.W$Freq)
Weight.spl.B <- smooth.spline(NGHS.B$agebin, NGHS.B$BMI, all.knots = T, cv=T, w=1/NGHS.B$Freq)

##bootstrap CI ##

IDlist.W <- unique(NGHS.W$ID)
IDlist.B <- unique(NGHS.B$ID)

nID.W <-  length(table(NGHS.W$ID))  #1164 girls.
nID.B <- length(table(NGHS.B$ID)) # 1213 girls

Bootsample.W <- function(){ 
  resample.ID.W <- sample(IDlist.W ,nID.W ,replace=T) 
  do.call("rbind", lapply(1:nID.W, function(i) subset(NGHS.W[,c('ID','agebin','BMI')] , ID==resample.ID.W[i])))
}


Bootsample.B <- function(){ 
  resample.ID.B <- sample(IDlist.B ,nID.B ,replace=T) 
  do.call("rbind", lapply(1:nID.B, function(i) subset(NGHS.B[,c('ID','agebin','BMI')], ID==resample.ID.B[i])))
}

set.seed(119)

Time.int<- 90:189
NN<- length(Time.int) ; NN#100

nBoot<- 1000
BootFIT.W <- matrix(NA, nrow= nBoot, ncol= NN)
BootFIT.B<- matrix(NA, nrow= nBoot, ncol= NN)

for (j in 1:nBoot)
{
  if ((j-floor(j/10)*10)== 0 ) print(j)
  
  BootdataW<- Bootsample.W()
  BootdataB<- Bootsample.B()
  
  # fit and predict, first estimate the constant r and convert to get the same lamda each time #
   testW <- smooth.spline(BootdataW$agebin, BootdataW$BMI, all.knots = T, spar= 1.023583, cv=NA)
   logrW <- log(testW$lambda)-(3*testW$spar-1) *log(256)
   sparW <- ((log(1.759866)- logrW )/log(256) +1)/3
  
  spl.W <- smooth.spline(BootdataW$agebin, BootdataW$BMI, all.knots = T, spar=sparW , cv=NA)
  BootFIT.W[j,] <-  predict(spl.W , (90:189))$y    
  
   testB <- smooth.spline(BootdataB$agebin, BootdataB$BMI, all.knots = T, spar= 1.035914, cv=NA)
   logrB <- log(testB$lambda)-(3*testB$spar-1) *log(256)
   sparB <- ((log(2.31726)- logrB )/log(256) +1)/3
  
  spl.B <- smooth.spline(BootdataB$agebin, BootdataB$BMI, all.knots = T, spar=sparB , cv=NA)
  BootFIT.B[j,] <-  predict(spl.B, (90:189))$y    
}

UpperCI.W <-  apply( BootFIT.W,  2, quantile,.975 )
LowerCI.W <-  apply( BootFIT.W,  2, quantile,.025 )

UpperCI.B <-  apply( BootFIT.B,  2, quantile,.975 )
LowerCI.B <-  apply( BootFIT.B,  2, quantile,.025 )


#------fig5.3 ---------#

postscript("fig5.3.BMI.ps", horizontal=T)

par(mar=c(28.5, 4.5, 3, 1),cex=2, mfrow=c(1, 3), pty='m', cex.lab=1.4, cex.axis=1.3,  cex.main=1.4)

# plot 1: Black
plot(NGHS.sub50B$agebin, NGHS.sub50B$BMI , main="BMI: African American girls",ylim=c(10,50),
     xlab = "", ylab = "", cex=0.9, pch=1,  axes=F)

  grid( lwd=2, col='gray60')
  axis(2)
  axis(1, at=seq(90, 190, by=20), labels=seq(9,19,2));box()
  lines(predict(LYM.spl.B , (90:189)), col = 1,lwd=2)
  mtext(expression(paste("Body Mass Index(kg/ ",m^2,")")), side=2, font=1,line=2.2,cex=1.1) # ylabel
  mtext("Age (years)", cex=1, side=1, line=2.4)
  mtext("A.", cex=1.5, side=3, line=0.7, font=2, at=72)

# plot 2: White

plot(NGHS.sub50w$agebin, NGHS.sub50w$BMI ,   main="BMI: Caucasian girls",
     xlab = "", ylab = "", cex=0.9, pch=1, ylim=c(10,50),  axes=F, col='gray40')
  grid( lwd=2 , col='gray60')
  axis(2)
  axis(1, at=seq(90, 190, by=20), labels=seq(9,19,2));box()
  lines(predict(LYM.spl.W , (90:189)), lwd=2 , col='gray40')
  
  mtext(expression(paste("Body Mass Index(kg/ ",m^2,")")), side=2, font=1,line=2.2,cex=1) # ylabel
  mtext("Age (years)", cex=1, side=1, line=2.4)
  mtext("B.", cex=1.5, side=3, line=0.7, font=2, at=72)

### plot3 : CI 

plot(NGHS.sub50B$agebin, NGHS.sub50B$BMI , main="mean BMI curves", type='n',ylim=c(16,28), 
     xlab = "", ylab ="", cex=0.8, pch=1,  axes=F)

   lines(predict(LYM.spl.W , (90:189)), lwd=2, col='gray40')
   grid(lwd=2,  col='gray60')
   lines(90:189, UpperCI.W, col='gray40',lwd=1.5, lty=2)
   lines(90:189, LowerCI.W, col='gray40',lwd=1.5, lty=2)
  
   lines(predict(LYM.spl.B , (90:189)), col = 1,lwd=2)
   lines(90:189, UpperCI.B, col = 1,lwd=1.5, lty=2)
   lines(90:189, LowerCI.B, col = 1,lwd=1.5, lty=2)
  
   axis(2)
   axis(1, at=seq(90, 190, by=20), labels=seq(9,19,2)); box()
  
   mtext(expression(paste("Body Mass Index(kg/ ",m^2,")")), side=2, font=1,line=2.2,cex=1) # ylabel
   mtext("Age (years)", cex=1, side=1, line=2.4)
  
   legend('bottomright', lty=1, lwd=2, col=c('black', 'gray40') , title="smoothing spline fit",
         legend=c("African American girls","Caucasian girls"), bty='n',   cex=1.3)
  
   mtext("C.", cex=1.5, side=3, line=0.7, font=2, at=72)


 dev.off()
#-------
 detach(NGHS)



















