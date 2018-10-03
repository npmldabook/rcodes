### Chapter 1###

## remove (almost) everything in the working environment.
rm(list = ls())
#----------------
#install.packages("npmlda")
library(npmlda) # install from CRAN

# Alternative or updated libary in Github
library(devtools)
install_github("npmldabook/npmlda")

library(splines)
library(lme4)

#------ Figure 1.1  -----------#
data(BMACS)
dim(BMACS) #1817    6
head(BMACS)
tail(BMACS)
names(BMACS) 
dimnames(BMACS)

length(table(BMACS$ID)) # 283 unique ID

summary(as.vector(table(BMACS$ID)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   1.00    3.00    6.00    6.42   10.00   14.00 

IDD<- unique(BMACS$ID) # 283 unique ID
BMACS$CD4p<- BMACS$CD4/100


postscript("fig1.1.ps", horizontal=T)

par(mar=c(4.5, 4.5, 3, 1),cex=2,mfrow=c(1,1), pty='m', cex.lab=1.8, cex.axis=1.5)

plot(CD4p ~ Time, data = BMACS,  xlab = "Time since infection (years)", ylab = "CD4 percentage", type = "n", ylim=c(0,0.8)) # type='n' , not plot dots
 for (i in unique(BMACS$ID)) { lines(CD4p ~ Time, data = BMACS[BMACS$ID  == i, ], lty = 1)}

dev.off()

#------ Figure 1.2  -----------#

#------ Figure 1.2 ----------#
str(NGHS)

NGHS <-  NGHS[!is.na(NGHS$BMI  ),]
dim(NGHS) # 19398    12

# pick first 150 girls
NGHS.sub100 <- cbind(NGHS[NGHS$ID<=150,])

### randomly select three girls ##

NGHS.sub<- subset(NGHS, ID %in% c(1,1916, 2225))


postscript("fig1.2-layout.ps", horizontal=F)

layout(rbind(c(1, 1, 1, 2, 2, 2) , c(3,3,4, 4,5, 5), c(6,6, 7 , 7, 8, 8)),  respect=F )

par( mar=c(4.3, 4.7, 4, 2), cex.lab=1.5, cex.axis=1.3,pty='m')

plot(BMI~ AGE, data =  NGHS.sub100 ,  xlab = "Age (years)", ylab = expression(paste("Body Mass Index(kg/ ",m^2,")")), cex=0.4) 
mtext("A.", cex=1.4, side=3, line=0.1, font=2, at=7)

plot(SBP ~ AGE, data =  NGHS.sub100 ,  xlab = "Age (years)", ylab = "Systolic Blood Pressure (mm Hg)", cex=0.4) 

par(mar=c(3, 3, 2, 1), pty='m', cex.lab=1.1, cex.axis=0.9, cex=1, cex.main=1.1, font.main=1)

##---BMI --- ## 
plot( BMI~AGE,  NGHS.sub[NGHS.sub$ID==1,] , type='o', ylab="BMI", xlab="", main="Girl 1", ylim=c(10, 50), axes=F)
axis(2,  mgp=c(3,0.5,0),tck=-0.04) 
axis(1, mgp=c( 3,0.5,0),tck=-0.04)
box()
mtext("BMI", side=2, font=1,line=2,cex=1.1) # ylabel
mtext("B.", cex=1.4, side=3, line=0.1, font=2, at=6.8)

plot( BMI~AGE,  NGHS.sub[NGHS.sub$ID==2225,]  , type='o', ylab="", xlab="", main="Girl 2", ylim=c(10, 50),axes=F)
axis(2,  mgp=c(3,0.5,0),tck=-0.04) 
axis(1, mgp=c( 3,0.5,0),tck=-0.04)
box()

plot( BMI~AGE,  NGHS.sub[NGHS.sub$ID==1916,]  , type='o', ylab="", xlab="", main="Girl 3", ylim=c(10, 50),axes=F)
axis(2,  mgp=c(3,0.5,0),tck=-0.04) 
axis(1, mgp=c( 3,0.5,0),tck=-0.04)
box()

##---SBP --- ## 

par(mar=c(4, 3, 2, 1), pty='m', cex.lab=1.1, cex.axis=0.9,cex=1, cex.main=1.1, font.main=1)

plot( SBP~AGE,  NGHS.sub[NGHS.sub$ID==1,] , type='o', ylab="SBP", xlab="Age (years)", main="Girl 1", ylim=c(80, 160),axes=F)
mtext("C.", cex=1.4, side=3, line=0.1, font=2, at=6.8)
mtext("SBP", side=2, font=1,line=2,cex=1.1) # ylabel
axis(2,  mgp=c(3,0.5,0),tck=-0.04) 
axis(1, mgp=c( 3,0.5,0),tck=-0.04)
box()

plot( SBP~AGE,  NGHS.sub[NGHS.sub$ID==2225,]  ,type='o',  ylab="", xlab="Age (years)", main="Girl 2", ylim=c(80, 160),axes=F)
axis(2,  mgp=c(3,0.5,0),tck=-0.04) 
axis(1, mgp=c( 3,0.5,0),tck=-0.04)
box() 


plot( SBP~AGE,  NGHS.sub[NGHS.sub$ID==1916,]  , type='o',  ylab="", xlab="Age (years)", main="Girl 3", ylim=c(80,160),axes=F)
axis(2,  mgp=c(3,0.5,0),tck=-0.04) 
axis(1, mgp=c( 3,0.5,0),tck=-0.04)
box()

dev.off()

#------ Figure 1.3  -----------#
## BDI before and after medication curves ##
## check BDI score, change, before and after med ##

data(BDIdata)
str(BDIdata)
BDIdata$Tijm <- BDIdata$time*12/365.25
BDIsub <- subset(BDIdata, med.time >=0 & med.time < 200)
dim(BDIsub)# 1465    6


length(unique(BDIdata$ID)) #557
length(unique(BDIsub$ID))  #92

summary(BDIsub$med.time[BDIsub$med.time>0] )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 10.0    28.0    76.0    73.4   111.0   172.0 


BDIsub$time <- BDIsub$Tij*365.25/12

BDIsub$visitday2<-  (BDIsub$time- BDIsub$med.time)


postscript("fig1.3.ps", horizontal=T)

par(mar=c(4.5, 4.5, 3, 1),cex=0.5,mfrow=c(1,2), pty='m', cex.lab=1.6, cex.axis=1.2)

 plot(BDIsub$time,  BDIsub$BDI,  xlab = "Days on study", ylab = "BDI score", ylim=c(0, 60), cex=0.5)
 for (i in unique(BDIsub$ID)) { lines(BDI ~ time, data = BDIsub[BDIsub$ID  == i, ], lty = 1)}

 # align before and after medication

 plot(BDIsub$visitday2,  BDIsub$BDI,  xlab = "Days since medication", ylab = "BDI score", ylim=c(0, 60), cex=0.5)

 for (i in unique(BDIsub$ID)) { lines(BDI ~ visitday2, data = BDIsub[BDIsub$ID  == i, ], lty = 1)}
 abline(v=0,lty=2)

dev.off()

#---------Figure 1.4 -------------#

data(HSCT)
str(HSCT)

summary(HSCT$Granu)
summary(HSCT$LYM)
summary(HSCT$MON)


HSCT$Granu.log <- log10(HSCT$Granu)
HSCT$LYM.log   <- log10(HSCT$LYM)
HSCT$MON.log   <- log10(HSCT$MON)

HSCT$'G-CSF'[ HSCT$'G-CSF'==0 ] <- NA
HSCT$'IL-15'[ HSCT$'IL-15'==0 ] <- NA
HSCT$'MCP-1'[ HSCT$'MCP-1'==0 ] <- NA

HSCT$GCSF.log <- log10(HSCT$'G-CSF')
HSCT$IL15.log   <- log10(HSCT$'IL-15')
HSCT$MCP1.log   <- log10(HSCT$'MCP-1')


postscript("fig1.4.ps", horizontal=T)

par(mar=c(4.5, 5.5, 3,1),cex=2,mfrow=c(2,3), pty='m', cex.lab=1.6, cex.axis=1.4, cex.main=1.6)

#-----------PMN-----------

plot(Granu.log ~  Days, data=HSCT,  main='Granulocyte',  col='black', pch=1, cex=1,
     xlab='', ylab='',  xlim=c(-8,30) , ylim=c(-3,1.5), axes=F)

axis(2,at= -3:1,  labels= c(1, 10,100,1000,10000), cex.axis=1.6 , las=2 )
axis(1, at=seq(-7, 28, by=7),cex.axis=1.6) 
box()

fit<- loess(Granu.log~ Days, data=HSCT )
lines(seq(-8, 35, 1), predict(fit, data.frame(Days  = seq(-8, 35, 1))),col=4, lwd=2) 

#-----------LYM  -------------

plot(LYM.log ~  Days,  data=HSCT, main='Lymphocyte',  col='black', pch=1, cex=1,
     xlab='', ylab='',  xlim=c(-8,30) , ylim=c(-3,0.5), axes=F)

axis(2,at= -3:0,  labels= c(1, 10,100,1000), cex.axis=1.6 , las=2 )
axis(1, at=seq(-7, 28, by=7),cex.axis=1.6) 
box()

fit<- loess(LYM.log~ Days, data=HSCT)
lines(seq(-8, 35, 1), predict(fit, data.frame(Days  = seq(-8, 35, 1))),col=4, lwd=2) 


#-----------MON -------------

plot(MON.log  ~  Days,  data=HSCT,    main='Monocyte',  col='black', pch=1, cex=1,
     xlab='', ylab='',  xlim=c(-8,30) , ylim=c(-3,0.5), axes=F)


axis(2,at= -3:0,  labels= c(1, 10,100,1000), cex.axis=1.6 , las=2 )
axis(1, at=seq(-7, 28, by=7),cex.axis=1.6) 
box()

fit<- loess(MON.log~ Days, data=HSCT)
lines(seq(-8, 35, 1), predict(fit, data.frame(Days  = seq(-8, 35, 1))),col=4, lwd=2) 


## now three cytokines ###

#-----------GCSF -------------

plot(GCSF.log   ~  Days,  data=HSCT,   main='G-CSF',  col='black', pch=1, cex=1,
     xlab='Days post-transplant', ylab='',  xlim=c(-8,30) , ylim=c(-0.5,4), axes=F)


axis(2,at= 0:4,  labels= c(1, 10,100,1000,10000), cex.axis=1.6 , las=2 )
axis(1, at=seq(-7, 28, by=7),cex.axis=1.6) 
box()


fit<- loess(GCSF.log  ~ Days, data=HSCT)
lines(seq(-8, 35, 1), predict(fit, data.frame(Days  = seq(-8, 35, 1))),col=4, lwd=2) 


#-----------IL15 -------------

plot(IL15.log  ~  Days,  data=HSCT,    main='IL-15',  col='black', pch=1, cex=1,
     xlab='Days post-transplant', ylab='',  xlim=c(-8,30) , ylim=c(0,3), axes=F)


axis(2,at= 0:3,  labels= c(1, 10,100,1000), cex.axis=1.6 , las=2 )
axis(1, at=seq(-7, 28, by=7),cex.axis=1.6) 
box()


fit<- loess(IL15.log~ Days, data=HSCT)
lines(seq(-8, 35, 1), predict(fit, data.frame(Days  = seq(-8, 35, 1))),col=4, lwd=2) 


#-----------MCP1 -------------

plot(MCP1.log   ~  Days,  data=HSCT,    main="MCP-1",  col='black', pch=1, cex=1,
     xlab='Days post-transplant', ylab='',  xlim=c(-8,30) , ylim=c(0,4), axes=F)


axis(2,at= 0:4,  labels= c(1, 10,100,1000,10000), cex.axis=1.6 , las=2 )
axis(1, at=seq(-7, 28, by=7),cex.axis=1.6) 
box()

fit<- loess(MCP1.log ~ Days, data=HSCT)
lines(seq(-8, 35, 1), predict(fit, data.frame(Days  = seq(-8, 35, 1))),col=4, lwd=2) 


dev.off()


### plot the smoothing line with ggplot2 ##
library(ggplot2)
library(grid)
library(gridExtra)

pdf(file="fig1.4b.pdf")
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,3 )))  

p1<- ggplot(HSCT, aes(x=Days  , y=Granu.log  )) + geom_point() +  geom_smooth( stat = "smooth", size=1)+ xlab("")+ylab("")  + ggtitle((title = 'Granulocyte'))
p2<- ggplot(HSCT, aes(x=Days  , y=LYM.log  )) + geom_point() +  geom_smooth(size=1)+ xlab("")+ylab("")  + ggtitle((title = "Lymphocyte"))
p3<- ggplot(HSCT, aes(x=Days  , y=MON.log  )) + geom_point() +  geom_smooth(size=1)+ xlab("")+ylab("")  + ggtitle((title = "Monocyte"))

p4<- ggplot(HSCT, aes(x=Days  , y=GCSF.log   )) + geom_point() +  geom_smooth(size=1)+ xlab("")+ylab("")  + ggtitle((title = 'G-CSF'))
p5<- ggplot(HSCT, aes(x=Days  , y=IL15.log  )) + geom_point() +  geom_smooth(size=1)+ xlab("")+ylab("")  + ggtitle((title = "IL-15"))
p6<- ggplot(HSCT, aes(x=Days  , y=MCP1.log  )) + geom_point() +  geom_smooth(size=1)+ xlab("")+ylab("")  + ggtitle((title = "MCP-1"))


print(p1+theme_bw(),  vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p2+theme_bw(),  vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p3+theme_bw(),  vp=viewport(layout.pos.row = 1, layout.pos.col = 3))

print(p4+theme_bw(),  vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(p5+theme_bw(),  vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
print(p6+theme_bw(),  vp=viewport(layout.pos.row = 2, layout.pos.col = 3))
dev.off()

#------ Figure 1.5  Kernel function -----------#
#1.1 uniform distribution KU(s)
S <- seq(-2, 2, by =0.01)
KU.S <- numeric(length(S))  
KU.S[S>= -1 & S<= 1] <- 1/2


#1.2 Epanechnikov
S <- seq(-2, 2, by =0.01)
KE.S <- numeric(length(S))  
KE.S[S>= -1 & S<= 1] <-  3/4*(1-S[S>= -1 & S<= 1]^2) 

# 1.3 trangular 
S <- seq(-2, 2, by =0.01)
KT.S <- numeric(length(S))  
KT.S[S>= -1 & S<= 1] <-  (1- abs(S[S>= -1 & S<= 1]))

#1.4  biweight 
S <- seq(-2, 2, by =0.01)
KQ.S <- numeric(length(S))  
KQ.S[S>= -1 & S<= 1] <-  15/16*(1- (S[S>= -1 & S<= 1])^2)^2

#1.5  Tricube

S <- seq(-2, 2, by =0.01)
KC.S <- numeric(length(S))  
KC.S[S>= -1 & S<= 1] <-  70/81* (1- (abs(S[S>= -1 & S<= 1]))^3)^3

#1.6 Gaussian

S2 <- seq(-3, 3, by =0.01)
KG.S <-  dnorm(S2, 0, 1)# Gaussian density

#--- Figure 1.5 ---- #
par(mfrow=c(3,2), mar=c(4.2, 4.4, 2.5, 1),cex=0.8, font.lab=1, cex.axis=1.1, font.axis=1, cex.main=1.2, cex.lab=1.4)
#1
plot(S, KU.S, type='s', xlim=c(-2, 2),  xlab='s', ylab=expression(K[U](s)),main='Uniform')

#2. Epanechnikov
plot(S, KE.S, type='s', xlim=c(-2, 2), xlab='s', ylab=expression(K[E](s)),main='Epanechnikov')


#3. Triangular
plot(S, KT.S, type='s', xlim=c(-2, 2), xlab='s', ylab=expression(K[T](s)),main='Triangular')


#4. Biweight
plot(S, KQ.S, type='s', xlim=c(-2, 2), xlab='s', ylab=expression(K[Q](s)),main='Quartic or Biweight')

# 5.  Tricube kernel
plot(S, KC.S, type='s', xlim=c(-2, 2), xlab='s', ylab=expression(K[C](s)),main='Tricube')


#6.  Gaussian, the equivalent bandwidth is the sigma= standard deviation=1
plot(S2, KG.S, type='s', xlim=c(-3, 3), xlab='s', ylab=expression(K[G](s)),main='Gaussian')
















