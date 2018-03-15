### Chapter 14  Ex2: transformation model ###
rm(list = ls())
#----------------
library(npmlda)
library(survival)

#str(NGHS ) #19071*12

NGHS <- data.frame(NGHS, agebin=numeric(nrow(NGHS))) 

for (i in 9:18){
  print(i)
  for (j in 1:10){
    NGHS$agebin[( NGHS$AGE*10 >= (i*10+ (j-1)) &  NGHS$AGE*10 < (i*10+ j))]  <-   (i*10+ (j-1)) 
    print(c(i+ (j-1)/10))
  } }

# summary(NGHS$agebin)

NGHS.sbp <- NGHS[!is.na(NGHS$SBP) & !is.na(NGHS$HTPCT ),] 
dim(NGHS.sbp) #19439 *13

####### From youth BP guidline paper for age 9-17 ###########

AA<- seq(9, 17, by=0.1)
M.SBP.Age<- 102.01027 +  1.94397 *(AA-10) + 0.00598*(AA-10)^2 - 0.00789*(AA-10)^3 -0.00059*(AA-10)^4
SDD<- 10.4855

Q75<-numeric(length(AA)) # Blood pressure quantile base on age for girls with median height

for (i in 1:length(AA))
{  Q75[i] <- qnorm(.75, mean= M.SBP.Age[i], sd = SDD )
}

#############################

NGHS.sbp <-  NGHS.sbp[ NGHS.sbp$agebin<=170,] #15750    13
NGHS.sbp$BP75IND <-NA 
for ( aa in 90: 170)
{
  kk<- aa-89
  NGHS.sbp$BP75IND[NGHS.sbp$agebin==aa ] <- (NGHS.sbp$SBP[NGHS.sbp$agebin==aa ]> Q75[kk])*1
}


# Estimate separately for each race
NGHS.B <- subset(NGHS.sbp, RACE==2)
dim(NGHS.B)  #  8189    13

S1<-10; S2 <- 10 + 3 ;

attach(NGHS.B)

# BP75IND is the indicator of Y(t)> y(75,50) quantile at t #
# agebin ranges from 90-189 to indicate age bins at 9.0, 9.1 ..

Pr.S12 <-   Kernel3D(ID, Y=BP75IND, Time=agebin, X=HTPCT,
            T1=S1*10, T2=S2*10, X0=75, Bndwdth1=15, 
            Bndwdth2= 15, Bndwdth3 =25)


Pr.S2 <-    Kernel3D.S2(ID, Y=BP75IND, Time=agebin, X=HTPCT, 
            T1=S1*10, T2=S2*10, X0=75, Bndwdth1=15,
            Bndwdth2= 15, Bndwdth3 =25)

#Pr.S1= Pr(Y(t)) can be caculated from Chap14_Ex1, then get RTP and RTPR as in Chapter 12.


