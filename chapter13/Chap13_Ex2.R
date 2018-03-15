### Chapter 13 Ex1: transformation model ###
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

#summary(NGHS$agebin)

# --use TVtrans.fit() in Chap13_EX1.R -- #

### LDL ###

NGHS.LDL <- NGHS[!is.na(NGHS$LDL) & !is.na(NGHS$BMI ),] 
dim(NGHS.LDL) #6733 *13

table(NGHS.LDL$agebin)
#  169 170 171  172 173 
#   46  24  16   8   3 
# 174 176 177 178 179 180 181 182 183 184 185 186 187 188 189 
# 1   2   7   11  18  54  44  53  49  83  62  66  75  87  71 
# need to combine some agebins.

NGHS.LDL$agebin[ NGHS.LDL$agebin >=170 & NGHS.LDL$agebin <= 171 ] <- 170.5
NGHS.LDL$agebin[ NGHS.LDL$agebin >=172 & NGHS.LDL$agebin <= 178 ] <- 175
NGHS.LDL$agebin[ NGHS.LDL$agebin >=179 & NGHS.LDL$agebin <= 180 ] <- 179.5

attach(NGHS.LDL)


AgebinLDL <- sort( unique(NGHS.LDL$agebin))
length(AgebinLDL) #92

OR.LDL <- TVtrans.fit(AgebinLDL,  Y=LDL, X1= (RACE==2)*1, X2=BMI)

LDL.beta1.lm <-LocalLm(AgebinTG/10, AgebinTG/10, OR.LDL[,1], bw=3)
LDL.beta2.lm <-LocalLm(AgebinTG/10, AgebinTG/10, OR.LDL[,2], bw=3)

detach(NGHS.LDL)


### TG ##
NGHS.TG <- NGHS[!is.na(NGHS$LDL) & !is.na(NGHS$BMI ),] 
dim(NGHS.TG) #6733 *13

NGHS.TG$agebin[ NGHS.TG$agebin >=170 & NGHS.TG$agebin <= 171 ] <- 170.5
NGHS.TG$agebin[ NGHS.TG$agebin >=172 & NGHS.TG$agebin <= 178 ] <- 175
NGHS.TG$agebin[ NGHS.TG$agebin >=179 & NGHS.TG$agebin <= 180 ] <- 179.5

attach(NGHS.TG)

AgebinTG <- sort( unique(NGHS.TG$agebin))
length(AgebinTG) #92

OR.TG <- TVtrans.fit(AgebinTG,  Y=TG, X1= (RACE==2)*1, X2=BMI)

TG.beta1.lm <-LocalLm(AgebinTG/10, AgebinTG/10, OR.TG[,1], bw=2.7)
TG.beta2.lm <-LocalLm(AgebinTG/10, AgebinTG/10, OR.TG[,2], bw=2.7)

detach(NGHS.TG)
# Generate resampling bootstrap by the same approach from earlier chapters #














