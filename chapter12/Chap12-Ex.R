### Chapter 12 Ex ###
rm(list = ls())
#----------------
library(npmlda)
#----------------
str(NGHS ) #19071*12

NGHS$Black <- (NGHS$RACE==2)*1
NGHS<- NGHS[!is.na(NGHS$BMI),]
dim(NGHS) #19398 obs

#--get 100 agebins --- #

NGHS <- data.frame(NGHS, agebin=numeric(nrow(NGHS))) 

for (i in 9:18){
  print(i)
  for (j in 1:10){
    NGHS$agebin[( NGHS$AGE*10 >= (i*10+ (j-1)) &  NGHS$AGE*10 < (i*10+ j))]  <-   (i*10+ (j-1)) 
    print(c(i+ (j-1)/10))
  } }
NGHS$agebin<- NGHS$agebin/10
summary(NGHS$agebin) 

NGHS.W <- NGHS[NGHS$Black ==0,]
NGHS.B <- NGHS[NGHS$Black ==1,]

IDD <- unique(NGHS.B$ID)
nID <- length(IDD)
nID # no. of subjects [1] 1213

nrow(NGHS.B) # no. of visits [1] 10028

Grid <- 0.2
S1cat <- seq(9.0, 18.8, by=Grid)
S12cat <- seq(9.0, 16.8, by=Grid)

## Compute the PA(x,s1)##
KK1 <- length(S1cat)
Yvec<- (NGHS.B$BMIPCT>= 85)*1
Xvec<-  NGHS.B$agebin

Prob.S1 <- numeric(KK1)

for(k in 1:KK1)
 {
   Prob.S1[k]<-NW.WtKernel(Xvec, Yvec, S1cat[k], Bndwdth =2.6)
}

## Compute the PA(x,s1,s2)##
KK2 <- length(S12cat) 
Prob.S1S2 <- numeric(KK2)
IDD <- NGHS.B$ID

for(k in 1:KK2)
 {
   #print(k)
   Prob.S1S2[k]<- Kernel2D(IDD, Xvec, Yvec, X01=S12cat[k],X02=S12cat[k]+2, Bndwdth1=0.9, Bndwdth2=0.9)
 }

## Compute RTP and RTPR ##
IND1 <- (S1cat >=9 & S1cat <=16.8)
IND2 <- (S1cat >=11 & S1cat <=18.8)
RTP <- Prob.S1S2/(Prob.S1[IND1])
RTPR <- RTP/(Prob.S1[IND2])

#  summary(RTP)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.7872  0.8709  0.8812  0.8768  0.8904  0.9115 
#  summary(RTPR)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.151   2.198   2.238   2.239   2.270   2.335 

### For white we can run similarly with H1=3.5,  hh1=hh2=0.8 ##
### Run bootstrap to get confidence interval ##




