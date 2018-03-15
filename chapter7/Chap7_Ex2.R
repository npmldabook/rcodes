### Chapter 7 Ex2 ###
rm(list = ls())
#----------------
library(npmlda)
library(nlme)
#----------------

NGHS$Black <- (NGHS$RACE==2)*1
NGHS<- NGHS[!is.na(NGHS$SBP) & !is.na(NGHS$BMIPCT) & !is.na(NGHS$HTPCT ),]      
nrow(NGHS) #19320

NGHS$HTPCTc<- NGHS$HTPCT-50
NGHS$BMIPCTc<- NGHS$BMIPCT-50

NGHS$RACE<- NULL
NGHS$Race <- NGHS$Black

NGHS$AGEc<- NGHS$AGE-9
NGHS.fit <- lme(SBP~AGEc+Race+Race:AGEc+HTPCTc+HTPCTc:AGEc
                + BMIPCTc + BMIPCTc:AGEc, random=~1|ID, data=NGHS)
summary(NGHS.fit)
# race effect is not significant #
# Fixed effects: SBP ~ AGEc + Race + Race:AGEc + HTPCTc + HTPCTc:AGEc + BMIPCTc +      BMIPCTc:AGEc 
#                   Value  Std.Error    DF  t-value p-value
# (Intercept)  100.12814 0.21507240 16938 465.5555  0.0000
# AGEc           0.87370 0.02517459 16938  34.7055  0.0000
# Race           0.30588 0.30215527  2374   1.0123  0.3115
# HTPCTc         0.04558 0.00478288 16938   9.5290  0.0000
# BMIPCTc        0.08750 0.00439813 16938  19.8955  0.0000
# AGEc:Race      0.07878 0.03540669 16938   2.2251  0.0261
# AGEc:HTPCTc   -0.00131 0.00065129 16938  -2.0166  0.0438
# AGEc:BMIPCTc   0.00252 0.00062072 16938   4.0556  0.0001


NGHS$AGEc<- NGHS$AGE-18.9
NGHS.fit19 <- lme(SBP~ AGEc+ Black+ Black:AGEc+  HTPCTc +HTPCTc:AGEc+ BMIPCTc+ BMIPCTc:AGEc, random=~1|ID, data=NGHS)
summary(NGHS.fit19)
# When move to 19 as center for age, the race effect is significant #
#                   Value  Std.Error    DF  t-value p-value
# (Intercept)  108.77773 0.22042661 16938 493.4873  0.0000
# AGEc           0.87370 0.02517459 16938  34.7055  0.0000
# Black          1.08582 0.30261959  2374   3.5881  0.0003
# HTPCTc         0.03257 0.00456320 16938   7.1383  0.0000
# BMIPCTc        0.11243 0.00449385 16938  25.0175  0.0000
# AGEc:Black     0.07878 0.03540669 16938   2.2251  0.0261
# AGEc:HTPCTc   -0.00131 0.00065129 16938  -2.0166  0.0438
# AGEc:BMIPCTc   0.00252 0.00062072 16938   4.0556  0.0001