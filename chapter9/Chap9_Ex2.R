### Chapter 9 Ex2 ###
rm(list = ls())
#----------------
library(npmlda)
library(splines)

# -- test stat, full model-- #
#-- Code as in Example 1 --#

 fit.WLS <- lm(SBP ~ 0 + Bs.age + Bs.race + Bs.HT + Bs.BMI,
                 weights=1/ni, data=NGHS)
 summary(fit.WLS)

 # -- test stat, full model-- #
 Ntotal<- nrow(Ct) #2376
 RSS1<- sum((1/NGHS$ni)* (fit.WLS$residuals)^2)/ Ntotal #70.02297
 
 #  model A:
 fit.WLS.A1  <-  lm(SBP ~ 1,  weights=1/ni , data=NGHS)  
 summary( fit.WLS.A1 )
 RSS.A1<- sum((1/NGHS$ni)* (fit.WLS.A1$residuals)^2)/ Ntotal # 88.34
(final.T.A1<- (RSS.A1-RSS1)/RSS1) # 0.2616372
 
 
 # model B:
 fit.WLS.A2  <-  lm(SBP ~ 0+ Bs.age,  weights=1/ni , data=NGHS)  
 summary( fit.WLS.A2 )
 RSS.A2<- sum((1/NGHS$ni)* (fit.WLS.A2$residuals)^2)/ Ntotal # 79.52248
 (final.T.A2<- (RSS.A2-RSS1)/RSS1) # [1] 0.1356628
  
 # model C:  beta1= race no effect  , 3 cov in the model
 #------------
 fit.WLS1  <-  lm(SBP ~ 0+ Bs.age+  Bs.HT  + Bs.BMI,  weights=1/ni , data=NGHS)  
 summary( fit.WLS1)
 RSS01<- sum((1/NGHS$ni)* (fit.WLS1$residuals)^2)/ Ntotal # 70.15608
 (final.Tn1<- (RSS01-RSS1)/RSS1) # 0.001900925
 
 
 ## model D:  beta2 = height no effect  , 3 cov in the model
 #------------
 fit.WLS2  <-  lm(SBP ~ 0+ Bs.age+ Bs.race+  Bs.BMI,  weights=1/ni , data=NGHS)  
 summary( fit.WLS2)
 RSS02<- sum((1/NGHS$ni)* (fit.WLS2$residuals)^2)/ Ntotal # 71.50392
 ( final.Tn2<- (RSS02-RSS1)/RSS1 ) # 0.02114957
 
 
 ## model E:  beta3 = BMI no effect  
 #------------
 fit.WLS3  <-  lm(SBP ~ 0+ Bs.age+ Bs.race+  Bs.HT ,  weights=1/ni , data=NGHS)  
 summary( fit.WLS3)
 RSS03<- sum((1/NGHS$ni)* (fit.WLS3$residuals)^2)/ Ntotal # 76.091
 ( final.Tn3<- (RSS03-RSS1)/RSS1 ) # 0.08665768
 
 
 ## model F:  beta1  = race is constant
 #------------
 fit.WLS4  <-  lm(SBP ~ 0+ Bs.age+  Black + Bs.HT +  Bs.BMI,  weights=1/ni , data=NGHS)  
 summary( fit.WLS4)
 ( RSS04<- sum((1/NGHS$ni)* (fit.WLS4$residuals)^2)/ Ntotal ) # 70.05739
 ( final.Tn4<- (RSS04-RSS1)/RSS1 )  #0.0004915256
 
 
 ## model G:  beta2  = height is constant
 #------------
 
 fit.WLS5  <-  lm(SBP ~ 0+ Bs.age+ Bs.race  + HTPCTc  + Bs.BMI ,  weights=1/ni , data=NGHS)  
 summary( fit.WLS5)
 ( RSS05<- sum((1/NGHS$ni)* (fit.WLS5$residuals)^2)/ Ntotal ) #70.12891
 ( final.Tn5<- (RSS05-RSS1)/RSS1 ) #0.001512951
 
 
 
 ## model H:  beta3  = BMI is constant:  
 #------------
 
 fit.WLS6  <-  lm(SBP ~  0+ Bs.age+ Bs.race  + Bs.HT  + BMIPCTc  ,  weights=1/ni , data=NGHS)  
 summary( fit.WLS6)
 ( RSS06<- sum((1/NGHS$ni)* (fit.WLS6$residuals)^2)/ Ntotal ) #70.12011
 
 ( final.Tn6<- (RSS06-RSS1)/RSS1 ) # 0.001387336
 
## Print the test statistics ##
  t(cbind(final.T.A1,  final.T.A2, final.Tn1,final.Tn2,final.Tn3,final.Tn4,final.Tn5,final.Tn6))
 
 NGHS2<- NGHS
 NGHS2<- data.frame(NGHS,
                    pseudoYA1 =  fit.WLS.A1$fitted +  fit.WLS$resid,
                    pseudoYA2 =  fit.WLS.A2$fitted +  fit.WLS$resid,
                    pseudoY1 =  fit.WLS1$fitted +  fit.WLS$resid,
                    pseudoY2 =  fit.WLS2$fitted +  fit.WLS$resid,
                    pseudoY3 =  fit.WLS3$fitted +  fit.WLS$resid,
                    pseudoY4 =  fit.WLS4$fitted +  fit.WLS$resid,
                    pseudoY5 =  fit.WLS5$fitted +  fit.WLS$resid,
                    pseudoY6 =  fit.WLS6$fitted +  fit.WLS$resid)
 
 head(NGHS2,3) ; NGHS2$SBP<- NULL
 
 #--------the bootstrap sampler -------------#
 N <- 2376
 Bootsample <- function(){ 
   resample.ID <- sample(x= unique(NGHS2$ID) ,size= N ,replace=T) 
   do.call("rbind", lapply(1:N , function(i) data.frame(subset(NGHS2, ID==resample.ID[i]),ID2=i ) ))
 }
 
 
 #--------the bootstrap sampler -------------#
 set.seed(827)
 nBoot<- 5000
 
 Boot.Coef  <- matrix(NA, nrow= nBoot, ncol= 8)
 for (i in 1:nBoot)
 {
   if ((i-floor(i/100)*100)== 0 ) print(i) 
   Bootdata <- Bootsample()
   
   Bs.age  <- bs(Bootdata$AGE, knots=KN0 , degree=3,  intercept=T, Boundary.knots=c(9, 19)) 
   Bs.race <- bs(Bootdata$AGE, knots=NULL,  degree=3,  intercept=T, Boundary.knots=c(9, 19))* Bootdata$Black  
   Bs.HT   <- bs(Bootdata$AGE, knots=NULL,  degree=3,  intercept=T, Boundary.knots=c(9, 19))* Bootdata$HTPCTc
   Bs.BMI  <- bs(Bootdata$AGE, knots=NULL,  degree=3,  intercept=T, Boundary.knots=c(9, 19))* Bootdata$BMIPCTc   
   
   
   ## model A:   intercept only model , no covariates
   
   fit.WLS.full  <-  lm(Bootdata$pseudoYA1 ~ 0+ Bs.age+ Bs.race  + Bs.HT  + Bs.BMI,  weights=1/Bootdata$ni)  
   RSS1.A1 <- sum((1/Bootdata$ni)* (fit.WLS.full$resid)^2)/ Ntotal #
   
   fit.WLS.A1  <-  lm(pseudoYA1 ~ 1,  weights=1/ni , data=Bootdata)  
   RSS0A1<- sum((1/Bootdata$ni)* (fit.WLS.A1$residuals)^2)/ Ntotal # 
   Tn.A1 <- (RSS0A1-RSS1.A1 )/RSS1.A1 
   
   
   ## model B:  baseline only but time-varing, no covariates
   #------------
   fit.WLS.full2  <-  lm(Bootdata$pseudoYA2 ~ 0+ Bs.age+ Bs.race  + Bs.HT  + Bs.BMI,  weights=1/Bootdata$ni)  
   RSS1.A2 <- sum((1/Bootdata$ni)* (fit.WLS.full2$resid)^2)/ Ntotal #
   
   fit.WLS2.A2  <-  lm(pseudoYA2 ~ 0+ Bs.age,  weights=1/ni , data=Bootdata)  
   RSS0A2<- sum((1/Bootdata$ni)* (fit.WLS2.A2$residuals)^2)/ Ntotal # 
   Tn.A2 <- (RSS0A2-RSS1.A2)/ RSS1.A2
   
   ## model C:  beta1= race no effect  , 3 cov in the model
   #------------
   fit.WLS  <-  lm(Bootdata$pseudoY1 ~ 0+ Bs.age+ Bs.race  + Bs.HT  + Bs.BMI,  weights=1/Bootdata$ni)  
   RSS1<- sum((1/Bootdata$ni)* (fit.WLS$residuals)^2)/ Ntotal #
   
   fit.WLS1  <-  lm(pseudoY1 ~ 0+ Bs.age+  Bs.HT  + Bs.BMI,  weights=1/ni , data=Bootdata)  
   RSS01<- sum((1/Bootdata$ni)* (fit.WLS1$residuals)^2)/ Ntotal # 
   Tn1<- (RSS01-RSS1)/RSS1 
   
   
   ## model D:  beta2 = height no effect  , 3 cov in the model
   #------------
   fit.WLS  <-  lm(Bootdata$pseudoY2 ~ 0+ Bs.age+ Bs.race  + Bs.HT  + Bs.BMI,  weights=1/Bootdata$ni)  
   RSS1<- sum((1/Bootdata$ni)* (fit.WLS$residuals)^2)/ Ntotal #
   
   
   fit.WLS2  <-  lm(pseudoY2 ~ 0+ Bs.age+ Bs.race+  Bs.BMI,  weights=1/ni , data=Bootdata)  
   RSS02<- sum((1/Bootdata$ni)* (fit.WLS2$residuals)^2)/ Ntotal # 
   Tn2<- (RSS02-RSS1)/RSS1  
   
   
   ## model E:  beta3 = BMI no effect  , 3 cov in the model
   #------------
   fit.WLS  <-  lm(Bootdata$pseudoY3 ~ 0+ Bs.age+ Bs.race  + Bs.HT  + Bs.BMI,  weights=1/Bootdata$ni)  
   RSS1<- sum((1/Bootdata$ni)* (fit.WLS$residuals)^2)/ Ntotal #
   
   
   fit.WLS3  <-  lm(pseudoY3 ~ 0+ Bs.age+ Bs.race+  Bs.HT ,  weights=1/ni , data=Bootdata)  
   RSS03<- sum((1/Bootdata$ni)* (fit.WLS3$residuals)^2)/ Ntotal # 
   Tn3<- (RSS03-RSS1)/RSS1  
   
   
   ## model F:  beta1  = race is constant: 4 cov in model 
   #------------
   fit.WLS  <-  lm(Bootdata$pseudoY4 ~ 0+ Bs.age+ Bs.race  + Bs.HT  + Bs.BMI,  weights=1/Bootdata$ni)  
   RSS1<- sum((1/Bootdata$ni)* (fit.WLS$residuals)^2)/ Ntotal #
   
   fit.WLS4  <-  lm(pseudoY4 ~ 0+ Bs.age+  Black + Bs.HT +  Bs.BMI,  weights=1/ni , data=Bootdata)  
   RSS04<- sum((1/Bootdata$ni)* (fit.WLS4$residuals)^2)/ Ntotal 
   Tn4<- (RSS04-RSS1)/RSS1  
   
   
   ## model G:  beta2  =HT is constant: 4 cov in model 
   #------------
   fit.WLS  <-  lm(Bootdata$pseudoY5 ~ 0+ Bs.age+ Bs.race  + Bs.HT  + Bs.BMI,  weights=1/Bootdata$ni)  
   RSS1<- sum((1/Bootdata$ni)* (fit.WLS$residuals)^2)/ Ntotal #
   
   fit.WLS5  <-  lm(pseudoY5 ~ 0+ Bs.age+ Bs.race  + HTPCTc  + Bs.BMI ,  weights=1/ni , data=Bootdata)  
   RSS05<- sum((1/Bootdata$ni)* (fit.WLS5$residuals)^2)/ Ntotal 
   Tn5<- (RSS05-RSS1)/RSS1  
   
   
   ## model H:  beta3  = BMI is constant: 4 cov in model 
   #------------
   fit.WLS  <-  lm(Bootdata$pseudoY6 ~ 0+ Bs.age+ Bs.race  + Bs.HT  + Bs.BMI,  weights=1/Bootdata$ni)  
   RSS1<- sum((1/Bootdata$ni)* (fit.WLS$residuals)^2)/ Ntotal #
   
   fit.WLS6  <-  lm(pseudoY6 ~  0+ Bs.age+ Bs.race  + Bs.HT  + BMIPCTc  ,  weights=1/ni , data=Bootdata)  
   
   RSS06<- sum((1/Bootdata$ni)* (fit.WLS6$residuals)^2)/ Ntotal  #70.12011
   Tn6<- (RSS06-RSS1)/RSS1 
   
   Boot.Coef[i,] <- c(Tn.A1,  Tn.A2, Tn1,Tn2,Tn3,Tn4,Tn5,Tn6)
   
 }
 
 
 final.T<- c(final.T.A1, final.T.A2, final.Tn1,final.Tn2,final.Tn3,final.Tn4,final.Tn5,final.Tn6)
 
 nBoot<-5000
 sum(final.T.A1 <= Boot.Coef[,1])/nBoot # 0
 sum(final.T.A2 <= Boot.Coef[,2])/nBoot # 0
 sum(final.Tn1  <= Boot.Coef[,1])/nBoot # 0.20
 sum(final.Tn2  <= Boot.Coef[,2])/nBoot # 0
 sum(final.Tn3  <= Boot.Coef[,3])/nBoot # 0
 sum(final.Tn4  <= Boot.Coef[,4])/nBoot # 0.23
 sum(final.Tn5  <= Boot.Coef[,5])/nBoot # 0.0076
 sum(final.Tn6  <= Boot.Coef[,6])/nBoot # 0
 
 #summary(Boot.Coef)
 
 