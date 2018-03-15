### Chapter 7 Ex3 ###
rm(list = ls())
#----------------
library(npmlda)
library(nlme)
#----------------


CD4fit <- lme(CD4~ Time+ preCD4 + Smoke +age, random=~Time|ID, data=BMACS)
summary(CD4fit)


Sample.Age <- as.numeric(unlist( with(BMACS, by(age ,ID , mean)))) 
mean.Age   <- mean(Sample.Age) #[1] 34.2
BMACS$ageC <- BMACS$age- mean.Age

Sample.preCD4 <- as.numeric(unlist( with(BMACS, by(preCD4,ID , mean)))) 
mean.preCD4 <- mean(Sample.preCD4 ) #[1] 42.91935
BMACS$preCD4C  <- BMACS$preCD4- mean.preCD4


Ct <-   data.frame(table(BMACS$ID))
Ct$IDD<- 1:nrow(Ct)

names(Ct)<- c('ID', 'ni','IDD')
BMACS<- merge(BMACS, Ct, by= 'ID') 


Time.grid<-  seq(0.1, 5.9 , by=0.2)
Beta <- with(BMACS, LocalLm.Beta( Time.grid, Time, X1=preCD4C, X2=Smoke, 
                                  X3= ageC, CD4, Bndwdth=1.8,  Weight = 1/ni))

##### Follow similar steps as Example 1 for bootstrap CI  ########
