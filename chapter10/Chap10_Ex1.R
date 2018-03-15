rm(list = ls())
#----------------

library(nlme) # or use lme4 
install.packages("geepack")
library(geepack)
#----------------
# A Synthetic Example of Concomitant Intervention #

n <-  24  # sample sizes
m1<-  10  # of visit for each subjects 

ii<- rep(m1,n); 
id1<-rep(1:n,ii) 
#length(id1) # 240 

s    <- c(rep(2,120), rep(8, 120))   # 12 obs in 2 and 12 obs in 8 
time <-  rep(1:10,24)                # time from 1 to 10
delta <- ifelse(time > s ,1,0)       # indicator of med   
y1    <- c(rep(c(rep(20,2),rep(19, 8)),12),rep(c(rep(19,8),rep(17, 2)),12))

set.seed(111)
y<- y1+ rnorm(10*24, 0, 3)

fulldata<-cbind(id1,s,time,delta,y)
fulldata<-data.frame(fulldata)

# make the sample plot of the data example ###
testdata1 <- cbind(c(1,2,3,3:10), c( rep(20,3),rep(19, 8)) )
testdata2 <- cbind(c(1:9, 9, 10), c(rep(19,9),rep(17, 2)) )

plot(testdata1[,1],  testdata1[,2], type="l",xlab="time", ylab="Y", main="A Simple Example", ylim=c(16,20),lwd=3 )
lines(testdata2[,1], testdata2[,2],  col=4, lwd=3)


############################################################
##### without S
###LM#: LM (working indep) 
summary(lm(y ~ delta, data=fulldata))

###LME:  random intercept
summary(lme(y~ delta, data=fulldata, random=~1 |id1))

###LME:  random intercept & slope
summary(lme(y~ delta, data=fulldata, random=~1+delta |id1))


############################################################
##### without S ## GEE ###
###GEE (working indep) 
summary(geese(y ~  delta , id=id1, data=fulldata, corstr="independence"))

###GEE : exchangable
summary(geese(y ~  delta , id=id1, data=fulldata, corstr="exchangeable"))

###GEE:  unstructured
summary(geese(y ~  delta , id=id1, data=fulldata, corstr="unstructured"))



############################################################
### adding S in the model 
###LM: LME (working indep) 
summary( lm(y ~ s+delta, data=fulldata) )

###LME:  random intercept
summary(lme(y~ s+ delta, data=fulldata, random=~1 |id1))

###LME:  random intercept & slope
summary(lme(y~ s+ delta, data=fulldata, random=~1+delta |id1))

#-----------------------------------------------------------------------------
## Run 1000 Simulation  to get average of the mean and SE estimate for beta, then
## get coverage probability for CI= mean+/-1.96*SE for above 9 methods (Table 10.1)
#-----------------------------------------------------------------------------


