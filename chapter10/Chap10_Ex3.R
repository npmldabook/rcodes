## Code for Ch 10, Sec 10.7.2: full model ##
# Adapted from code authored by: Dr. Wenhua Jiang
# Reference: Wu, C. O., Tian, X. and Jiang, W. A shared parameter model for the
#   estimation of longitudinal concomitant intervention effects Biostatistics,
#   12(4):737–749, 2011.

rm(list = ls())
#----------------

library(npmlda)
library(nlme)
library(MASS)
library(splines)

# functions of the full model #
source("ENRICHD_full.R")

#####################
Ct <-   data.frame(table(BDIdata$ID))
names(Ct)<- c('ID', 'count')
BDIdata<- merge(BDIdata, Ct, by= 'ID')

medall<- BDIdata[, c(1,3,2,5,6)]
names(medall)<- c("ID",	"BDIrec",	"recday",	"medday",	"count")
medall <- as.matrix(medall)

s <- 1
medall.id <- rep(NA, 1000)
medall.id[1] <- medall[1, 1]
measurement <- rep(NA, 1000)
measurement[1] <- sum(medall[, 1] == medall.id[1])

for (i in 2:dim(medall)[1]) {
  if (medall[i, 1] != medall[i - 1, 1]) {
    s <- s + 1
    medall.id[s] <- medall[i, 1]
    measurement[s] <- sum(medall[, 1] == medall.id[s])
  }
}
medall.id <- medall.id[1:s]   # ID list
measurement <- measurement[1:s] # Frequency per ID
length(medall.id) #557
max(measurement) # 36
min(measurement) # 5
summary(measurement)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 5.00    8.00   12.00   12.78   16.00   36.00 

## Results for the full model ##
set.seed(101)
BDIdata <- medall

Results.full <- SPM.full(data=BDIdata, n.iter=10, n.boot=1000)
# Note: To save running time, we can divide up no. of bootstraps 
# and run  mutiple programs in computer clusters 


## Results for the Sub model ##
source("ENRICHD_sub.R")
set.seed(101)

Results.sub <- SPM.sub(data=BDIdata, n.iter=10, n.boot=1000)









