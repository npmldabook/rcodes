### Chapter 1###

## remove (almost) everything in the working environment.
rm(list = ls())


# library ("XX")  # need to set up later
# data("BMACS")   # can make sure it can be loaded later.

library(splines)
library(lme4)

setwd("C:/Users/tianx/Documents/1tianx/1R-code/Data")

BMACS<-  read.table("BMACS.txt"  ,header=T)[,2:7]

dim(BMACS) #1817    6

head(BMACS)
tail(BMACS)

   ID Time Smoke   age preCD4 CD4
1 1022  0.2     0 26.25     38  17
2 1022  0.8     0 26.25     38  30
3 1022  1.2     0 26.25     38  23
4 1022  1.6     0 26.25     38  15
5 1022  2.5     0 26.25     38  21
6 1022  3.0     0 26.25     38  12

names(BMACS) 
#[1] "ID"     "Time"   "Smoke"  "age"    "preCD4" "CD4" (CD4 percentage)

#dimnames(BMACS)

length(table(BMACS$ID)) # 283 unique ID

# check # of repeasted 

summary(as.vector(table(BMACS$ID)))

 Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1.00    3.00    6.00    6.42   10.00   14.00 
> 

IDD<- unique(BMACS$ID) # 283 unique ID

BMACS$CD4p<- BMACS$CD4/100


postscript("fig1.1.ps", horizontal=T)

par(mar=c(4.5, 4.5, 3, 1),cex=2,mfrow=c(1,1), pty='m', cex.lab=1.8, cex.axis=1.5)

plot(CD4p ~ Time, data = BMACS,  xlab = "Time since infection (years)", ylab = "CD4 percentage", type = "n", ylim=c(0,0.8)) # type='n' , not plot dots
 for (i in unique(BMACS$ID)) { lines(CD4p ~ Time, data = BMACS[BMACS$ID  == i, ], lty = 1)}

dev.off()

#or above: make it brief, also plot(y~x) or plot(x,y) are the same 

plot(CD4/100 , Time, data = BMACS,  xlab = "Time since infection (years)", ylab = "CD4 percentage", type = "n", ylim=c(0,0.8))
 for (i in unique(BMACS$ID)) { lines(CD4p ~ Time, data = BMACS[BMACS$ID  == i, ], lty = 1)}


######### the following code has no use since pre-CD4 is not baseline#######
k<-1
BMACS2<- BMACS
for (i in unique(BMACS$ID))
{
  print(k)
  datai<- BMACS[BMACS$ID  == i, ]
  newid<-  data.frame(ID=i,  Time=0, Smoke= datai$Smoke[1], age= datai$age[1], preCD4=NA, CD4= datai$preCD4[1])
  BMACS2<- rbind(BMACS2, newid)
  k<- k+1
}

BMACS2[1800:2100,]

#order by ID and time, drop preCD4.

BMACS2<- BMACS2[order(BMACS2$ID, BMACS2$Time),]

BMACS2<- BMACS2[,-5]


BMACS2$CD4p<- BMACS2$CD4/100

# plot from time 0.
postscript("fig1.1b.ps", horizontal=T)

par(mar=c(4.5, 4.5, 3, 1),cex=2,mfrow=c(1,1), pty='m', cex.lab=1.8, cex.axis=1.5)

plot(CD4p ~ Time, data = BMACS2,  xlab = "Time since infection (years)", ylab = "CD4 percentage", type = "n", ylim=c(0,0.8)) # type='n' , not plot dots
 for (i in unique(BMACS2$ID)) { lines(CD4p ~ Time, data = BMACS2[BMACS2$ID  == i, ], lty = 1)}

dev.off()

save(BMACS2, BMACS, file = "BMACS.CD4.RData") 


######### the following code has no use since pre-CD4 is not baseline#######


### Rice & Wu CD4 counts data ##


CD4<- read.table("CD4data.txt", header=F)
names(CD4) <- c("time", "CD4", "age", "packs", "DrugUse", "NoSex", "cesd", "ID")

dim(CD4) #2376 values 8 
head(CD4)
tail(CD4)
length(table(CD4$ID)) # 369


      time CD4  age packs DrugUse NoSex cesd    ID
1 -0.741958 548 6.57     0       0     5    8 10002
2 -0.246407 893 6.57     0       1     5    2 10002
3  0.243669 657 6.57     0       1     5   -1 10002
4 -2.729637 464 6.95     0       1     5    4 10005
5 -2.250513 845 6.95     0       1     5   -4 10005
6 -0.221766 752 6.95     0       1     5   -5 10005

#time since seroconversion
#CD4 count
#age (relative to arbitrary origin)
#packs of cigarettes smoked per day
#recreational drug use yes/no
#number of sexual partners
#cesd (mental illness score)
#subject ID


with(CD4[CD4$time>=0,],   plot(time*12,CD4,  xlab = "time since seroconversion (month)", ylab = "CD4 Counts"  ,type='n' ))

CD4after<- CD4[CD4$time>=0,] 
CD4after$month <- CD4after$time*12

#par(mar=c(4.5, 4.5, 3, 1),cex=2,mfrow=c(1,2), pty='m', cex.lab=1.8, cex.axis=1.5)

for (i in unique(CD4after$ID)) { lines(CD4 ~ month , data = CD4after[CD4after$ID  == i, ], lty = 1)}

#########################################################################################
####  Example 1.2 NGHS : BMI  and SBP example, overall data and selected Sample ###

 # all NGHS has 20900 obs, missing 11 age, then cut age at 9-19 to get 19701. (1188 obs). min age was 8.96 continous.
 # NGHS <- NGHS[NGHS$AGE >=9  & NGHS$AGE <= 19,] 
 # dim(NGHS) 
 # [1] 19701     10


 NGHS <- read.table("NGHS10var.txt", header= T)
  dim(NGHS ) #19071*11
  head(NGHS) 
 row.names     ID RACE      AGE SYSAV DIA4AV DIA5AV HEIGHT WEIGHT   BMIPCT    HTPCT agebin      BMI
1         1 100011    2  9.35000    92     63     56 140.20  34.80 69.42891 80.28091     93 17.70448
2         2 100011    2 10.34634    84     72     68 148.50  38.95 59.74340 88.92655    103 17.66260
3         3 100011    2 11.28268   101     62     57 157.10  44.10 53.78723 93.38943    112 17.86842
4         4 100011    2 12.33676   100     NA     68 163.90  51.15 60.23078 92.70794    123 19.04091
5         5 100011    2 13.38809    99     81     77 166.45  57.40 70.56797 87.48605    133 20.71783
6         6 100011    2 14.42847   100     67     56 165.95  57.80 66.59429 76.79805    144 20.98811

  names(NGHS)

  sapply(NGHS, class)
row.names        ID      RACE       AGE     SYSAV    DIA4AV    DIA5AV    HEIGHT    WEIGHT    BMIPCT     HTPCT 
"integer" "integer" "integer" "numeric" "numeric" "numeric" "numeric" "numeric" "numeric" "numeric" "numeric" 
   agebin       BMI 
"numeric" "numeric"

  length(table(NGHS$ID)) #2378.
   summary(NGHS$AGE) #9.00 - 19.00 to get exacly 100 bins ###

 # plot BMI and SBP scatter plot

 summary(NGHS$AGE) #9.00 - 19.00 to get exacly 100 bins ###
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   # 9.00   11.65   13.90   14.05   16.50   19.00 

 AgeBins<- seq(9, 18.9, by=.1) #100 bins

 NGHS <- data.frame(NGHS, agebin=numeric(nrow(NGHS))) 

 for (i in 9:18){
	 print(i)
    for (j in 1:10){
     NGHS$agebin[( NGHS$AGE*10 >= (i*10+ (j-1)) &  NGHS$AGE*10 < (i*10+ j))]  <-   (i*10+ (j-1)) 
	  print(c(i+ (j-1)/10))
	} }

   summary(NGHS$agebin)  
  # The agebin/10= agebin name, *10 makes it interger to avoid floating point error later
  #Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 90     116     139     140     165     189 

 
### BMI information ###
## Colin pediatric paper:   Normal (<85th percentile); overweight(85<= . < 95th); obese (>95 th percentile)

NGHS$BMI  <-  NGHS$WEIGHT/(NGHS$HEIGHT/100)^2
NGHS <-  NGHS[!is.na(NGHS$BMI  ),]
 dim(NGHS) # 19398    13
 NGHS.sub100 <- cbind(NGHS[1:1228,])
 length(table(NGHS.sub100$ID))  #150 girls.

NGHS$Black <- (NGHS$RACE==2)*1
with(NGHS, table(Black , RACE))


# sort NGHS by ID and age
NGHS <- NGHS[order(NGHS$ID,NGHS$AGE),] 

dim(NGHS) #19398*14


IDD <- as.numeric(names(table(NGHS$ID)))
nID <-length(IDD) ; nID 


NGHS$IDD<- NA #2377 unique ID 

j<-1
 for ( ID1 in IDD ){
   if ((j-floor(j/200)*200)== 0 ) print(j) 
    NGHS$IDD[NGHS$ID == ID1] <- j
    j<-j+1
}

### select three girls with different BMI

#1: 17-24, 92-102
#2224, 30-44 , 113-116
#511: 28-33, 110-118 

NGHS[NGHS$SYSAV>140,]

NGHS[NGHS$IDD %in%  1916, c('ID', 'IDD', 'BMI', 'SYSAV') ]

NGHS.sub<- subset(NGHS, IDD %in% c(1,1916, 2225))

# can simplify code later to make it shorter , NOT JUST MAKE THE PLOT ##

postscript("fig1.3.ps", horizontal=T)

# select certain individual to make plot 3*3 , individual 3 on BMI on top
 and 3 on SBP on the bottoms ## 


############### make one plot #################
postscript("fig1.2-layout.ps", horizontal=F)

#par(mai =c(0,0, 20, 0))

layout(rbind(c(1, 1, 1, 2, 2, 2) , c(3,3,4, 4,5, 5), c(6,6, 7 , 7, 8, 8)),  respect=F )

par( mar=c(4.3, 4.7, 4, 2), cex.lab=1.5, cex.axis=1.3,pty='m')


plot(BMI~ AGE, data =  NGHS.sub100 ,  xlab = "Age (years)", ylab = expression(paste("Body Mass Index(kg/ ",m^2,")")), cex=0.4) 
mtext("A.", cex=1.4, side=3, line=0.1, font=2, at=7)

plot(SYSAV ~ AGE, data =  NGHS.sub100 ,  xlab = "Age (years)", ylab = "Systolic Blood Pressure (mm Hg)", cex=0.4) 

par(mar=c(3, 3, 2, 1), pty='m', cex.lab=1.1, cex.axis=0.9, cex=1, cex.main=1.1, font.main=1)

##---BMI --- ## 
 plot( BMI~AGE,  NGHS.sub[NGHS.sub$IDD==1,] , ylab="BMI", xlab="", main="Girl 1", ylim=c(10, 50), axes=F)
 lines(BMI~AGE,NGHS.sub[NGHS.sub$IDD==1,])
 axis(2,  mgp=c(3,0.5,0),tck=-0.04) 
 axis(1, mgp=c( 3,0.5,0),tck=-0.04)
 box()

 mtext("BMI", side=2, font=1,line=2,cex=1.1) # ylabel
 mtext("B.", cex=1.4, side=3, line=0.1, font=2, at=6.8)


 plot( BMI~AGE,  NGHS.sub[NGHS.sub$IDD==2225,]  , ylab="", xlab="", main="Girl 2", ylim=c(10, 50),axes=F)
  lines(BMI~AGE,NGHS.sub[NGHS.sub$IDD==2225,])
  axis(2,  mgp=c(3,0.5,0),tck=-0.04) 
 axis(1, mgp=c( 3,0.5,0),tck=-0.04)
 box()

 plot( BMI~AGE,  NGHS.sub[NGHS.sub$IDD==1916,]  , ylab="", xlab="", main="Girl 3", ylim=c(10, 50),axes=F)
  lines(BMI~AGE,NGHS.sub[NGHS.sub$IDD== 1916,])
  axis(2,  mgp=c(3,0.5,0),tck=-0.04) 
 axis(1, mgp=c( 3,0.5,0),tck=-0.04)
box()

##---SBP --- ## 

par(mar=c(4, 3, 0, 1), pty='m', cex.lab=1.1, cex.axis=0.9,, cex=1, cex.main=1.1, font.main=1)


 plot( SYSAV~AGE,  NGHS.sub[NGHS.sub$IDD==1,] , ylab="SBP", xlab="Age (years)", main="Girl 1", ylim=c(80, 160),axes=F)
 lines(SYSAV~AGE,NGHS.sub[NGHS.sub$IDD==1,])
mtext("C.", cex=1.4, side=3, line=0.1, font=2, at=6.8)



 mtext("SBP", side=2, font=1,line=2,cex=1.1) # ylabel
  mtext("Age (years)", cex=1, side=1, line=2)
 axis(2,  mgp=c(3,0.5,0),tck=-0.04) 
 axis(1, mgp=c( 3,0.5,0),tck=-0.04)
box()

 plot( SYSAV~AGE,  NGHS.sub[NGHS.sub$IDD==2225,]  , ylab="", xlab="", main="Girl 2", ylim=c(80, 160),axes=F)
  lines(SYSAV~AGE,NGHS.sub[NGHS.sub$IDD==2225,])
axis(2,  mgp=c(3,0.5,0),tck=-0.04) 
 axis(1, mgp=c( 3,0.5,0),tck=-0.04)
box() 
 mtext("Age (years)", cex=1, side=1, line=2)


 plot( SYSAV~AGE,  NGHS.sub[NGHS.sub$IDD==1916,]  , ylab="", xlab="Age (years)", main="Girl 3", ylim=c(80,160),axes=F)
  lines(SYSAV~AGE,NGHS.sub[NGHS.sub$IDD== 1916,])
axis(2,  mgp=c(3,0.5,0),tck=-0.04) 
 axis(1, mgp=c( 3,0.5,0),tck=-0.04)
box()
  mtext("Age (years)", cex=1, side=1, line=2)

dev.off()


### Figure 3: enrichd ##
## BDI before and after medication curves ##
## check BDI score, change, before and after med ##
 

 BDIdata<-read.table(file="BDImed.txt", head=T)
 dim(BDIdata); #[1] 1446    6
 
length(unique(BDIdata$ID)) #91
 BDIdata[1:40,];
 BDIdata <- data.frame(BDIdata, BDIdata[,3]*BDIdata[,2])
 names(BDIdata)<- c("ID","Tij","Si","Deltaij","Rij" , "Y", "SiTij")
  # Tij- month since study 
  # Si -medication time(month)
  # Rij- month after medication
  #  Y BDI score

summary(BDIdata$Si)*365.25/12
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
  0.00000   0.00000  20.00048  44.01262  78.98531 172.00231

summary(BDIdata$Si[BDIdata$Si>0])*365.25/12

  Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
  7.00063  27.99946  76.00244  73.08044 111.00556 172.00231


BDIdata$visitday<- BDIdata$Tij*365.25/12


BDIdata$visitday2<-  (BDIdata$Tij- BDIdata$Si)*365.25/12




postscript("fig1.3.ps", horizontal=T)

par(mar=c(4.5, 4.5, 3, 1),cex=0.5,mfrow=c(1,2), pty='m', cex.lab=1.6, cex.axis=1.2)

plot(BDIdata$visitday,  BDIdata$Y,  xlab = "Days on study", ylab = "BDI score", ylim=c(0, 60), cex=0.5)

  for (i in unique(BDIdata$ID)) { lines(Y ~ visitday, data = BDIdata[BDIdata$ID  == i, ], lty = 1)}

  # align before and after medication

  
plot(BDIdata$visitday2,  BDIdata$Y,  xlab = "Days since medication", ylab = "BDI score", ylim=c(0, 60), cex=0.5)

for (i in unique(BDIdata$ID)) { lines(Y ~ visitday2, data = BDIdata[BDIdata$ID  == i, ], lty = 1)}
abline(v=0,lty=2)

dev.off()

 
##  Figure 4:  Cytokines: short curve but 6-9 variables :Lykocyte and cyto ##
 # lowess or spline fit , on the cytokine change ##
 # Plot log m then add lowess line. 6 plots #

# I will plot 3 blood counts (gran, Lym and Mono) and three cytokine GCSF, IL-15 and MCP-1.
# took at all zeros, non-detectable for easy fitting.  
# suggest a cytokine storm, conditioning, and recovering.

 # try to use ggplot or lowess to generate 2*3 or 3*3 plot ##
 #lowess: just add a thick blue line 
  # study ggplot (or ask Dihua later, how to change to white, grid , add smooth line and CI ##
   # describe the data 
   
   
load(file="Allcytokine.RData")

dim(Pat) # 545* 40 var

names(Pat)

[1] "ID"       "DAYS"     "CRP"      "LYM"      "MON"      "PMN"      "Temp"     "GCSF"     "IL15"     "TNFa"     "IL6"     
[12] "IL7"      "IL1"      "IL2"      "IFNg"     "MCP1"     "SCF"      "RANTES"   "HGF"      "SCGFb"    "MIP1b"    "IP10"    
[23] "IL18"     "GMCSF"    "PDGFbb"   "IL9"      "Eotaxin"  "IFNa2"    "IL1ra"    "LIF"      "FGFbasic" "IL8"      "IL17"    
[34] "IL4"      "IL13"     "VEGF"     "IL12p70"  "IL5"      "IL10"     "MIP1a"   


Lst <- unique(Pat[,1])#names(table(Pat[,1]))

names(Pat)

#NEED DAYS, PMN, LMN, MON, GCSF,IL15, MCP1 # all reached at peak during cytopenia

## PMN 
Varno<- 6 ; names(Pat)[Varno]

PMNall <- data.frame(Pat[c(1,2,Varno)])

PMNall<- PMNall[!is.na(PMNall[,3]) & PMNall$DAYS  <=35 & PMNall$DAYS >=-8 ,] #333

PMNall[ PMNall[,3]==0 ,3]<- NA

PMNall<- PMNall[!is.na(PMNall[,3]) & PMNall$DAYS  <=35 & PMNall$DAYS >=-8 ,] #333


 PMNall$log10.X  <- log10(PMNall[,3])
 summary( PMNall$log10.X)

## LYM ##


Varno<- 4 ; names(Pat)[Varno]

LYMall <- data.frame(Pat[c(1,2,Varno)])

LYMall<- LYMall[!is.na(LYMall[,3]) & LYMall$DAYS  <=35 & LYMall$DAYS >=-8 ,] #333

LYMall[ LYMall[,3]==0 ,3]<- NA

LYMall<- LYMall[!is.na(LYMall[,3]) & LYMall$DAYS  <=35 & LYMall$DAYS >=-8 ,] #333


 LYMall$log10.X  <- log10(LYMall[,3])
 summary( LYMall$log10.X)

##  MON ##


Varno<- 5 ; names(Pat)[Varno]

MONall <- data.frame(Pat[c(1,2,Varno)])

MONall<- MONall[!is.na(MONall[,3]) & MONall$DAYS  <=35 & MONall$DAYS >=-8 ,] #333

MONall[ MONall[,3]==0 ,3]<- NA

MONall<- MONall[!is.na(MONall[,3]) & MONall$DAYS  <=35 & MONall$DAYS >=-8 ,] #333


 MONall$log10.X  <- log10(MONall[,3])
 summary( MONall$log10.X)

####  GCSF ##

Varno<- 8 ; names(Pat)[Varno]

GCSFall <- data.frame(Pat[c(1,2,Varno)])

GCSFall<- GCSFall[!is.na(GCSFall[,3]) & GCSFall$DAYS  <=35 & GCSFall$DAYS >=-8 ,] #333

GCSFall[ GCSFall[,3]<= 0.12 ,3] <- NA

GCSFall<- GCSFall[!is.na(GCSFall[,3]) & GCSFall$DAYS  <=35 & GCSFall$DAYS >=-8 ,] #333


 GCSFall$log10.X  <- log10(GCSFall[,3])
 summary( GCSFall$log10.X)


####  IL-15 ##

Varno<- 9 ; names(Pat)[Varno]

IL15all <- data.frame(Pat[c(1,2,Varno)])

IL15all<- IL15all[!is.na(IL15all[,3]) & IL15all$DAYS  <=35 & IL15all$DAYS >=-8 ,] #333

IL15all[ IL15all[,3]<=0.84  ,3 ]   <- NA

IL15all<- IL15all[!is.na(IL15all[,3]) & IL15all$DAYS  <=35 & IL15all$DAYS >=-8 ,] #333


 IL15all$log10.X  <- log10(IL15all[,3])
 summary( IL15all$log10.X)

####  MCP1 ##

Varno<- 16 ; names(Pat)[Varno]

MCP1all <- data.frame(Pat[c(1,2,Varno)])

MCP1all<- MCP1all[!is.na(MCP1all[,3]) & MCP1all$DAYS  <=35 & MCP1all$DAYS >=-8 ,] #333

MCP1all[ MCP1all[,3]<=0.41 ,3]<- NA

MCP1all[ MCP1all[,3]> 10000 ,3]<- NA


MCP1all<- MCP1all[!is.na(MCP1all[,3]) & MCP1all$DAYS  <=35 & MCP1all$DAYS >=-8 ,] #333


 MCP1all$log10.X  <- log10(MCP1all[,3])
 summary( MCP1all$log10.X)


## ?loess()
## first loess return the model or function call from loess
## then use predict loess to get new prediction
##  then lines to plot.



postscript("fig1.4a.ps", horizontal=T)# case example

 par(mar=c(4.5, 5.5, 3,1),cex=2,mfrow=c(1,1), pty='m', cex.lab=1.6, cex.axis=1.4, cex.main=1.6, cex=1.1)

 plot(log10.X~DAYS  ,data=PMNall[PMNall$ID=='VALJB',], pch=16, ylim=c(-4,3), xlab="Days post-transplant", ylab="", xlim=c(-8,30) ,
  axes=F, main="An example of one subject mutiple series" )
 lines(log10.X~DAYS,data=PMNall[PMNall$ID=='VALJB',])


 axis(1, at=seq(-7, 28, by=7),cex.axis=1.6) 
 box(, lwd=1)

 points(log10.X~DAYS  ,data=LYMall[LYMall$ID=='VALJB',], pch=16, col=2)
 lines(log10.X~DAYS,data=LYMall[LYMall$ID=='VALJB',],col=2)


 points(log10.X~DAYS  ,data=MONall[MONall$ID=='VALJB',], pch=16, col=3)
 lines(log10.X~DAYS,data=MONall[MONall$ID=='VALJB',],col=3)


 points(log10.X~DAYS  ,data=GCSFall[GCSFall$ID=='VALJB',], pch=16, col=4)
 lines(log10.X~DAYS,data=GCSFall[GCSFall$ID=='VALJB',],col=4)


 points(log10.X~DAYS  ,data=IL15all[IL15all$ID=='VALJB',], pch=16, col='orange')
 lines(log10.X~DAYS,data=IL15all[IL15all$ID=='VALJB',],col='orange')


 points(log10.X~DAYS  ,data=MCP1all[MCP1all$ID=='VALJB',], pch=16, col='brown')
 lines(log10.X~DAYS,data=MCP1all[MCP1all$ID=='VALJB',],col='brown')

 legend('bottomright', legend= c("PMN", "LYM", "MON", "GCSF", "IL-15", "MCP-1" ),  pch=16, lty=1,
   col=c("black","red","green","blue","orange","brown") ,cex=1.2)


postscript("fig1.4.ps", horizontal=T)


par(mar=c(4.5, 5.5, 3,1),cex=2,mfrow=c(2,3), pty='m', cex.lab=1.6, cex.axis=1.4, cex.main=1.6)


#-----------PMN-----------

plot(PMNall$DAYS  , PMNall$log10.X,   main='Granulocyte',  col='black', pch=1, cex=1,
  xlab='', ylab='',  xlim=c(-8,30) , ylim=c(-3,1.5), axes=F)


axis(2,at= -3:1,  labels= c(1, 10,100,1000,10000), cex.axis=1.6 , las=2 )

axis(1, at=seq(-7, 28, by=7),cex.axis=1.6) 
box(, lwd=1)

fit<- loess(log10.X~ DAYS, data=PMNall )
lines(seq(-8, 35, 1), predict(fit, data.frame(DAYS  = seq(-8, 35, 1))),col=4, lwd=2) 

#-----------LYM  -------------

plot(LYMall$DAYS  , LYMall$log10.X,   main='Lymphocyte',  col='black', pch=1, cex=1,
  xlab='', ylab='',  xlim=c(-8,30) , ylim=c(-3,0.5), axes=F)


axis(2,at= -3:0,  labels= c(1, 10,100,1000), cex.axis=1.6 , las=2 )

axis(1, at=seq(-7, 28, by=7),cex.axis=1.6) 
box(, lwd=1)


fit<- loess(log10.X~ DAYS, data=LYMall)
lines(seq(-8, 35, 1), predict(fit, data.frame(DAYS  = seq(-8, 35, 1))),col=4, lwd=2) 


#-----------MON -------------

plot(MONall$DAYS  , MONall$log10.X,   main='Monocyte',  col='black', pch=1, cex=1,
  xlab='', ylab='',  xlim=c(-8,30) , ylim=c(-3,0.5), axes=F)


axis(2,at= -3:0,  labels= c(1, 10,100,1000), cex.axis=1.6 , las=2 )

axis(1, at=seq(-7, 28, by=7),cex.axis=1.6) 
box(, lwd=1)

fit<- loess(log10.X~ DAYS, data=MONall)
lines(seq(-8, 35, 1), predict(fit, data.frame(DAYS  = seq(-8, 35, 1))),col=4, lwd=2) 


## now three cytokines ###


#-----------GCSF -------------

plot(GCSFall$DAYS  , GCSFall$log10.X,   main='G-CSF',  col='black', pch=1, cex=1,
  xlab='Days post-transplant', ylab='',  xlim=c(-8,30) , ylim=c(-0.5,4), axes=F)


axis(2,at= 0:4,  labels= c(1, 10,100,1000,10000), cex.axis=1.6 , las=2 )

axis(1, at=seq(-7, 28, by=7),cex.axis=1.6) 
box(, lwd=1)


fit<- loess(log10.X~ DAYS, data=GCSFall)
lines(seq(-8, 35, 1), predict(fit, data.frame(DAYS  = seq(-8, 35, 1))),col=4, lwd=2) 


#-----------IL15 -------------

plot(IL15all$DAYS  , IL15all$log10.X,   main='IL-15',  col='black', pch=1, cex=1,
  xlab='Days post-transplant', ylab='',  xlim=c(-8,30) , ylim=c(0,3), axes=F)


axis(2,at= 0:3,  labels= c(1, 10,100,1000), cex.axis=1.6 , las=2 )

axis(1, at=seq(-7, 28, by=7),cex.axis=1.6) 
box(, lwd=1)


fit<- loess(log10.X~ DAYS, data=IL15all)
lines(seq(-8, 35, 1), predict(fit, data.frame(DAYS  = seq(-8, 35, 1))),col=4, lwd=2) 


#-----------MCP1 -------------

plot(MCP1all$DAYS  , MCP1all$log10.X,   main="MCP-1",  col='black', pch=1, cex=1,
  xlab='Days post-transplant', ylab='',  xlim=c(-8,30) , ylim=c(0,4), axes=F)


axis(2,at= 0:4,  labels= c(1, 10,100,1000,10000), cex.axis=1.6 , las=2 )

axis(1, at=seq(-7, 28, by=7),cex.axis=1.6) 
box(, lwd=1)



fit<- loess(log10.X~ DAYS, data=MCP1all)
lines(seq(-8, 35, 1), predict(fit, data.frame(DAYS  = seq(-8, 35, 1))),col=4, lwd=2) 


dev.off()




#----------------------------------------------------------#
### how to how to use ggplot to add a smooth spline and CI#
library(ggplot2)
library(grid)
library(gridExtra)



postscript("fig1.4b.ps", horizontal=T)

pdf(file="fig1.4b.pdf")
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,3 )))  

p1<- ggplot(PMNall, aes(x=DAYS  , y=log10.X  )) + geom_point() +  geom_smooth( stat = "smooth", size=1)+ xlab("")+ylab("")  + ggtitle((title = 'Granulocyte'))
p2<- ggplot(LYMall, aes(x=DAYS  , y=log10.X  )) + geom_point() +  geom_smooth(size=1)+ xlab("")+ylab("")  + ggtitle((title = "Lymphocyte"))
p3<- ggplot(MONall, aes(x=DAYS  , y=log10.X  )) + geom_point() +  geom_smooth(size=1)+ xlab("")+ylab("")  + ggtitle((title = "Monocyte"))

p4<- ggplot(GCSFall, aes(x=DAYS  , y=log10.X  )) + geom_point() +  geom_smooth(size=1)+ xlab("")+ylab("")  + ggtitle((title = 'G-CSF'))
p5<- ggplot(IL15all, aes(x=DAYS  , y=log10.X  )) + geom_point() +  geom_smooth(size=1)+ xlab("")+ylab("")  + ggtitle((title = "IL-15"))
p6<- ggplot(MCP1all, aes(x=DAYS  , y=log10.X  )) + geom_point() +  geom_smooth(size=1)+ xlab("")+ylab("")  + ggtitle((title = "MCP-1"))


print(p1+theme_bw(),  vp=viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p2+theme_bw(),  vp=viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p3+theme_bw(),  vp=viewport(layout.pos.row = 1, layout.pos.col = 3))

print(p4+theme_bw(),  vp=viewport(layout.pos.row = 2, layout.pos.col = 1))
print(p5+theme_bw(),  vp=viewport(layout.pos.row = 2, layout.pos.col = 2))
print(p6+theme_bw(),  vp=viewport(layout.pos.row = 2, layout.pos.col = 3))
dev.off()



 #------------------------------------------
 # Another data source :later check Gail with IPF lung function data , can I request it from biolincc. ##
 # check how to request
 #panther IPF NAC vs. placebo should be good example, FVC primary 0,15,30,45,60, check treatment difference
  #also check STEP IPF data , what is primary?  






























































