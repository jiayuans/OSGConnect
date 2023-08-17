#!/usr/bin/env Rscript

dirg <- "C:/UCHealth/RA/Project/EPIC-CF/Code/jags/"

library(mvtnorm)
library(rjags)
library(R2jags)
library(runjags)
library(Hmisc)
library(mcmcplots)
library(tidyverse)

#############################################################################
##read in data#################
##################Reading in main dataset##########################
#################################################################
dat <- read.csv(file=paste(dirg,"/Data-PA-cohort.csv",sep=""))
dat <- dat[,-1]
head(dat, n=10)

N <- length(dat$cffidno)
M <- length(unique(dat$cffidno))
length(dat$cffidno)
dat.uniq <- dat[!duplicated(dat$cffidno), ]
dat.uniq$num <- c(1:M)
dat.uniq1 <-dat.uniq[,c(1,12)]
dat1 <- merge(dat, dat.uniq1, by=c('cffidno'), all=TRUE)

count <- dat %>% count(cffidno)
 
##########################Assigning unique number to each subject##########################
dat1 <- dat1[order(dat$cffidno) , ]

#################Readingin data for VisitAge matrix#############################

X <- dat1$VisitAge 

#################Readingin data for Pa+ matrix#############################

Y <- dat1$cltpa
num <- dat1$num

#################input variables for simulation#####################
#################################################################
####### model finding change point, fixed change points
##---------------------------------------------------------------         
#################################################################
#### 

modelrancp <- "model { 
  for(i in 1:N){ 
  ### PA model
        Y[i] ~ dbin(p2[i],1)
        logit(p2[i]) <-  c0 + c[1] * (X[i]-cp[1]) + c[2] * (X[i]-cp[1]) * step(X[i]-cp[1]) + c[3] * (X[i]-cp[2]) * step(X[i]-cp[2]) + u[num[i]]
  }
  for(j in 1:M){ 
        u[j] ~ dnorm(0,tau)
  }
  c0 ~ dnorm(0,0.0001); 
	for (k in 1:3){
	      c[k] ~ dnorm(0,0.0001);		
	}
	tau <- pow(sig,-2); ## precision
	sig ~ dunif(0,100);  ## variance
	for (k in 1:2){
	      cp.temp[k] ~ dunif(min,max);		
	}
	cp[1:2]<-sort(cp.temp)
}"

####Observed DATA 
data <- dump.format(list(X=X, Y=Y, M=M, N=N, min=min(X), max=max(X), num=num)) 
###initial Values
inits1 <- dump.format(list(c0=-4.2, c=c(0,0,0),
                           sig=1.4, cp.temp=c(8,14),
                           .RNG.name="base::Super-Duper", .RNG.seed=1))
inits2 <- dump.format(list(c0=-4.3, c=c(0,0,0)+0.01, 
                           sig=1.5, cp.temp=c(8,14)+0.1,
                           .RNG.name="base::Super-Duper", .RNG.seed=2))
#### Run the model and produce plots
res1 <- run.jags(model=modelrancp, burnin=2000, sample=2000, 
                 monitor=c("c0", "c", "sig", "cp", "cp.temp", "u"), 
                 data=data, n.chains=2, inits=c(inits1,inits2), thin=10, plots = FALSE)

sum <- summary(res1)
sum
#Lower95       Median      Upper95          Mean         SD Mode       MCerr MC%ofSD SSeff        AC.100      psrf
#c0         -3.869522497 -3.703125935 -3.552665184 -3.702334e+00 0.08069332   NA 0.005946015     7.4   184  2.845252e-01 1.0190303
#c[1]       -0.279467592 -0.147086798 -0.068675640 -1.555021e-01 0.05573610   NA 0.014955607    26.8    14  8.784821e-01 1.1086884
#c[2]        0.208923520  0.294680117  0.420187574  3.029601e-01 0.05552456   NA 0.014035683    25.3    16  8.692032e-01 1.0917959
#c[3]        0.107313094  0.158401453  0.213194774  1.605906e-01 0.03410654   NA 0.001238686     3.6   758  9.157736e-02 1.0180476
#sig         1.511698956  1.586029035  1.673861881  1.586990e+00 0.04187450   NA 0.001113267     2.7  1415 -3.896497e-02 1.0029415
#cp[1]       3.488181358  4.231414390  5.159705571  4.272391e+00 0.43471904   NA 0.070870796    16.3    38  5.135209e-01 1.0495947
#cp[2]      13.261157968 14.277713579 15.493423509  1.432591e+01 0.60328626   NA 0.022405221     3.7   725  5.877915e-02 1.0268652
