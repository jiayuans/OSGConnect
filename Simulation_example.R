#!/usr/bin/env Rscript 

##============================
library(coda)
library(rjags)
library(runjags)
library(tidyverse)

##################################################################
##    Functions to calculate the WAIC and Read data
## 
##################################################################
dat <- read.csv(file="Data-PA-cohort.csv")
dat <- dat[,-1]
head(dat, n=10)

M <- length(unique(dat$cffidno))
length(dat$cffidno)

count <- dat %>% count(cffidno)

dat1 <-aggregate(dat$VisitAge, by=list(dat$cffidno),
                 FUN=max, na.rm=TRUE)
names(dat1) <- c('cffidno','age.max')
head(dat1)
nrow(dat1)

dat2 <-aggregate(dat$age.min, by=list(dat$cffidno),
                 FUN=max, na.rm=TRUE)
names(dat2) <- c('cffidno','age.min')
head(dat2)
nrow(dat2)

first.t<-dat2$age.min
last.t<-dat1$age.max
first.tt<-first.t[1:100]
last.tt<-last.t[1:100]

####time of first visit and last visit#######
N<-length(last.tt)
###set number of iterations#################################
I=11

###############set true values#########################################
b0=2
b1=0
b2=-4
Verror=1
cp.true=15
s=23###starting seed####
#############################################################

B1.mean<-rep(NA,I-1)
B2.mean<-rep(NA,I-1)
c0.mean<-rep(NA,I-1)
c1.mean<-rep(NA,I-1)
c2.mean<-rep(NA,I-1)
cp.mean<-rep(NA,I-1)
pi1.mean<-rep(NA,I-1)
pi2.mean<-rep(NA,I-1)
z.mean<-rep(NA,I-1)
tau.mean<-rep(NA,I-1)
u.mean<-rep(NA,I-1)
cp.mu.mean<-rep(NA,I-1)
cp.tau.mean<-rep(NA,I-1)

#participant ID
ID<-rep(1:N)
length(ID)

#######################################################
for (r in 2:I){
  set.seed(s+100*(r-1))
  cpset<-rnorm(N,mean = 15, sd = 3)
  t<-round(first.tt)
  tt<-round(last.tt)
  kt<-(tt-t)*4
  kk=max(kt)
  
  U<-runif(N)
  a<-rnorm(N,0,1)
  a2<-rnorm(N,0,1)
  
  I1<-matrix(NA, nrow=N, ncol=kk, byrow=TRUE)
  I2<-matrix(NA, nrow=N, ncol=kk, byrow=TRUE)
  p2<-matrix(NA, nrow=N, ncol=kk, byrow=TRUE)
  Y<-matrix(NA, nrow=N, ncol=kk, byrow=TRUE)
  X<-matrix(NA, nrow=N, ncol=kk, byrow=TRUE)
  z<-rep(NA,N)
  for (i in 1:N){
    X[i,1:kt[i]]<-c(seq(t[i],tt[i]-0.25,0.25))
  }
  #fix(X)
  
  for (i in 1:N){
    for (j in 1:kt[i]){
      if (U[i]<0.80) {
        I1[i,j]<-ifelse(X[i,j]<=cpset[i],1,0)
        I2[i,j]<-ifelse(X[i,j]>=cpset[i],1,0)
        p2[i,j]=exp(b0+b1*(X[i,j]-cpset[i])*I1[i,j]+b2*(X[i,j]-cpset[i])*I2[i,j]+a[i])/(1+exp(b0+b1*(X[i,j]-cpset[i])*I1[i,j]+b2*(X[i,j]-cpset[i])*I2[i,j]+a[i]))
        Y[i,j]=rbinom(Verror, 1, p2[i,j])
        z[i]=1
      } else {
        p2[i,j]=exp(b0+b1*(X[i,j]-10)+a[i])/(1+exp(b0+b1*(X[i,j]-10)+a[i]))
        Y[i,j]=rbinom(Verror, 1, p2[i,j])
        z[i]=2
      }
    }
  }
  
  alpha = c(1,1)
  
  ############Model in the JAGS format#####################
  ############One fixed mixture CP#####################  
  modelrancp <- "model { 
    for(i in 1:N){ 
         for(j in 1:kt[i]){
         ### PA model
            Y[i,j] ~ dbin(p2[i,j],1)
            logit(p2[i,j]) <- logit(p[i,j,z[i]])
            logit(p[i,j,1]) <-  c0 + c[1] * (X[i,j]-cp[i]) + c[2] * (X[i,j]-cp[i]) * (2*step(X[i,j]-cp[i])-1) + u[i]
            logit(p[i,j,2]) <-  c0 + B1*(X[i,j]-10) + u[i]
         } 
         z[i]~dcat(pi[1:2])
         u[i] ~ dnorm(0,tau)
         cp[i] ~ dnorm(cp.mu,cp.tau)
    }
    pi[1:2] ~ ddirch(alpha[])
    c0 ~ dnorm(0,0.0001); 
  	for (k in 1:2){
  	     c[k] ~ dnorm(0,0.0001)		
  	}
  	tau ~ dgamma(0.001,0.001)
    cp.mu ~ dnorm(0,0.001)
	  cp.tau ~ dgamma(0.001,0.001)
    B1<- c[1]-c[2]
    B2<- c[1]+c[2]
}"
  
  ####Observed DATA 
  data <- dump.format(list(X=X, Y=Y, N=N, alpha=alpha, kt=kt)) 
  ###initial Values
  inits1 <- dump.format(list(c0=2, c=c(-2,-2), cp.mu=15, cp.tau=1/9, tau=1, pi=c(0.79,0.21),
                             .RNG.name="base::Super-Duper", .RNG.seed=1))
  inits2 <- dump.format(list(c0=2.1, c=c(-2,-2)+0.01, cp.mu=15.1, cp.tau=1/9, tau=1, pi=c(0.8,0.2),
                             .RNG.name="base::Super-Duper", .RNG.seed=2))
  
  #### Run the model and produce plots
  res1 <- run.jags(model=modelrancp, burnin=8000, sample=5000, 
                   monitor=c("B1", "B2", "c0", "c", "cp", "pi", "z", "tau","u","cp.mu","cp.tau"), 
                   data=data, n.chains=2, inits=c(inits1,inits2), thin=10, plots = FALSE)
  
  summary <- summary(res1)
  
  B1.mean[r]<-summary[1,4] 
  B2.mean[r]<-summary[2,4] 
  c0.mean[r]<-summary[3,4] 
  c1.mean[r]<-summary[4,4] 
  c2.mean[r]<-summary[5,4] 
  cp.mean[r]<-mean(summary[6:105,4]) 
  pi1.mean[r]<-summary[106,4] 
  pi2.mean[r]<-summary[107,4] 
  z.mean[r]<-summary[108:207,4] 
  tau.mean[r]<-summary[208,4]
  u.mean[r]<-mean(summary[209:308,4])
  cp.mu.mean[r]<-summary[309,4] 
  cp.tau.mean[r]<-summary[310,4]
}

Sim.results=cbind(B1.mean,B2.mean,c0.mean,c1.mean,c2.mean,cp.mean,pi1.mean,pi2.mean,z.mean,tau.mean,u.mean,cp.mu.mean,cp.tau.mean)
Sim.results=Sim.results[-1,]
round(colMeans(Sim.results),2)

#write.csv(summary,"C:/UCHealth/RA/Project/EPIC-CF/Analysis_Jiayuan/Result/Simulation_mixture_1random_cp_summary.csv")
#write.csv(Sim.results,"C:/UCHealth/RA/Project/EPIC-CF/Analysis_Jiayuan/Result/Simulation_mixture_1random_cp.csv")


# B1.mean  B2.mean  c0.mean  c1.mean  c2.mean  cp.mean pi1.mean pi2.mean   z.mean tau.mean   u.mean 
#0.00    -4.07     2.02    -2.04    -2.04    15.00     0.82     0.18     1.10     1.07     0.00 

## true values:
#b1=0
#b2=-4
#c0=2
#c1=-2
#c2=-2
#cp1.true=15
#tau=1

