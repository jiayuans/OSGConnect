#!/usr/bin/env Rscript

library(runjags)


first.tt <- c(12.06,10.13,9.95,12.74,9.03,10.48,10.28,10.28,10.39,9.95,13.97,10.21,11.94,10.39,10.36,10.79,9.81,9.87,10.51,11.35,11.64,12.16,12.63,11.69,13.00,11.72,12.49,12.94,11.52,12.34,12.38,11.56,12.23,12.31,11.60,11.54,11.02,11.82,10.56,12.99,10.75,13.91,11.63,10.68,12.14,11.56,10.35,11.25,11.63,11.52,12.91,0.51,0.45,3.68,0.44,0.82,0.44,0.69,5.50,0.15,1.42,0.31,0.77,0.17,0.35,1.10,2.11,1.07,0.07,12.52,1.03,1.67,0.28,0.38,1.18,1.39,0.53,0.59,0.76,2.12,1.44,5.00,0.92,2.56,3.29,0.20,2.49,10.64,1.40,3.95,2.04,6.23,0.47,1.22,10.23,2.60,1.35,3.41,1.11,0.23,10.34,1.57,5.53,4.66,0.88,8.02,1.10,1.67,1.99,0.28,1.72,1.55,0.64,7.37,0.44,2.47,1.97,1.02,0.58,6.37,0.09,2.29,3.05,3.00,0.27,4.12,0.65,3.18,1.28,2.03,0.22,1.30,6.90,0.27,3.08,3.31,0.21,2.97,4.49,4.62,2.20,2.73,1.69,1.25,0.15,0.33,3.65,0.36,1.41,1.46,0.29,5.89,0.81,3.68,2.39,0.12,4.45,0.11,0.87,1.00,9.09,1.51,8.94,1.29,2.60,2.41,2.41,1.23,3.01,2.54,1.88,0.85,1.39,3.55,1.30,1.71,3.41,3.92,4.82,10.96,1.83,2.53,11.85,1.21,6.80,1.30,1.41,0.49,1.00,0.78,1.33,0.34,0.93,7.97,1.86,1.45,0.87,9.29,1.22,3.74)
last.tt <- c(19.11,18.44,18.39,21.45,17.96,18.34,18.68,18.46,18.82,18.52,18.86,18.61,19.03,19.23,13.85,18.73,13.34,18.42,14.10,17.85,20.02,20.68,20.56,20.37,21.35,15.36,20.28
             ,18.46,15.27,19.59,20.51,18.84,20.45,16.15,19.55,19.72,19.64,15.16,18.96,19.04,19.30,19.86,19.22,19.14,14.48,19.67,19.24,19.95,19.43,19.47,15.11,8.56,8.20,11.79
             ,8.16,8.33,8.30,8.30,11.38,8.27,9.46,8.31,9.02,8.30,8.36,9.00,10.20,8.57,8.42,19.06,8.22,8.65,8.45,8.33,9.41,8.28,8.36,8.35,8.39,8.93,8.28
             ,12.64,7.99,8.33,10.23,7.98,4.30,12.98,8.28,10.83,8.04,11.75,8.30,3.47,12.49,9.52,8.31,11.32,9.19,8.18,18.16,9.47,11.37,8.59,8.57,14.77,9.14,9.35,8.80,8.22,8.58,8.53,3.69,14.49,6.80,5.53,4.08,7.80,3.58,14.37,8.02,8.20,6.25,8.89,8.13,11.75,8.34,4.39,8.09,9.70,7.93,8.26,14.74,7.95,10.74,8.05,7.96,8.70,8.21,8.10,10.41,10.64,3.94,8.40,8.44,3.27,11.83,3.96,4.40,8.38,8.50,13.44,3.91,10.43,10.49,8.53,11.17,8.54,9.49,4.44,12.49,9.44,16.88,9.40,10.25,9.46,9.50,4.66,11.16,9.56,10.02,9.36,9.80,9.59,9.78,4.82,11.43,12.44,12.98,14.45,9.52,9.45,15.35,9.56,14.76,9.40,9.18,9.10,8.96,4.23,9.18,8.91,8.96,16.57,9.77,9.56,9.54,17.15,9.36,11.43)

####time of first visit and last visit#######
N<-length(last.tt)
###set number of iterations#################################
I=21

###############set true values#########################################
c0=-3.7
c1=-0.12
c2=0.27
c3=0.16
Verror=1
cp1.true=4.5
cp2.true=14.4
s=23###starting seed####
#############################################################

B1.mean<-rep(NA,I-1)
B2.mean<-rep(NA,I-1)
B3.mean<-rep(NA,I-1)
c0.mean<-rep(NA,I-1)
c1.mean<-rep(NA,I-1)
c2.mean<-rep(NA,I-1)
c3.mean<-rep(NA,I-1)
cp1.mean<-rep(NA,I-1)
cp2.mean<-rep(NA,I-1)
tau.mean<-rep(NA,I-1)
u.mean<-rep(NA,I-1)
cp1.mu<-rep(NA,I-1)
cp1.tau<-rep(NA,I-1)
cp2.temp<-rep(NA,I-1)


#participant ID
ID<-rep(1:N)
length(ID)


#######################################################
for (r in 2:I){
  set.seed(s+100*(r-1))
  t<-round(first.tt)
  tt<-round(last.tt)
  kt<-(tt-t)*4
  kk=max(kt)
  
  U<-runif(N)
  a<-rnorm(N,0,0.4)
  
  I1<-matrix(NA, nrow=N, ncol=kk, byrow=TRUE)
  I2<-matrix(NA, nrow=N, ncol=kk, byrow=TRUE)
  p2<-matrix(NA, nrow=N, ncol=kk, byrow=TRUE)
  Y<-matrix(NA, nrow=N, ncol=kk, byrow=TRUE)
  X<-matrix(NA, nrow=N, ncol=kk, byrow=TRUE)
  for (i in 1:N){
    X[i,1:kt[i]]<-c(seq(t[i],tt[i]-0.25,0.25))
  }
  
  for (i in 1:N){
    for (j in 1:kt[i]){
      I1[i,j]<-ifelse(X[i,j]< cp1.true,-1,1)
      I2[i,j]<-ifelse(X[i,j]< cp2.true,-1,1)
      p2[i,j]=exp(c0+c1*(X[i,j]-cp1.true)+c2*(X[i,j]-cp1.true)*I1[i,j]+c3*(X[i,j]-cp2.true)*I2[i,j]+a[i])/(1+exp(c0+c1*(X[i,j]-cp1.true)+c2*(X[i,j]-cp1.true)*I1[i,j]+c3*(X[i,j]-cp2.true)*I2[i,j]+a[i]))
      Y[i,j]=rbinom(Verror, 1, p2[i,j])
    }
  }
  
  ############Model in the JAGS format#####################
  ############Two fixed CP#####################  
  modelrancp <- "model { 
    for(i in 1:N){ 
         for(j in 1:kt[i]){
         ### PA model
             Y[i,j] ~ dbin(p2[i,j],1)
             logit(p2[i,j]) <- c0 + c[1] * (X[i,j]-cp1) + c[2] * (X[i,j]-cp1) * (2*step(X[i,j]-cp1)-1) + c[3] * (X[i,j]-cp2) * (2*step(X[i,j]-cp2)-1) + u[i]
         } 
         u[i] ~ dnorm(0,tau)
     }
  c0 ~ dnorm(0,0.0001) 
	for (k in 1:3){
	     c[k] ~ dnorm(0,0.0001);		
	}
	tau ~ dgamma(0.001,0.001)
  cp1 ~ dnorm(cp1.mu,cp1.tau)	
	cp2.temp ~ dunif(0,max)
	cp2 <- cp1 + cp2.temp
	cp1.mu ~ dnorm(0,0.001)
	cp1.tau ~ dgamma(0.001,0.001)
  B1<-c[1]-c[2]-c[3]
  B2<-c[1]+c[2]-c[3]
  B3<-c[1]+c[2]+c[3]
}"
  
  ####Observed DATA 
  data <- dump.format(list(X=X, Y=Y, N=N, kt=kt, max=max(tt))) 
  ###initial Values
  inits1 <- dump.format(list(c0=-3.7, c=c(-0.1,0.3,0.2),  cp1=4, cp2.temp=8, tau=1,
                             .RNG.name="base::Super-Duper", .RNG.seed=1))
  inits2 <- dump.format(list(c0=-3.6, c=c(-0.1,0.3,0.2)+0.01, cp1=4.1, cp2.temp=8.1, tau=1,
                             .RNG.name="base::Super-Duper", .RNG.seed=2))
  
  #### Run the model and produce plots
  res1 <- run.jags(model=modelrancp, burnin=12000, sample=12000, 
                   monitor=c("B1", "B2","B3", "c0", "c", "cp1", "cp2", "tau","u","cp1.mu","cp1.tau", "cp2.temp"), 
                   data=data, n.chains=2, inits=c(inits1,inits2), thin=10, plots = FALSE)
  
  summary <- summary(res1)
  
  B1.mean[r]<-summary[1,4] 
  B2.mean[r]<-summary[2,4] 
  B3.mean[r]<-summary[3,4] 
  c0.mean[r]<-summary[4,4] 
  c1.mean[r]<-summary[5,4] 
  c2.mean[r]<-summary[6,4] 
  c3.mean[r]<-summary[7,4]
  cp1.mean[r]<-summary[8,4] 
  cp2.mean[r]<-summary[9,4] 
  tau.mean[r]<-summary[10,4]
  u.mean[r]<-mean(summary[11:110,4])
  cp1.mu[r]<-summary[111,4]
  cp1.tau[r]<-summary[112,4]
  cp2.temp[r]<-summary[113,4]
  
}

Sim.results=cbind(B1.mean,B2.mean,B3.mean,c0.mean,c1.mean,c2.mean,c3.mean,cp1.mean,cp2.mean,tau.mean,u.mean)
Sim.results=Sim.results[-1,]
