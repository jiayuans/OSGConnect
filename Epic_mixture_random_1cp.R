#############################################################################
##################Mixture Model #################
##################read in data#################
##################Reading in main dataset##########################
#################################################################
library(coda)
library(rjags)
library(runjags)
library(tidyverse)

dat0 <- read.csv(file="Data-PA-cohort.csv")
set.seed(2256)
dat<- dat0[dat0$ru < 0.2, -1] #subset of the data
length(unique(dat0$cffidno)) #1734
length(unique(dat$cffidno)) #323
head(dat, n=10)

M <- length(unique(dat$cffidno))
length(dat$cffidno)

count <- dat %>% count(cffidno)
max(count$n)

##########################Assigning unique number to each subject##########################
dat1 <- dat[order(dat$cffidno) , ]
dat2 <- dat1 %>% group_by(cffidno) %>% mutate(time = c(1:length(cffidno)))

Yd.temp <- data.frame(cffidno = rep(unique(dat$cffidno),each=74), time = 1:74)
head(Yd.temp,n=12)

Y.epic <- merge(dat2,Yd.temp,by=c('cffidno','time'),all.y=TRUE)
nrow(dat)
nrow(Y.epic)
str(Y.epic)
head(Y.epic[,1:6],150 )

#################Readingin data for VisitAge matrix#############################

X <- matrix(Y.epic$VisitAge, M, 74, byrow=TRUE)
X[1,]##need to remove first row of column headers
dim(X)
X[1,]

#################Readingin data for Pa+ matrix#############################

Y <- matrix(Y.epic$cltpa, M, 74, byrow=TRUE)
Y[1:2,]##need to remove first row of column headers
dim(Y)

#################input variables for simulation#####################
#### checking for how many individuals we have NAs in the middle of followup
sum.na <- rep(NA,M)
k=rep(NA,M)

ids <- unique(Y.epic$cffidno) ## 103104 103125 103129 103145 103147
for (i in 1:M){
  na.indices <- which(Y.epic$VisitAge[Y.epic$cffidno==ids[i]] %in% NA)
  if (length(na.indices)==0){
    k[i] <- 74} else{
      k[i] <- min(na.indices)-1}
}

alpha = c(1,1)


#################################################################
####### model finding 1 fixed change point (wide format data), same c0 and u,  cp1 ~ dunif(min,max)
##---------------------------------------------------------------         
#################################################################
#### 

modelrancp <- "model { 
  for(i in 1:M){ 
        for(j in 1:k[i]){
  ### PA model
        Y[i,j] ~ dbin(p2[i,j],1)
        logit(p2[i,j]) <- logit(p[i,j,z[i]])
        logit(p[i,j,1]) <-  c0 + c[1] * (X[i,j]-cp1[i]) + c[2] * (X[i,j]-cp1[i]) * (2*step(X[i,j]-cp1[i])-1) + u[i]
        logit(p[i,j,2]) <-  c0 + B1*(X[i,j]-10) + u[i]
        }
        z[i]~dcat(pi[1:2])
        u[i] ~ dnorm(0,tau.u)
        cp1[i] ~ dnorm(cp1.mu,cp1.tau)
  }
  pi[1:2] ~ ddirch(alpha[])
  c0 ~ dnorm(0,0.0001)
	for (k in 1:2){
	   c[k] ~ dnorm(0,0.0001)	
	}
	B1<-c[1]-c[2]
  B2<-c[1]+c[2]
	tau.u ~ dgamma(0.001,0.001)
  cp1.mu ~ dnorm(0,0.001);
	cp1.tau ~ dgamma(0.001,0.001)
}"

####Observed DATA 
data <- dump.format(list(X=X, Y=Y, M=M, k=k, alpha=alpha)) 

###initial Values
inits1 <- dump.format(list(c0=-1, c=c(0.5,-0.5), pi=c(0.79,0.21), tau.u=1, cp1.mu=5, cp1.tau=1,
                           .RNG.name="base::Super-Duper", .RNG.seed=1)) 
inits2 <- dump.format(list(c0=-1.1, c=c(0.5,-0.5)-0.1, pi=c(0.8,0.2), tau.u=1, cp1.mu=5.1, cp1.tau=1.1,
                           .RNG.name="base::Super-Duper", .RNG.seed=2))

#### Run the model and produce plots
res1 <- run.jags(model=modelrancp, burnin=18000, sample=12000, 
                 monitor=c("c0", "c", "cp1", "B1", "B2", "u", "pi", "z", "tau.u"), 
                 data=data, n.chains=2, inits=c(inits1,inits2), thin=20, plots = FALSE)

sum1 <- summary(res1)
sum.df <- as.data.frame(sum1)
cp1 <- sum.df[4:326,4]
mean(cp1) # 2.911205
median(cp1) # 3.138959





