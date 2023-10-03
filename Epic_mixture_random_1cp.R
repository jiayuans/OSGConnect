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
