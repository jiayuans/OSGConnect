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
dat0