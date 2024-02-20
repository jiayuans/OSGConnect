#!/usr/bin/env Rscript

library(rjags)
library(coda)
library(R2jags)

# JAGS
set.seed(42) # Set a random seed for reproducibility of the simulation

samplesize <- 100 # Number of data points
b_x <- sort(rnorm(samplesize)) # x (explanatory variable)

int_true <- 30 # True intercept
slope_true <- 10 # True slope
mu <- int_true + slope_true * b_x # True means of normal distributions
sigma <- 5 # True standard deviation of normal distributions

b_y <- rnorm(samplesize, mean = mu, sd = sigma) # y (response variable)

df <- data.frame(b_x = b_x, b_y = b_y)
head(df)

jagsdata <- with(df, list(b_y = b_y, b_x = b_x, N = length(b_y)))

lm1_jags <- function(){
  # Likelihood:
  for (i in 1:N){
    b_y[i] ~ dnorm(mu[i], tau) # tau is precision (1 / variance)
    mu[i] <- alpha + beta * b_x[i]
  }
  # Priors:
  alpha ~ dnorm(0, 0.001) # intercept
  beta ~ dnorm(0, 0.001) # slope
  sigma ~ dunif(0, 10) # standard deviation
  tau <- 1 / (sigma * sigma) # sigma^2 doesn't work in JAGS
}

init_values <- function(){
  list(alpha = rnorm(1), beta = rnorm(1), sigma = runif(1))
}

params <- c("alpha", "beta", "sigma")

fit_lm1 <- jags(data = jagsdata, inits = init_values, parameters.to.save = params, model.file = lm1_jags,
                n.chains = 3, n.iter = 12000, n.burnin = 2000, n.thin = 10, DIC = F)

fit_lm1
