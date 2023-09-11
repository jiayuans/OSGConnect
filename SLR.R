#!/usr/bin/env Rscript
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

#      mu.vect sd.vect   2.5%    25%    50%    75%  97.5%  Rhat n.eff
# alpha  29.543   0.462 28.634 29.227 29.544 29.869 30.447 1.001  2700
# beta   10.101   0.447  9.215  9.804 10.110 10.398 10.988 1.001  2900
# sigma   4.607   0.332  4.024  4.378  4.583  4.828  5.303 1.002  1600

# 2.5% and 97.5% quantiles: 95% credible interval
# n.eff: effective sample size
# Rhat: measure of how well the Markov chains have mixed and should ideally have a value very close to 1
