#!/usr/bin/env Rscript

library(ggplot2)
library(lme4)

set.seed(222)

N <- 150
nobs <- 6
sigma <- 20
sig00 <- 10
b0 <- 80
b1 <- -6

ID = rep(1:N, each = nobs)
Time = rep(0:5, 150)

dat <- data.frame(ID, Time)

u0 <- rnorm(N, 0, sig00)
dat$u0 <- rep(u0, each=nobs)

dat$err <- rnorm(900, 0, sigma)

dat$Stress <- (b0 + u0)  + b1*dat$Time + dat$err

dat <- dat[,-c(3:4)]
summary(dat) 
ggplot(dat, aes(x = Time, y = Stress, group=ID, color = ID)) +
  geom_line(alpha=.3) + labs(title="Stress Levels of Students Across Pet Therapy Intervention Program") +
  theme(legend.position = "none")

## Model 0: Linear regression (fixed-effects only)
model0 <- lm(Stress ~ Time, data=dat)
summary(model0)

## Model 1: Random intercept model
model1 <- lmer(Stress ~ Time + (1|ID), data=dat, REML=FALSE)
summary(model1)
