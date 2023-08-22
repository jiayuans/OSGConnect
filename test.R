#!/usr/bin/env Rscript

data(trees)
head(trees)
mod <- lm(Girth  ~ Height , data = trees)
summary(mod)
