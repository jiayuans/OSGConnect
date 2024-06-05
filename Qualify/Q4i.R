#!/usr/bin/env Rscript

library(tidyverse)
library(dplyr)
library(Matrix)
library(broom)
library(MASS)


fmri_phys_dat <- read.csv("C:/UCHealth/Qualify/takehome_QE_2024/fmri_phys_func.csv")


# Function of the proposed method from (g)
simulation <- function(outcome.df) {
  # Predictors (measures of functional connectivity)
  predictor.df <- fmri_phys_dat[,17:46]
  
  # Number of outcomes (M) and predictors (P)
  M <- ncol(outcome.df)
  P <- ncol(predictor.df) 
  N <- nrow(fmri_phys_dat)
  
  # Fit regression models
  model_results_list <- list()
  cov_matrices <-  list()
  for (m in 1:M) {
    Y <- outcome.df[, m]
    model <- lm(Y ~ ., data = predictor.df)
    model_results <- tidy(model, conf.int = TRUE) %>%
      mutate(Y = paste0("Y.", m))
    
    model_results_list[[m]] <- model_results
    cov_matrices[[m]] <- vcov(model)
  }
  
  model_results_df <- bind_rows(model_results_list)
  
  comb_cov_matrices <- bdiag(cov_matrices)
  comb_cov_matrices <- as.matrix(comb_cov_matrices)
  
  # Extract beta hat and se hat
  beta_hat <- model_results_df %>% dplyr::select(Y, term, estimate)
  beta_hat_matrix <- beta_hat %>% mutate(term = factor(term, levels = unique(term))) %>%
    spread(key = term, value = estimate) %>% dplyr::select(-Y)
  beta_hat1 <- as.numeric(unlist(beta_hat[,3]))
  
  se_hat <- model_results_df %>% dplyr::select(Y, term, std.error)
  se_hat_matrix <- se_hat %>% mutate(term = factor(term, levels = unique(term))) %>%
    spread(key = term, value = std.error) %>% dplyr::select(-Y)
  
  # Ensure beta_hat_matrix rows match cov_matrices length
  beta_hat_list <- split(beta_hat_matrix, 1:nrow(beta_hat_matrix))
  
  # Number of simulations
  B <- 100000
  set.seed(123)
  # Function to simulate betas
  simulate_max_deviation <- function(beta_hat, cov_matrix, B) {
    simulated_betas <- mvrnorm(B, mu = beta_hat, Sigma = cov_matrix)
    return(simulated_betas)
  }
  
  # Simulate max deviations for each model
  max_devs_list <- simulate_max_deviation(beta_hat1, comb_cov_matrices, B)
  
  # Calculate nu
  nu <- apply(max_devs_list,1, function(x)
    max(abs(x-beta_hat1)/sqrt(diag(comb_cov_matrices))))
  
  # Quantile calculation
  alpha_grid <- seq(1/B, 1, length.out = B)
  quantiles <- as.numeric(quantile(nu, probs = 1 - 0.05))
  
  beta_hat_matrix <- as.matrix(beta_hat_matrix)
  se_hat_matrix <- as.matrix(se_hat_matrix)
  
  # Adjusted p-value calculation
  adj_p <- matrix(NA, nrow = M, ncol = P+1)
  for (i in 1:M) {
    for (k in 1:(P+1)) {
      alpha_i <- which((beta_hat_matrix[i, k] - quantiles * se_hat_matrix[i,k] < 0) &
                         (beta_hat_matrix[i, k] + quantiles * se_hat_matrix[i,k] > 0))
      if (length(alpha_i) > 0) {
        adj_p[i, k] <- 1 - min(alpha_grid[alpha_i])
      } else {
        adj_p[i, k] <- 1 / B
      }
    }
  }
  adj_p_num <- as.numeric(adj_p)
  return(adj_p_num)
}



# Function to set up the simulation
simulate_and_fit <- function(num_simulations) {
  # Outcomes (measures of physical function)
  outcome.df <- fmri_phys_dat[,2:16]
  
  # Predictors (measures of functional connectivity)
  predictor.df <- fmri_phys_dat[,17:46]
  
  FP <- 0
  M <- ncol(outcome.df)
  P <- ncol(predictor.df) 
  N <- nrow(fmri_phys_dat)
  
  # Fit regression models for the residual standard error
  res.se <- rep(0,M )
  for (m in 1:M) {
    Y <- outcome.df[, m]
    model <- lm(Y ~ ., data = predictor.df)
    res.se[m] <- summary(model)$sigma
  }
  
  for (i in 1:num_simulations) {
    outcome.df.simulated <- matrix(NA, N, M)
    # Simulate outcomes using the normal distribution of residuals 
    # with sd estimated from residual standard error from the regressions
    for (j in 1:M) {
      outcome.df.simulated[, j] <- rnorm(N,0,res.se[j])
    }
    adj_p_num <- simulation(outcome.df.simulated) 
    
    # Check for false positives
    if (any(adj_p_num < alpha)) {
      FP <- FP + 1
    }
  }
  # Family-wise error rate
  return(FP / num_sim)
  
}

num_sim <- 5000
alpha <- 0.05
FWER <- simulate_and_fit(num_sim)
FWER

cat("Family-wise Type-I Error Rate:", FWER, "\n")
