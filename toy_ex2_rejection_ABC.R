# Load required libraries
library(EasyABC)
library(ggplot2)

set.seed(0)

##### Simulate Observed Data #####
n <- 250
true_mu <- 0
true_sigma <- 5
observed_data <- rnorm(n, mean = true_mu, sd = true_sigma)

##### Metropolis-Hastings Sampling #####
likelihood <- function(mu1, sigma1, mu2, sigma2) {
  prod(dnorm(observed_data, mean = mu1, sd = sigma1) /
         dnorm(observed_data, mean = mu2, sd = sigma2))
}

mh_sampler <- function(iterations, init) {
  samples <- matrix(NA, nrow = iterations, ncol = 2)
  samples[1, ] <- init
  accept_count <- 0
  
  for (i in 2:iterations) {
    if (i %% (iterations / 100) == 0) {
       cat(sprintf("Progress: %d%%\n", i / (iterations / 100)))
    }
    proposal <- c(rnorm(1, mean = samples[i-1, 1], sd = 0.5),
                  runif(1, max(0.001, samples[i-1, 2] - 0.5), min(10, samples[i-1, 2] + 0.5)))
    
    if (proposal[2] <= 0) next  # Ensure positive sigma
    
    likelihood_ratio <- likelihood(proposal[1], proposal[2], samples[i - 1, 1], samples[i - 1, 2])
    acceptance_prob <- min(1, likelihood_ratio)
    
    if (runif(1) < acceptance_prob) {
      samples[i, ] <- proposal
      accept_count <- accept_count + 1
    } else {
      samples[i, ] <- samples[i - 1, ]
    }
  }
  
  cat("MH Acceptance Rate:", accept_count / iterations, "\n")
  return(samples)
}

iterations <- 5000000
init <- c(mean(observed_data), sd(observed_data))
mh_samples <- mh_sampler(iterations, init)
mh_samples <- mh_samples[5001:nrow(mh_samples), ]  # Burn-in

mu_samples_mh <- mh_samples[, 1]
sigma_samples_mh <- mh_samples[, 2]

##### Plot MH Posterior #####
par(mfrow = c(1, 2))
hist(mu_samples_mh, breaks = 50, col = "skyblue", main = "Posterior μ (MH)", xlab = "μ", probability = TRUE)
hist(sigma_samples_mh, breaks = 50, col = "salmon", main = "Posterior σ (MH)", xlab = "σ", probability = TRUE)




####ABC posterior

# Define prior distributions
prior <- list(c("normal", 0, 5), c("unif", 0, 10))

# Define the simulation model
model <- function(params) {
  mu <- params[1]
  sigma <- params[2]
  y <- rnorm(n, mean = mu, sd = sigma)
  return(c(mean(y), sqrt(var(y))))  # Summary statistics: mean and std deviation
}

# Perform ABC rejection
abc_results <- ABC_rejection(
  model = model,
  prior = prior,
  summary_stat_target = c(mean(observed_data), sqrt(var(observed_data))),
  nb_simul = 100000,
  tol = 1,
  use_seed = TRUE,
  progress_bar = TRUE
)

# Extract posterior samples and remove burn-in
posterior_samples <- abc_results$param
mu_samples_rejection <- posterior_samples[,1]
sigma_samples_rejection <- posterior_samples[,2]


hist(mu_samples_rejection, breaks = 50, col = "skyblue", main = "Estimated Posterior Distribution of μ", xlab = "μ", probability = TRUE)
hist(sigma_samples_rejection, breaks = 50, col = "lightcoral", main = "Estimated Posterior Distribution of σ", xlab = "σ", probability = TRUE)





grid_size <- 40000
grid_mu <- seq(-10, 10, length.out = grid_size)
grid_sigma <- seq(0, 10, length.out = grid_size)  # sigma must be non-negative

# Initialize probability distributions
mh_mu_prob <- numeric(grid_size)
rejection_mu_prob <- numeric(grid_size)
mh_sigma_prob <- numeric(grid_size)
rejection_sigma_prob <- numeric(grid_size)

# Compute histogram probabilities for Metropolis-Hastings (mu)
mu_indices_mh <- round((mu_samples_mh + 10) / 20* (grid_size - 1)) + 1
mu_indices_mh <- pmax(1, pmin(grid_size, mu_indices_mh))  # Clip to valid range
mh_mu_prob[mu_indices_mh] <- mh_mu_prob[mu_indices_mh] + 1 / length(mu_samples_mh)

# Compute histogram probabilities for Estimation (mu)
mu_indices_rejection <- round((mu_samples_rejection + 10) / 20 * (grid_size - 1)) + 1
mu_indices_rejection <- pmax(1, pmin(grid_size, mu_indices_rejection))
rejection_mu_prob[mu_indices_rejection] <- rejection_mu_prob[mu_indices_rejection] + 1 / length(mu_samples_rejection)

# Compute histogram probabilities for Metropolis-Hastings (sigma)
sigma_indices_mh <- round((sigma_samples_mh+0) / 10* (grid_size - 1)) + 1
sigma_indices_mh <- pmax(1, pmin(grid_size, sigma_indices_mh))
mh_sigma_prob[sigma_indices_mh] <- mh_sigma_prob[sigma_indices_mh] + 1 / length(sigma_samples_mh)

# Compute histogram probabilities for Estimation (sigma)
sigma_indices_rejection <- round((sigma_samples_rejection+0) / 10 * (grid_size - 1)) + 1
sigma_indices_rejection <- pmax(1, pmin(grid_size, sigma_indices_rejection))
rejection_sigma_prob[sigma_indices_rejection] <- rejection_sigma_prob[sigma_indices_rejection] + 1 / length(sigma_samples_rejection)

# Compute L1 errors (no need to multiply by bin width)
error_mu <- sum(abs(mh_mu_prob - rejection_mu_prob))
error_sigma <- sum(abs(mh_sigma_prob - rejection_sigma_prob))

# Final combined error
total_error <- sqrt(error_mu^2 + error_sigma^2)

# Print errors
print(paste("Error in mu:", error_mu))
print(paste("Error in sigma:", error_sigma))
print(paste("Total error:", total_error))

