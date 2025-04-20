source("H.r")
source("classify.r")
source("rhathat.r")

mu = 2.0
v = 4.0
n = 250

xObs = rnorm(n,mean=mu,sd=sqrt(v))

nIt = 100000
sigmaProposal = c(0.1,0.1) # proposal standard deviations
theta = c(3,3.9)  # start value
muCorrect = c()
vCorrect = c()

for (i in c(1:nIt)) 
{
  thetaProposal = theta
  thetaProposal[1] = theta[1] + sigmaProposal[1] * rnorm(1)
  thetaProposal[2] = theta[2] + sigmaProposal[2] * rnorm(1)
  
  mean = theta[1]
  var = theta[2]
  meanProposal = thetaProposal[1]
  varProposal = thetaProposal[2]
  
  
  if (varProposal > 0.0) # otherwise reject
  { 
    logr = log(var) - log(varProposal) - (n/2)*log(varProposal) + (n/2)*log(var)
    logr = logr - (1/(2*varProposal))*sum((xObs-meanProposal)^2) + (1/(2*var))*sum((xObs-mean)^2)
    
    if (logr > 0 || runif(1) <= exp(logr)) # accept proposal
    {  
      theta = thetaProposal
    }
  }
  
  muCorrect = c(muCorrect,theta[1])
  vCorrect = c(vCorrect,theta[2])
}

# simulate from the approximate posterior 
library(class)  

theta = c(3.0,3.9)  # start value
classifyParam = NULL # Control parameters for logistic regression 

q = 0.5
sigma0 = 0.0065
sigma1 = 0.0065

y0Obs = 0.5
y1Obs = 0.5


muApprox = c()
vApprox = c()
signApprox = c()


xGenerated = rnorm(n,mean = theta[1],sd = sqrt(theta[2]))
r0 = rhathat(0,y0Obs,classifyParam,q,sigma0,xObs,theta[1],sqrt(theta[2]))
r1 = rhathat(1,y1Obs,classifyParam,q,sigma1,xObs,theta[1],sqrt(theta[2]))
fhat = abs(r0 * r1)




for (i in c(1:nIt)) 
{
  thetaProposal = theta
  thetaProposal[1] = theta[1] + sigmaProposal[1] * rnorm(1)
  thetaProposal[2] = theta[2] + sigmaProposal[2] * rnorm(1)
  
  mean = theta[1]
  sd = sqrt(theta[2])
  
  if (thetaProposal[2] > 0.0) # otherwise reject
  { 
    meanProposal = thetaProposal[1]
    sdProposal = sqrt(thetaProposal[2])
    
    r0Proposal = rhathat(0,y0Obs,classifyParam,q,sigma0,xObs,meanProposal,sdProposal)
    r1Proposal = rhathat(1,y1Obs,classifyParam,q,sigma1,xObs,meanProposal,sdProposal)
    fhatProposal = abs(r0Proposal * r1Proposal)
    
    rr = (sd * sd) / (sdProposal * sdProposal)
    rr = rr * fhatProposal / fhat
    
    print(c(i,rr,mean,sd))
    
    if (runif(1) <= rr) 
    { 
      theta = thetaProposal
      r0 = r0Proposal
      r1 = r1Proposal
      fhat = fhatProposal
    }
  }
  
  muApprox = c(muApprox,theta[1])
  vApprox = c(vApprox,theta[2])
  signApprox = c(signApprox,sign(r0 * r1))
}


# --- Posterior summaries ---
burnin <- 500
idx <- (burnin + 1):nIt

muExact_mean <- mean(muCorrect[idx])
muApprox_mean <- sum(muApprox[idx] * signApprox[idx]) / sum(signApprox[idx])

varExact_mean <- mean(vCorrect[idx])
varApprox_mean <- sum(vApprox[idx] * signApprox[idx]) / sum(signApprox[idx])

sdExact_mean <- mean(sqrt(vCorrect[idx]))
sdApprox_mean <- sum(sqrt(vApprox[idx]) * signApprox[idx]) / sum(signApprox[idx])

muExact_sd <- sd(muCorrect[idx])
muApprox_sd <- sqrt(sum((muApprox[idx] - muApprox_mean)^2 * signApprox[idx]) / sum(signApprox[idx]))

varExact_sd <- sd(vCorrect[idx])
varApprox_sd <- sqrt(sum((vApprox[idx] - varApprox_mean)^2 * signApprox[idx]) / sum(signApprox[idx]))

sdExact_sd <- sd(sqrt(vCorrect[idx]))
sdApprox_sd <- sqrt(sum((sqrt(vApprox[idx]) - sdApprox_mean)^2 * signApprox[idx]) / sum(signApprox[idx]))

# --- Print results ---
cat("Posterior Mean of mu:\nExact  :", muExact_mean, "\nApprox:", muApprox_mean, "\n\n")
cat("Posterior Mean of variance:\nExact  :", varExact_mean, "\nApprox:", varApprox_mean, "\n\n")
cat("Posterior Mean of std dev:\nExact  :", sdExact_mean, "\nApprox:", sdApprox_mean, "\n\n")

cat("Posterior SD of mu:\nExact  :", muExact_sd, "\nApprox:", muApprox_sd, "\n\n")
cat("Posterior SD of variance:\nExact  :", varExact_sd, "\nApprox:", varApprox_sd, "\n\n")
cat("Posterior SD of std dev:\nExact  :", sdExact_sd, "\nApprox:", sdApprox_sd, "\n\n")


# --- Grid-based posterior estimation using rhathat ---
grid_vals = seq(1.0, 3.0, by = 0.2)
grid_points = expand.grid(mu = grid_vals, sigma = grid_vals)
posterior_estimate = numeric(nrow(grid_points))
nrow(grid_points)
for (i in 1:nrow(grid_points)) 
{
  print(i)
  mu_i = grid_points$mu[i]
  sigma_i = grid_points$sigma[i]
  
  r0 = rhathat(0, y0Obs, classifyParam, q, sigma0, xObs, mu_i, sigma_i)
  r1 = rhathat(1, y1Obs, classifyParam, q, sigma1, xObs, mu_i, sigma_i)
  
  posterior_estimate[i] = abs(r0 * r1)
}

# Normalize estimated posterior
posterior_estimate = posterior_estimate / sum(posterior_estimate)

# --- True posterior using exact likelihood + prior ---
true_posterior = numeric(nrow(grid_points))
for (i in 1:nrow(grid_points)) 
{
  print(i)
  mu_i = grid_points$mu[i]
  sigma_i = grid_points$sigma[i]
  
  if (sigma_i > 0 && sigma_i < 10)
  {
    loglik = sum(dnorm(xObs, mean = mu_i, sd = sigma_i, log = TRUE))
    logprior_mu = dnorm(mu_i, mean = 0, sd = 5, log = TRUE)           # Prior: N(0,25)
    logprior_sigma = log(1 / 10)                                       # Prior: Uniform(0,10)
    
    true_posterior[i] = exp(loglik + logprior_mu + logprior_sigma)
  }
}
true_posterior = true_posterior / sum(true_posterior)

# --- L1 Error (Total Variation Distance) ---
l1_error = sum(abs(posterior_estimate - true_posterior))
cat("Grid-based Posterior L1 Error (Total Variation Distance):", l1_error, "\n")

# PLot the densities
library(ggplot2)

grid_points$Posterior_Estimate <- posterior_estimate
grid_points$True_Posterior <- true_posterior

ggplot(grid_points, aes(mu, sigma)) +
  geom_tile(aes(fill = Posterior_Estimate)) +
  labs(title = "Estimated Posterior (Classifier-based)", fill = "Density") +
  theme_minimal()

ggplot(grid_points, aes(mu, sigma)) +
  geom_tile(aes(fill = True_Posterior)) +
  labs(title = "True Posterior (Likelihood + Prior)", fill = "Density") +
  theme_minimal()

