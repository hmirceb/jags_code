#==================================================================================#
# Bayesian meta-analysis model to calculate average growth rates between populations
#=================================================================================#

model{
  
  for(i in 1:N){
    P[i] <- 1/(SD[i]^2)
    Mean[i] ~ dnorm(regression_fitted[i], P[i])
    regression_fitted[i] <- intercept + mu_randomeffect[mu[i]] + taxon_randomeffect[taxon[i]]
  }
  
  # Nested random effect of taxon inside monitoring unit:
  # MU:
  for(i in 1:n_mu){mu_randomeffect[i] ~ dnorm(0, mu_precision)}
  # Taxon:
  for(i in 1:n_taxon){taxon_randomeffect[i] ~ dnorm(mu_randomeffect[mu[i]], taxon_precision)}
  
  
  # Prior distributions for each parameter. Normal distributions are parameterized as
  # dnorm(mean, precision) and gamma distributions as dgamma(shape, rate).
  intercept ~ dnorm(0, 1e-06)
  
  # Monitoring unit variance:
  mu_precision ~ dgamma(0.001, 0.001)
  mu_sd <- sqrt(1/mu_precision)
  
  # Taxon variance:
  taxon_precision ~ dgamma(0.001, 0.001)
  taxon_sd <- sqrt(1/taxon_precision)
  
}

# End of model