#==================================================================================#
# Bayesian meta-analysis model to calculate average quasi-extinction probabilities
# between populations
#=================================================================================#

model{
  
  for(i in 1:N){
    # Precision of quasi-extinction probabilities
    P[i] <- SD[i]^2
    
    phi[i] <- log(1/P[i]) 
    Mean[i] ~ dbeta(alpha[i], beta[i])
    
    alpha[i] <- mu[i] * phi[i]
    beta[i] <- (1-mu[i]) * phi[i]
    
    # Logit link function
    logit(mu[i]) <- intercept + MU_DYN_randomeffect[MU_DYN[i]] + TAXON_randomeffect[TAXON[i]]
    
  }
  
  # Nested random effect of taxon inside monitoring unit:
  # MU:
  for(i in 1:nMU){MU_DYN_randomeffect[i] ~ dnorm(0, MU_DYN_precision)}
  # Taxon:
  for(i in 1:nTAXON){TAXON_randomeffect[i] ~ dnorm(MU_DYN_randomeffect[MU_DYN[i]], TAXON_precision)}
  
  # Prior distributions for each parameter. Normal distributions are parameterized as
  # dnorm(mean, precision) and gamma distributions as dgamma(shape, rate).
  intercept ~ dnorm(0, 1e-06)
  
  # Monitoring unit variance:
  MU_DYN_precision ~ dgamma(0.001, 0.001)
  MU_DYN_sd <- sqrt(1/MU_DYN_precision)
  
  # Taxon variance:
  TAXON_precision ~ dgamma(0.001, 0.001)
  TAXON_sd <- sqrt(1/TAXON_precision)
}

# End of model

