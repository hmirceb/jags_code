# Appendix C for: Sampling design and estimates of observation error greatly reduce quasi-extinction probability in plant populations

#==================================================================================#
# JAGS stochastic exponential population model to calculate
# long-term population growth rates while considering the variance associated with
# observation error. This code also allows estimation of stochastic growth rates with data 
# sets that includ missing years of data. 
#=================================================================================#

#### START OF SCRIPT ####

model{
  # Estimates of population size and growth rates and comparison with observed population data:
  for (ii in 1:n_good_counts) {
    # Estimation of 
    log_N[good_data[ii]] <- sum(mn.loglams[(good_data[ii] - lags[good_data[ii]]):(good_data[ii]-1)])+
      log_M[good_data[ii] - lags[good_data[ii]]]+
      epsilon_t[good_data[ii] - lags[good_data[ii]]] - epsilon_t[good_data[ii]]
    
    log_M[good_data[ii]] ~ dnorm(log_N[good_data[ii]],
                                  (process_precision/lags[good_data[ii]])) # The process variance increases linearly with the size of gaps between data
  }
  
  # estimation of mean annual growth rates including plot effects, and estimation of annual observation errors
  for (ii in 1:n_years) {
    mn.loglams[ii] <- mean_log_lambda + 
      plot_randomeffect[plot_level[ii]]
    
    epsilon_t[ii] ~ dnorm(0, epsilon_oe_precision)
  }
  
  # generating plot random effects
  for(plot_iterator in 1:n_plots){
    plot_randomeffect[plot_iterator] ~ dnorm(0, plot_precision)
  }
  
  oe_gamma_rate <- 0.5*epsilon_oe_precision # conversion of OE precision to gamma parameter
  
  # use of duplicate count data to estimate OE precision
  for (ii in 1:n_dups) { 
    dupvars[ii] ~ dgamma(1/2, oe_gamma_rate)
  }  
  
  # Prior distributions for each parameter. Normal distributions are parameterized as
  # dnorm(mean, precision) and gamma distributions as dgamma(shape, rate).
  
  # Growth rate:
  mean_log_lambda ~ dnorm(0, 0.001) 
  
  # Observation error variance:
  epsilon_oe_precision ~ dgamma(0.001, 0.001) 
  
  # Process variance:
  process_precision ~ dgamma(0.001, 0.001)
  
  # Plot variance:
  # plot_precision ~ dgamma(0.001, 0.001)
  
  #  convert precisions to SD for easier interpretation.
  epsilon_oe_sd <- sqrt(1/epsilon_oe_precision)
  process_sd <- sqrt(1/process_precision)
  plot_sd <- sqrt(1/plot_precision)
}

#### END OF SCRIPT ####
