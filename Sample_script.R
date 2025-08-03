# Appendix D for: Sampling design and estimates of observation error greatly reduce quasi-extinction probability in plant populations

#### START OF SCRIPT ####

#====================================================================================================#
# Script for fitting a Bayesian stochastic exponential population model to calculate
# long-term population growth rates while considering the variance associated with observation error
# and also allowing for missing years of data
#====================================================================================================#

#-----------------#
#### Libraries ####
#----------------#

library(tidyverse) # Data manipulation
library(rjags) # JAGS interface for R
library(runjags) # For easier interaction with JAGS

#------------#
#### Data ####
#------------#
# Example data from a real population of Acatea spicata under monitoring in the Adopt a Plant citizen science program:
count_data <- data.frame(taxon = 'Acatea_spicata',
                         N = c(62, 51, 52, 41, 59, 59, 98, 21, 46, 25, 17, 19, 20, 18, 37, 20, 29, 22, 26, 14, 35, 14, 9, 13, 
                               11, 15, 14, 12, 25, 29, 27, 71, 27, 30, 49, 25, 13, 17, 15, 13, 13, 14),
                         N2 = c(NA, NA, NA, 42, 57, NA, 102, 24, 49, NA, NA, NA, 19, 18, NA, NA, NA, 22, 28, NA,
                                42, 12, 13, NA, NA, NA, 13, 12, NA, NA, NA, 75, 26, NA, 60, 28, 15, NA, NA, NA, 13, 14),
                         plot = rep(c('plot_1', 'plot_2', 'plot_3'), each = 14),
                         lags = c(-1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, # Tags for starting year (-1), years with valid counts (>=1) and years with missing data (0) 
                                  -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, # Positive numbers indicate the number of transitions between years with valid data,
                                  -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) # 1 corresponds to consecutive valid years, 2 to two transition, etc.
)

# Prepare data for model:
log_M <- log(count_data$N+1) # Add 1 to count data to avoid possible problems with zeros and log transform
n_years <- length(log_M) # Number of monitoring years across all plots
lags <- count_data$lags # Get tags for starting year, years with valid counts and years with missing data 
good_data <- which(lags > 0) # Get years with valid count data
n_good_counts <- length(good_data) # Number of years with valid count data
plot_level <- as.factor(count_data$plot) # Unique plots as factor
n_plots <- length(unique(count_data$plot)) # Number of unique plots
      
#----------------------------------#
#### Observation error variance ####
# ---------------------------------#
# DF with years with a double census:
doubles <- count_data %>% 
  filter(!is.na(N2)) %>%
  mutate(across(N:N2, \(x) log(x+1))) # Log transform doubles

# DF with double censuses for that specific population:
doubles_pop <- doubles %>% 
  filter(!is.na(N2))

# DF with double censuses for that species overall:
doubles_sp <- doubles %>%
  filter(taxon == unique(count_data$taxon))

# If that particular population has at least a double census, use it:
if(dim(doubles_pop)[1] != 0){ # Check that the DF with double censuses has rows
  dups <- doubles_pop %>%
    dplyr::select(N, N2) # Get only the count data
  
  dupvars <- apply(dups, 1, var) # Variance between each census
  n_dups <- length(dupvars) # Number of double censuses
}

# If that particular population does not have a double census, but the species
# has at least one elsewhere, use them:
if(dim(doubles_pop)[1] == 0 & # Check that DF of the population is empty
   dim(doubles_sp)[1] != 0){ # and DF of species has data
  dups <- doubles_sp %>%
    dplyr::select(N, N2)
  
  dupvars <- apply(dups, 1, var)
  n_dups <- length(dupvars)
}

#---------------------#
#### Model fitting ####
#---------------------#
      
# JAGS MCMC settings:
nc <- 4 # Number of MCMC chains
nb <- 100000 # Number burn-in samples
nt <- 10 # Thinning interval
ns <- 1000000 # Number samples by chain
nad <- 10000 # Number of adaptive samples

# Variables to monitor in model:
monitor <- c('mean_log_lambda',
             'process_sd', 
             'epsilon_oe_sd', 
             'plot_sd')

# Bundle data for JAGS:
data = list(log_M = log_M,
            n_years = n_years,
            lags = lags, 
            good_data = good_data,
            n_good_counts = n_good_counts,
            plot_level = plot_level,
            n_plots = n_plots,
            n_dups = n_dups,
            dupvars = dupvars)

# Function to set up random initial values (inits) for each chain
init.vals <- function(){
  mean_log_lambda <- rnorm(1, 0, 1) # Init for mean growth rate
  process_precision <- runif(1, 0, 10) # Init for process precision (variance)
  plot_precision <- runif(1, 0, 10) # Init for plot precision (variance)
  epsilon_oe_precision <- runif(1, 0, 10) # Init for OE precision (variance)
  seed <- sample.int(1000000, 1) # Random seed to initialize each pseudo random number generator (PRNG)
  name <- sample(c('base::Wichmann-Hill', # Randomly choose one of these PRNG
                    'base::Marsaglia-Multicarry',
                    'base::Super-Duper',
                    'base::Mersenne-Twister'),
                 1)
  out <- list(.RNG.name = name,
              .RNG.seed = seed, 
              mean_log_lambda = mean_log_lambda, 
              process_precision = process_precision,
              plot_precision = plot_precision,
              epsilon_oe_precision = epsilon_oe_precision)
}

# Fit model:
population_model <- run.jags('SI_Model.r', 
                            n.chains = nc,
                            burnin = nb,
                            thin = nt,
                            sample = ns,
                            adapt = nad,
                            method = 'parallel',
                            monitor = monitor,
                            data = data,
                            inits = init.vals)
population_model # Model summary
plot(population_model) # Convergence plots

#### END OF SCRIPT ####
