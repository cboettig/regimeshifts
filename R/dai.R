#' dai
#' 
#' Growth of yeast population under a serial dilution regime, based on Dai et al 2012 Science paper
#' @param n_t The current population size
#' @param t Number of hours between serial dilutions
#' @param epsilon random shock for dilution. Note: only applied
#'  once per dilution (e.g once per function call)
#' @param theta parameters for the model; defaults to estimates from Dai et al; see details
#' @param DF Dilution factor, the environmental variable manipulated to cause the bifurcation
#' @return the population size the next day, after serial dilution and growth
#' @details 
#' Dai et al model:
#' 
#' $n_{t+1} = n_t g(n_t + \\epsilon n_t, \\theta)$
#' $g(n_t) = n_{t+1} / n_t$
#' 
#' From the supplement of Dai et al (2012), Science:
#' "This model is based on two phases of daily
#' growth: a slow exponential growth phase at low cell densities, followed by a logistic growth
#' phase with a higher per capita growth rate at intermediate cell densities. This model has 5
#' parameters: T lag is the lag time before yeast cells start to grow after being transferred into new
#' media (the total time for daily growth is 23 hours). In the slow exponential phase, the population
#' grows with a constant per capita growth rate gamma_low . After the population reaches a threshold
#' density N_c , the subsequent logistic growth is determined by gamma_high (gamma_high > gamma_low) and the carrying
#' capacity K"
#' 
#' @examples
#' 
#' max_days <- 30
#' DF <- seq(0, 2000, length=max_days) # schedule for env degredation (increased dilution)
#' y <- numeric(max_days)
#' 
#' y[1] <- 1.76e5 # initial density
#' 
#' for(day in 1:(max_days-1)){
#' y[day+1] <- dai(y[day], DF = DF[day])
#' }
#' 
#' plot(seq_along(y), y)
#' 
#' @export
#' 
dai <- function(n_t, 
              t = 23, # Hours between serial dilutions
              epsilon = rnorm(1, 0, 0.15), 
              theta = c(gamma_high = 0.439, # hr^-1
                        gamma_low = 0.309, # hr^-1
                        T_lag = 2.97, # hr
                        N_c = 2.76e2, # cells/μl
                        K = 1.76e5),  # cells/μl
              DF = 600)  # dilution factor
{
  ## avoid referencing these repeatedly for readability and speed
  gamma_high <- theta[["gamma_high"]]
  gamma_low <- theta[["gamma_low"]]
  T_lag <- theta[["T_lag"]]
  N_c <- theta[["N_c"]]
  K <- theta[["K"]]
  
  
  ## Dilute and start growing. (Stochasticity enters only via the dilution process)
  n_t <- n_t  * (1 + epsilon) / DF
  
  # Lag phase, could have been scaled out of the model
  if(t < T_lag)
    n_t1 <- n_t
  
  ## Numerical happiness
  if(n_t < 1e-20)
    n_t1 <- 0
  
  ## Actual model
  else {
    
    if(n_t < N_c){ # Needs to switch once n > N_c
      
      ## Analytically find out how long before we leave the low-growth regime
      t_c <- log(N_c / n_t) / gamma_low + T_lag
      
      ## Um, now this should just be equal to N_c
      n_low <- n_t * exp((t_c - T_lag) * gamma_low) 
      
      ## Spend remaining time in > N_c growth regime:      
      tau <- t - t_c        
      n_t <- n_low
      n_t1 <- K  / (1 + (K / n_t - 1) * exp(- gamma_high * tau))
      
      
    } else if(n_t >= N_c){
      ## simpler if we're always in the high-growth regime:
      tau <- t - T_lag
      n_t1 <- K / (1 + (K / n_t - 1) * exp(- gamma_high * tau))
      
    }
  }
  # And now we can return n_{t+1}
  n_t1
}
