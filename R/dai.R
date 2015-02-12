#' @details 
#' Dai et al model:
#' 
#' $n_{t+1} = n_t g(n_t + \epsilon n_t, \theta)$
#' $g(n_t) = n_{t+1} / n_t$
#' 
#' This model is based on two phases of daily
#' growth: a slow exponential growth phase at low cell densities, followed by a logistic growth
#' phase with a higher per capita growth rate at intermediate cell densities. This model has 5
#' parameters: T lag is the lag time before yeast cells start to grow after being transferred into new
#' media (the total time for daily growth is 23 hours). In the slow exponential phase, the population
#' grows with a constant per capita growth rate γ low . After the population reaches a threshold
#' density N c , the subsequent logistic growth is determined by γ high (γ high >γ low ) and the carrying
#' capacity K

g <- function(n_t, 
              t = 23, # Hours between serial dilutions
              epsilon = rnorm(1, 0, 0.15), 
              theta = c(gamma_high = 0.4, # hr^-1
                        gamma_low = 0.3, # hr^-1
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
  # DF <- theta[["DF"]]
  
  ## Dilute and start growing. Noise enters in the dilution process
  n_t <- n_t  * (1 + epsilon) / DF
    
  if(t < T_lag)
    n_t1 <- n_t
  if(n_t < 1e-20)
    n_t1 <- 0
  else {
    
    if(n_t < N_c){ # Needs to switch once n > N_c
      
      t_c <- log(N_c / n_t) / gamma_low + T_lag
      
      n_low <- n_t * exp((t_c - T_lag) * gamma_low) 
      
      ## NO! must integrate
      #n_high <- n_low * exp((t - t_c) * gamma_high * (1 - n_low / K ))
      
      # x_t = K exp(r t + K c_1) / (exp(r t + K c_1) - 1)
      # N_0 =  K exp( K c_1) / (exp( K c_1) - 1)
      # N_0 (exp( K c_1) - 1) =  K exp( K c_1) 
      # K exp( K c_1) - N_0 exp( K c_1) + N_0 = 0
      # exp( K c_1)( K  - N_0) + N_0 = 0
      # exp(K c_1) = N_0 / (K - N_0) 
      # = 1/(K/N_0 - 1) = B
      # 
      # x_t = K B exp(r t )  / (B exp(r t) - 1)
      # x_t = K exp(r t )  / (exp(r t) - 1/B)
      # x_t = K exp(r t )  / (exp(r t) - 1 + K/N_0)      
      # x_t = K / (1 +   (K/N_0 - 1) exp( - r t) )

      # Why are these then not the same??
      
      # n_t1 <- K /(1 + (K/n_low - 1) * exp(- gamma_high * tau))
      # n_t1 <- n_t * K * exp(gamma_high * tau )  / (n_t * exp(gamma_high * tau) - n_t + K)
      
      
      
      tau <- t - t_c        
      n_t <- n_low
      n_t1 <- K /(1 + (K/n_t - 1) * exp(- gamma_high * tau))
      #n_t1 <- n_t * K * exp(gamma_high * tau )  / (n_t * exp(gamma_high * tau) - n_t + K)
      
      
      
    } else if(n_t >= N_c){
      #n_t1 <- n_t * exp((t - T_lag) * gamma_high * (1 - n_t / K ))
      ## NO! must integrate this instead: 
      
      tau <- t - T_lag
      n_t1 <- K /(1 + (K/n_t - 1) * exp(- gamma_high * tau))
      #n_t1 <- n_t * K * exp(gamma_high * tau )  / (n_t * exp(gamma_high * tau) - n_t + K)
      
      
    }
  }
  # return n_{t+1}
  n_t1
}





max_days <- 30
n <- numeric(max_days)
x <- numeric(max_days)
n[1] <- 2e3
x[1] <- 1e3
for(day in 1:(max_days-1)){
  n[day+1] <- g(n[day])
  x[day+1] <- g(x[day])
}


df <- data.frame(t = 1:max_days, n = n, x = x)
library(ggplot2)
ggplot(df) + 
  geom_line(aes(t, n), col = 1) + 
  geom_line(aes(t, x), col = 2) +
  scale_y_log10()


# Stepwise changes
DF <- as.numeric(sapply(seq(0, 2000, length=9), rep, 40))

# continuous linear increase
DF <- seq(0, 2000, length=1e3)

max_days <- length(DF)
y <- numeric(max_days)

y[1] <- 1.76e5
for(day in 1:(max_days-1)){
  y[day+1] <- g(y[day], DF = DF[day])
}

qplot(seq_along(y), y)

