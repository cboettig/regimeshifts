# may
#
# stochastic oscillator which calculates mean growth based on a deterministic growth model
# implements mean (mu) in a random draw from norm distribution 
# returns a population density at next time step
##
may <- function(r = .5,         #growth rate
                x_t   ,         #x value at time t
                K = 2, 
                Q = 5, 
                H = .38, 
                sigma = .04,    #standard deviation
                a = 0.245)
{
  # Determinstic mean looks like standard R
  mu_t <- x_t + x_t * r * (1 - x_t / K)  - a * x_t ^ Q / (x_t ^ Q + H ^ Q)
  
  # stochastic implementation, now in base R
  #draws single value from distribution with mean mu and sd param sigma
  y_t1 <- rnorm(1, mu_t, sd = sigma) 
  
  #x value for n+1 time step
  x_t1 <- max(y_t1,0)
  
  #return x at new individual time step
  x_t1
}

