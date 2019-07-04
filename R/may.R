# may
#
# stochastic oscillator which calculates mean growth based on a deterministic growth model
# implements mean (mu) in a random draw from norm distribution 
# returns a population density at next time step
##
may <- function(r = .5, #growth rate
                x0 = 1, #initial x value
                K = 2, 
                Q = 5, 
                H = .38, 
                sigma = .04, #standard deviation
                a = 0.245, 
                N = 1e4)
{
  #sets x, y, and mu as same length N
  x <- numeric(N)  
  y <- numeric(N)
  mu <- numeric(N)
  x[1] <- x0
  for(t in 1:(N-1)){
    
    # Determinstic mean looks like standard R
    mu[t] <- x[t] + x[t] * r * (1 - x[t] / K)  - a * x[t] ^ Q / (x[t] ^ Q + H ^ Q)
    
    # stochastic implementation, now in base R
    #draws single value from distribution with mean mu and sd param sigma
    y[t+1] <- rnorm(1, mu[t], sd = sigma) 
    
    #x value for n+1 time step
    x[t+1] <- max(y[t+1],0)
  }
  #return x
  x 
}