#'
#' @param x_t current state of stochastic oscillator
#' ######## currently unsure what $x_t$ and other params represent
#'

#setwd("C:/Users/Amy/src/Boettiger/regimeshifts")
#getwd()

may <- function(x_t,
                r = .5,
                K = 2, 
                Q = 5,
                H = .38,
                sigma = .04,
                a = 0.245) {

  # Deterministic mean
  d_mean <- x_t + x_t * r * (1 - x_t / K) - a * x_t ^ Q / (x_t ^ Q + H ^ Q)
  
  # Normal distribution
  #y_t1 <- ((sigma / (2 * pi)) ^ .5) * exp((-sigma(x - d_mean) ^ 2) / 2)
  ######## this doesn't quite work for saving a random variable distribution func...
  y_t1 <- function(z, mean, sd) {pnorm(z, mean, 1/(sd^2))}

  # the random variable x is truncated at 0
  ######## must find an interpretation of y_t1 such that some max() function can apply
  x_t1 <- max(y_t1, 0)
  
  # Return the newly-calculated x_{t+1}
  x_t1
}




## Now run tests with the function
######## will be moved to a different file later

num_days <- 1000

# Initialize variables
x <- numeric(num_days)
x[1] <- 1.2

# Simulate for num_days
for (day in 1:(num_days)) {
  x[day+1] <- may(x[day])
}


# For testing purposes
x[2] <- may(x[1])
