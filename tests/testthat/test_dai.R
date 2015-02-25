test_that("dai model runs", {
  max_days <- 30
  DF <- seq(0, 2000, length=max_days) # schedule for env degredation (increased dilution)
  y <- numeric(max_days)
   
  y[1] <- 1.76e5 # initial density
   
  for(day in 1:(max_days-1)){
   y[day+1] <- dai(y[day], DF = DF[day])
  }
  
  expect_is(y, "numeric")
  
})