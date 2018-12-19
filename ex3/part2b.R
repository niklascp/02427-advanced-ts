library(ctsmr)
library(ggplot2)
load("Exercise3.RData")

sdeImproved <- function(data){
  # Generate a new object of class ctsm
  model = ctsm()
  # Add a system equation and thereby also a state
  # Gv in Ti: Aw/Ci*Gv or Tm: Aw/Cm*Gv
  model$addSystem(dTi1 ~  1/Ci*(1/Rin*(yTi2-Ti1) + 1/Ria*(Ta-Ti1) + 1/Rim*(Tm-Ti1)  + Ph1 + Aw*Gv)*dt + exp(sigma1)*dw1)
  model$addSystem(dTm ~  1/Cm*(1/Rim*(Ti1-Tm))*dt + exp(sigma2)*dw2)
  # Set the names of the inputs
  model$addInput(Ta,yTi2,bs1,bs5,Gv,Ph1)
  # Set the observation equation: Ti is the state, yTi is the measured output
  model$addObs(yTi1 ~ Ti1)
  # Set the variance of the measurement error
  model$setVariance(yTi1 ~ exp(e1))
  ##----------------------------------------------------------------
  # Set the initial value (for the optimization) of the value of the state at the starting time point
  model$setParameter(Ti1 = c(init = 15, lb = 0, ub = 40))
  model$setParameter(Tm = c(init = 15, lb = 0, ub = 40))
  ##----------------------------------------------------------------
  # Set the initial value for the optimization
  model$setParameter(Ci = c(init = 1, lb = 1E-5, ub = 1E5))
  model$setParameter(Cm = c(init = 1000, lb = 1E-5, ub = 1E5))
  model$setParameter(Ria = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(Rim = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(Rin = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(Aw = c(init = 6, lb = 1E-2, ub = 7.5+4.8+5)) 
  model$setParameter(sigma1 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(sigma2 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(e1 = c(init = -1, lb = -50, ub = 10))
  
  ##----------------------------------------------------------------    

  # Run the parameter optimization

  fit = model$estimate(data,firstorder = TRUE)
  return(fit)
}

idx <- (Hour>8 & Hour < 23) # It is impossible to fit a window area for the hours without any sun, so we limit the window area estimation to the hours with sun.
bs = bs(Hour[idx],df=5,intercept=TRUE) 
lengt(seq(9,22))

bs1 <- bs2 <- bs3 <- bs4 <- bs5 <- bs6 <- numeric(dim(AllDat)[1])

# Idea: Use azimuth of the sun

bs1[idx] = bs[,1]
bs2[idx] = bs[,2]
bs3[idx] = bs[,3]
bs4[idx] = bs[,4]
bs5[idx] = bs[,5]

AllDat$bs1 = bs1
AllDat$bs2 = bs2
AllDat$bs3 = bs3
AllDat$bs4 = bs4
AllDat$bs5 = bs5

fit1 <- sdeImproved(AllDat)

summary(fit1,extended=TRUE)

pred <- predict(fit1)
mean((pred[[1]]$state$pred$Ti1 - AllDat$yTi1)^2)


sdeImproved2 <- function(data){
  # Generate a new object of class ctsm
  model = ctsm()
  # Add a system equation and thereby also a state
  # Gv in Ti: Aw/Ci*Gv or Tm: Aw/Cm*Gv
  model$addSystem(dTi1 ~  1/Ci*(1/Rin*(yTi2-Ti1) + 1/Ria*(Ta-Ti1) + 1/Rim*(Tm-Ti1)  + Ph1 + Aw*Gv)*dt + exp(sigma1)*dw1)
  model$addSystem(dTm ~  1/Cm*(1/Rim*(Ti1-Tm) + 1/Rma*(Ta-Tm))*dt + exp(sigma2)*dw2)
  # Set the names of the inputs
  model$addInput(Ta,yTi2,Gv,Ph1)
  # Set the observation equation: Ti is the state, yTi is the measured output
  model$addObs(yTi1 ~ Ti1)
  # Set the variance of the measurement error
  model$setVariance(yTi1 ~ exp(e1))
  ##----------------------------------------------------------------
  # Set the initial value (for the optimization) of the value of the state at the starting time point
  model$setParameter(Ti1 = c(init = 15, lb = 0, ub = 40))
  model$setParameter(Tm = c(init = 15, lb = 0, ub = 40))
  ##----------------------------------------------------------------
  # Set the initial value for the optimization
  model$setParameter(Ci = c(init = 1, lb = 1E-5, ub = 1E5))
  model$setParameter(Cm = c(init = 1000, lb = 1E-5, ub = 1E5))
  model$setParameter(Ria = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(Rim = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(Rin = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(Rma = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(Aw = c(init = 6, lb = 1E-2, ub = 7.5+4.8+5)) 
  model$setParameter(sigma1 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(sigma2 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(e1 = c(init = -1, lb = -50, ub = 10))
  
  ##----------------------------------------------------------------    

  # Run the parameter optimization

  fit = model$estimate(data,firstorder = TRUE)
  return(fit)
}

fit2 <- sdeImproved2(AllDat)

pred <- predict(fit2)
mean((pred[[1]]$state$pred$Ti1 - AllDat$yTi1)^2)
