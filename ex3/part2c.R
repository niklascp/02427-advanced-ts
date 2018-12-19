library(ctsmr)
library(ggplot2)
load("Exercise3.RData")

sdeMultiRoom <- function(data){
  # Generate a new object of class ctsm
  model = ctsm()
  # Add a system equation and thereby also a state
  # Gv in Ti: Aw/Ci*Gv or Tm: Aw/Cm*Gv
  model$addSystem(dT1 ~ 1/Co*(1/Ro*(Ta-T1) + 1/Rn12*(T2-T1) + 1/Rm*(Tm-T1) + Ph1 + Aw1*Gv)*dt + exp(sigma1)*dw1)
  model$addSystem(dT2 ~ 1/Ch*(1/Rh*(Ta-T2) + 1/Rn12*(T1-T2) + 1/Rn23*(T3-T2) + 1/Rm*(Tm-T2) + Ph1 + Aw2*Gv)*dt + exp(sigma2)*dw2)
  model$addSystem(dT3 ~ 1/Ch*(1/Rh*(Ta-T3) + 1/Rn23*(T2-T3) + 1/Rn34*(T4-T3) + 1/Rm*(Tm-T3) + Ph2 + Aw3*Gv)*dt + exp(sigma3)*dw3)
  model$addSystem(dT4 ~ 1/Co*(1/Ro*(Ta-T4) + 1/Rn34*(T3-T4) + 1/Rm*(Tm-T4) + Ph2 + Aw4*Gv)*dt + exp(sigma4)*dw4)
  model$addSystem(dTm ~ 1/Cm*(1/Rm*(T1-Tm) + 1/Rm*(T2-Tm) + 1/Rm*(T3-Tm) + 1/Rm*(T4-Tm))*dt + exp(sigma5)*dw5)
  # Set the names of the inputs
  model$addInput(Ta,Ph1,Ph2,Gv)
  # Set the observation equation: Ti is the state, yTi is the measured output
  model$addObs(yTi1 ~ T1)
  model$addObs(yTi2 ~ T2)
  model$addObs(yTi3 ~ T3)
  model$addObs(yTi4 ~ T4)
  # Set the variance of the measurement error
  model$setVariance(yTi1 ~ exp(e1))
  model$setVariance(yTi2 ~ exp(e2))
  model$setVariance(yTi3 ~ exp(e3))
  model$setVariance(yTi4 ~ exp(e4))
  ##----------------------------------------------------------------
  # Set the initial value (for the optimization) of the value of the state at the starting time point
  model$setParameter(T1 = c(init = 15, lb = 0, ub = 40))
  model$setParameter(T2 = c(init = 15, lb = 0, ub = 40))
  model$setParameter(T3 = c(init = 15, lb = 0, ub = 40))
  model$setParameter(T4 = c(init = 15, lb = 0, ub = 40))
  model$setParameter(Tm = c(init = 15, lb = 0, ub = 40))
  ##----------------------------------------------------------------
  # Set the initial value for the optimization
  model$setParameter(Co = c(init = 1, lb = 1E-5, ub = 1E5))
  model$setParameter(Ch = c(init = 1, lb = 1E-5, ub = 1E5))
  model$setParameter(Cm = c(init = 1000, lb = 1E-5, ub = 1E5))
  model$setParameter(Ro = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(Rh = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(Rm = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(Rn12 = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(Rn23 = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(Rn34 = c(init = 20, lb = 1E-4, ub = 1E5))
  model$setParameter(Aw1 = c(init = 6, lb = 1E-2, ub = 7.5+4.8+5)) 
  model$setParameter(Aw2 = c(init = 6, lb = 1E-2, ub = 7.5+4.8+5)) 
  model$setParameter(Aw3 = c(init = 6, lb = 1E-2, ub = 7.5+4.8+5)) 
  model$setParameter(Aw4 = c(init = 6, lb = 1E-2, ub = 7.5+4.8+5)) 
  model$setParameter(sigma1 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(sigma2 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(sigma3 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(sigma4 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(sigma5 = c(init = 1, lb = -30, ub = 10))
  model$setParameter(e1 = c(init = -1, lb = -5, ub = 5))
  model$setParameter(e2 = c(init = -1, lb = -5, ub = 5))
  model$setParameter(e3 = c(init = -1, lb = -5, ub = 5))
  model$setParameter(e4 = c(init = -1, lb = -5, ub = 5))

  ##----------------------------------------------------------------    

  # Run the parameter optimization

  fit = model$estimate(data,firstorder = TRUE)
  return(fit)
}

fit_multi_room <- sdeMultiRoom(AllDat)

summary(fit_multi_room,extended=FALSE)

pred_multi_room <- predict(fit_multi_room)
mean((pred_multi_room[[1]]$state$pred$T1 - AllDat$yTi1)^2)
mean((pred_multi_room[[1]]$state$pred$T2 - AllDat$yTi2)^2)
mean((pred_multi_room[[1]]$state$pred$T3 - AllDat$yTi3)^2)
mean((pred_multi_room[[1]]$state$pred$T4 - AllDat$yTi4)^2)
