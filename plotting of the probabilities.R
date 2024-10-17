source("Monte_carlo_calculation.R")
source("reserve_calculation_2d.R")
source("two_dimensional_estimation.R")
load("simulated_Data_Aalen_Johansen_1000.Rdata")
#load("simulated_Data_Aalen_Johansen_5000.Rdata")
load("Monte_Carlo_Data.Rdata")
#Setting:
#sim is the censored simulated data for Aalen-Johansen estimation
#sim_mc is the simulated non-censored data for the Monte-Carlo estimation
#initial payment
b_0 <- 100000
#Premium payment
premium <- 10000
#interest function
interest <- 3
interest_f <- function(t1, t2) {
  return(1 / (log(100 + interest) - log(100)) *
           (-(1 + interest / 100)^(-t2) + (1 + interest / 100)^(-t1)) *
           (1 + interest / 100)^(40))
}
#pension (normally calculated via the technical basis)
p_11_technical <- function(x) {
  exp(-(10^((19 * x) / (500) - 159 / (125))) / (38 * log(10))
      - x / (2000)
      + 1 / (38 * 10^(159 / (125)) * log(10)))
}
pension <- (b_0 + premium * integrate(p_11_technical, 40, 65)$value) /
  (integrate(p_11_technical, 65, Inf)$value)
## rho (normally calculated via the technical basis)
V_1_plus <- function(t) {
  if (t > 65) {
    result <- integrate(p_11_technical, t, Inf)$value *
      pension
  } else {
    result <- integrate(p_11_technical, 65, Inf)$value *
      pension
  }
  return(result)
}
V_1_minus <- function(t) {
  if (t > 65) {
    result <- 0
  } else if (t > 40) {
    result <- integrate(p_11_technical, t, 65)$value * premium
  } else {
    result <- b_0 + integrate(p_11_technical, 40, 65)$value * premium
  }
  return(result)
}
rho <- function(t) {
  if (t >= 65) {
    return(1)
  } else if (t < 40) {
    return(0)
  } else if (V_1_plus(t)>10^(-12)){
    return((V_1_plus(t) - V_1_minus(t)) / V_1_plus(t))
  } else{return(0)}
}
#estimation 2-D
estimation<-aalen_johansen_2d(sim)
#reserve calculation
reserve_2<-reserve_2d(estimation,b_0,premium,
                      interest,interest_f,pension,rho,V_1_plus,V_1_minus)
monte_carl<-monte_carlo_2d(sim_mc)

#graphical analysis
index_2d <- which(estimation$ordered_times > 100)[1]
index_mc <- which(monte_carl$ordered_times > 100)[1]
if (is.na(index_2d)){index_2d=length(estimation$ordered_times)}
if (is.na(index_mc)){index_mc=length(monte_carl$ordered_times)}
plot(estimation$ordered_times[1:index_2d],reserve_2[1:index_2d],col = "red",type = "l", lty = 2,xlab="",ylab="")
lines(monte_carl$ordered_times[1:index_mc],monte_carl$reserve_cont[1:index_mc],col= "blue")
legend("bottomright",                
       legend = c("Monte Carlo estimation", "as-if Markov estimation"),  
       col = c("blue", "red"),   
       lty = c(1, 2),             
       lwd = 1,
       cex = 0.7) 



