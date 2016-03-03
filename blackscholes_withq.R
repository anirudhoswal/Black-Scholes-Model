black_scholes <- function(S, k, Tm, r, sigma,q){
  values <- c(1)
  d1 <- (log(S/k) + (r+-q+(sigma^2)/2)*(Tm))/(sigma*sqrt(Tm))
  'd2 <- (log(S/k) + (r-(sigma^2)/2)*(Tm))/(sigma*sqrt(Tm))'
  d2 <- d1 - (sigma*sqrt(Tm))
  values[1] <- S*exp(-q*Tm)*pnorm(d1) - k*exp(-r*Tm)*pnorm(d2)
  values[2] <- k*exp(-r*Tm)*pnorm(-d2)- S*exp(-q*Tm)*pnorm(-d1)
  
  values
}

