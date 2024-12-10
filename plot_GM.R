library(GeoModels)
library(hypergeo)

# Generalized Wendland function parameters and definition
generalized_wendland = function(r, nu = 1, g = 1, psi = 1) {
  mu = (g + 2 * nu + 1) / 2
  delta_vgp = psi * (gamma(mu) / gamma(g))^(1 / (1 + 2 * nu))
  if(nu == 0){
    K = ( gamma(2 * nu + g + 1)) / 
      (gamma(nu + g + 1) * 2^(g+ 1))
  }else{
    K = (gamma(nu) * gamma(2 * nu + g + 1)) / 
      (gamma(2 * nu) * gamma(nu + g + 1) * 2^(g+ 1))
  }
  
  
  wendland_value = numeric(length(r))
  
  for (i in seq_along(r)) {
    if (r[i] <= delta_vgp) {
      rs = (r[i] / delta_vgp)^2
      wendland_value[i] = K * (1 - rs)^(nu + g) * 
        hypergeo(g/2, (g + 1)/2,
                 nu + g + 1, 1 - rs)
    } else {
      wendland_value[i] = 0
    }
  }
  return(wendland_value)
}

# Plotting the Generalized Wendland function
distance_seq = seq(0, 2, length.out = 500)

cov_values_GenWend = generalized_wendland(distance_seq)
cov_values_exp = exp(-3*distance_seq)

plot(distance_seq, cov_values_GenWend, type = "l", 
     xlab = "Distance", ylab = expression(rho(h)), 
     lty = 1,
     col = "steelblue")
lines(distance_seq, cov_values_exp, type = "l", 
     xlab = "Distance", ylab = "Covariance", 
     lty = 2,
     col = "pink", lwd = 2)
legend('topright',
       legend = c("Generalized Wendland",
                  "Exponential"),
       cex = 0.8, lwd = 2,
       col = c('steelblue', 'pink'), 
       lty = c(1, 2))
