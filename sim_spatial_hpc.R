args = commandArgs(trailingOnly = TRUE)
i = as.integer(Sys.getenv("LSB_JOBINDEX"))
# set the simulation parameters
NN = as.integer(args[1])
path = "/gpfs_common/share01/research1/wzhao24/written"
code_path = paste0(path, "/code")
res_path = paste0(path, "/resultsNN_",NN)

library(GeoModels)
library(fields)
library(spam)

# specify spatial locations
set.seed(997)

# specify spatial locations
x = runif(NN, 0, 1)
y = runif(NN, 0, 1)
coords = cbind(x, y)

# correlation parameters
corrmodel = "GenWend_Matern"
GeoModels::CorrParam(corrmodel)
mu = 3
scale = 0.05
power2 = 1/mu
smooth = 0

# regression parameters
model = "Gaussian"  # switch to "logGaussian" for non-Gaussian random field
NuisParam(model)
mean = 5
nugget = 0
sill = 1

param = list(nugget = nugget, mean = mean, scale = scale, sill = sill, power2 = power2, smooth = smooth)
start = list(mean = mean, scale = scale, sill = sill, power2 = power2)
fixed = list(nugget = nugget, smooth = smooth)
lower = list(mean = -1e4, scale = 1e-4, sill = 1e-4, power2 = 1e-4)
upper = list(mean = 1e4, scale = 1e4, sill = 1e4, power2 = 1/1.5)
truth = c(mean = mean, power2 = power2, scale = scale, sill = sill)

optimizer = "nlminb"
  
set.seed(i)
cat("NN =", NN, ", Simulation =", i, "\n")
    
# simulate the spatial Gaussian random field
data = GeoSim(coordx = coords, corrmodel = corrmodel, 
              sparse = TRUE, model = model, param = param)$data
    
# Full likelihood fit

GeoFit2(data = data, coordx = coords, corrmodel = corrmodel,
        model = model, varest = TRUE,
        optimizer = optimizer, 
        lower = lower, upper = upper,
        likelihood = "Full", type = "Standard",
        start = start, fixed = fixed)

t0 = Sys.time()
fit_F = tryCatch(GeoFit2(data = data, coordx = coords, corrmodel = corrmodel,
                        model = model, varest = TRUE,
                        optimizer = optimizer, 
                        lower = lower, upper = upper,
                        likelihood = "Full", type = "Standard",
                        sensitivity = TRUE, start = start, fixed = fixed),
                 error = function(e){
                   GeoFit(data = data, coordx = coords, corrmodel = corrmodel,
                          model = model, varest = TRUE,
                          optimizer = optimizer, 
                          lower = lower, upper = upper,
                          likelihood = "Full", type = "Standard",
                          sensitivity = TRUE, start = start, fixed = fixed)
                 })
if(fit_F$convergence!="Successful"){
  cat("Full likelihood fit did not converge\n")
}
t_F = difftime(Sys.time(), t0, units = "secs")
est_F = unlist(fit_F$param)
se_F = fit_F$stderr
cat("Full likelihood fit done\n")
    
# Marginal likelihood fit
t0 = Sys.time()
fit_M = GeoFit2(data = data, coordx = coords, corrmodel = corrmodel,
                model = model, varest = TRUE,
                optimizer = optimizer, 
                lower = lower, upper = upper,
                likelihood = "Marginal", type = "Pairwise", neighb = 3,
                sensitivity = TRUE, start = start, fixed = fixed)
if(fit_M$convergence!="Successful"){
  cat("Marginal likelihood fit did not converge\n")
}
t_M = difftime(Sys.time(), t0, units = "secs")
    
GeoVarestbootstrapRepeat <- function(fit, ...) {
  flag <- FALSE
  cat("flag is ", flag, "\n")
  rst <- NULL
  while (!flag) {
    tryCatch({
      # sample(1, sise = 5)
      rst <- GeoVarestbootstrap(fit, ...)
      flag <- TRUE
      cat("flag is ", flag, "\n")
    }, error = function(e) {
      e
      # cat("flag is ", flag, "\n")
      # cat('Retrying...\n')
    })
  }
  return(rst)
}
    
vr_M = GeoVarestbootstrapRepeat(fit_M, ncores = 1, parallel = TRUE)
est_M = unlist(fit_M$param)
se_M = vr_M$stderr
cat("Marginal likelihood fit done\n")
    
# Conditional likelihood fit
t0 = Sys.time()
fit_C = GeoFit2(data = data, coordx = coords, corrmodel = corrmodel,
                model = model, varest = TRUE,
                optimizer = optimizer, 
                lower = lower, upper = upper,
                likelihood = "Conditional", type = "Pairwise", neighb = 3,
                sensitivity = TRUE, start = start, fixed = fixed)
if(fit_C$convergence!="Successful"){
  cat("Conditional likelihood fit did not converge\n")
}
t_C = difftime(Sys.time(), t0, units = "secs")
vr_C = GeoVarestbootstrapRepeat(fit_C, ncores = 1, parallel = TRUE)
est_C = unlist(fit_C$param)
se_C = vr_C$stderr
cat("Conditional likelihood fit done\n")
    
# tapered likelihood 
t0 = Sys.time()
fit_T = tryCatch(GeoFit2(data = data, coordx = coords, corrmodel = corrmodel,
                         model = model, varest = TRUE,
                         optimizer = optimizer, 
                         lower = lower, upper = upper,
                         likelihood = "Full", type = "Standard",
                         taper = 'Tapering', sparse = TRUE,
                         sensitivity = TRUE, start = start, fixed = fixed),
                 error = function(e){
                   GeoFit(data = data, coordx = coords, corrmodel = corrmodel,
                           model = model, varest = TRUE,
                           optimizer = optimizer, 
                           lower = lower, upper = upper,
                           likelihood = "Full", type = "Standard",
                           taper = 'Tapering', sparse = TRUE,
                           sensitivity = TRUE, start = start, fixed = fixed)
                 })
if(fit_T$convergence!="Successful"){
  cat("Tapered likelihood fit did not converge\n")
}
t_T = difftime(Sys.time(), t0, units = "secs")
est_T = unlist(fit_T$param)
se_T = fit_T$stderr
cat("Tapering likelihood fit done\n")

save.image(file = paste0(res_path,"/Sim_",i, "_NN_", NN,".RData"))

