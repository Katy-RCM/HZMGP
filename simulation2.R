#-------------------------------#
# Fitting HZMGP data with .stan #
#-------------------------------#

rm(list=ls(all=TRUE))

require(rstan)

# Rutas
data_address   <- "data_HZMGP.RData"
stan_address   <- "zmpg_exp.stan"  # change with 'zmp_exp.stan' or 'sm_exp.stan'
output_address <- " "

load(data_address)

B <- 100

for (i in 1:B) {
  
  cat("Fit:", i, "started", format(Sys.time(), "%H:%M"), "\n") 
  
  data <- data_gen[[i]]
  
  t     <- data$t
  delta <- data$cens
  a     <- data$x_1
  
  x1    <- model.matrix(~ a, data = data)
  x2    <- model.matrix(~ a, data = data)
  
  set.seed(123+i)
  
  fit_hzmgp <- stan( file = stan_address,
                         data = list(n = length(t), Nbetas_omega = ncol(x1), Nbetas_mu = ncol(x2),
                                     t = t, delta = delta, x1 = x1, x2 = x2),
                         warmup = 300, iter = 1000,
                         chains = 3, seed = 123, cores = getOption("mc.cores", 3))
  
  sim <- extract(fit_hzmgp)
  rmv <- c("omega", "mu", "lp__")
  sim_final <- sim[!(names(sim) %in% rmv)]
  
  output_name <- paste0("sim_",i, ".RData")
  output_address_i <- file.path(output_address, output_name)
  save(sim_final, file = output_address_i)
  
} 


#------------------------------#
# Fitting HZMP data with .stan #
#------------------------------#

rm(list=ls(all=TRUE))

require(rstan)

# Rutas
data_address   <- "data_HZMP.RData"
stan_address   <- "zmp_exp.stan"  # change with 'zmpg_exp.stan' or 'sm_exp.stan'
output_address <- " "

load(data_address)

B <- 100

for (i in 1:B) {
  
  cat("Fit:", i, "started", format(Sys.time(), "%H:%M"), "\n") 
  
  data <- data_hzmp[[i]]
  
  t     <- data$t
  delta <- data$cens
  a     <- data$x_1
  
  x1    <- model.matrix(~ a, data = data)
  x2    <- model.matrix(~ a, data = data)
  
  set.seed(123+i)
  
  fit_hzmp <- stan( file = ruta_stan,
                    data = list(n = length(t), Nbetas_omega = ncol(x1), Nbetas_mu = ncol(x2),
                                t = t, delta = delta, x1 = x1, x2 = x2),
                    warmup = 300, iter = 1000,
                    chains = 3, seed = 123, cores = getOption("mc.cores", 3))
  
  sim       <- extract(fit_hzmp)
  rmv       <- c("omega", "mu", "lp__")
  sim_final <- sim[!(names(sim) %in% rmv)]
  
  output_name <- paste0("sim_final_",i, ".RData")
  output_address_i <- file.path(ruta_salida, output_name)
  save(sim_final, file = output_address_i)
  
} 


#----------------------------#
# Fitting SM data with .stan #
#----------------------------#

rm(list=ls(all=TRUE))

require(rstan)

# Rutas
data_address   <- "data_SM.RData"
stan_address   <- "sm_exp.stan"  # change with 'zmp_exp.stan' or 'zmpg_exp'
output_address <- " "

load(data_address)

B <- 100

for (i in 1:B) {
  
  cat("Fit:", i, "started", format(Sys.time(), "%H:%M"), "\n") 
  
  data <- datos_gera_mp_x1[[i]]
  
  t     <- data$t
  delta <- data$delta
  a     <- data$x_1
  
  XC    <- model.matrix(~ a, data = data)
  XU    <- XC[,-1, drop = FALSE]
  
  set.seed(123+i)
  
  fit_sm <- stan(file = stan_address,
                 data = list(n = length(t), Nbetas = ncol(XC),
                             t = t, delta = delta, XC = XC, XU = XU),
                 warmup = 300, iter = 1000, chains = 3, 
                 seed = 123, cores = getOption("mc.cores", 3))
  
  sim_sm   <- extract(fit_sm)
  excluir  <- c("lp__")
  sim_sm_1 <- sim_sm[!(names(sim_sm) %in% excluir)]
  
  output_name <- paste0("sim_sm_",i, ".RData")
  output_address_i <- file.path(output_address, output_name)
  save(sim_sm_1, file = output_address_i)
  
} 