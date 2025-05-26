require(survival)
require(lamW)
require(dplyr)
require(purrr) 

# HZMGP - survival function
st <- function(t,x_1,x_2,phi,lambda,gamma,b00,b10,b20,b01,b11,b21){
  omega <- (exp(b00+b10*x_1+b20*x_2)/(1+exp(b00+b10*x_1+b20*x_2)))
  mu    <- exp(b01+b11*x_1+b21*x_2)
  aux0 <- (mu*phi)/(1+(mu*phi))
  aux1 <- exp(-lambda*(t^gamma)) 
  aux2 <- -aux0*aux1*exp(-aux0)
  aux3 <- lambertW0(aux2)
  aux4 <- -(1/phi)
  p    <- (omega)/(1-exp((-mu)/(1+(mu*phi))))
  
  st0 <- 1-p+(p*exp(aux4*(aux3+aux0)))
  return(st0)
}

# Fixed parameters
phi    <- 0.13
lambda <- 0.10
gamma  <- 1.08
b00    <- 1.5 # omega
b10    <- 0.6  
b20    <- 3.5  
b01    <- 1.4 # mu
b11    <- 0.1  
b21    <- 1.2  

# Generating data
B <- 100   # Replicates
n <- 20000 # Sample size
data_gen <- list()

for(b in 1:B){
  cat("Data:", b, "/", B, "\n")
  set.seed(123+b)
  # Generate covariates
  x_1 <- rnorm(n = n, 0, 1)                     # Continuous var
  x_2 <- rbinom(n = n, size = 1, prob = 0.8)    # Dummy var
  
  data <- tibble(x_1 = x_1, x_2 = x_2) %>% 
    mutate(phi = phi, lambda = lambda, gamma = gamma, b00 = b00, b10 = b10, b20 = b20, b01 = b01, b11 = b11, b21 = b21) %>% 
    mutate(t = purrr::pmap_dbl(.l = list(phi, lambda, gamma, b00, b10, b20, b01, b11, b21, x_1, x_2), .f = function(p, l, g, b_0, b_1, b_2, b_3, b_4, b_5, x, z) {
      u <- runif(1)
      yAux <- tryCatch(
        uniroot(f = function(y) st(y, x, z, p, l, g, b_0, b_1, b_2, b_3, b_4, b_5) - u, lower = 0, upper = 20),
        error = function(e) NULL
      )
      if (is.null(yAux)) Inf else yAux$root
    })) %>% 
    mutate(cens = ifelse(is.infinite(t), 0, 1)) %>% 
    mutate(t = ifelse(is.infinite(t), max(t[is.finite(t)]), t))
  
  # Final data frame for this dataset
  t_final     <- (data$t) + 0.001 
  delta_final <- data$cens      
  x_final     <- data$x_1       # Continuous var
  z_final     <- data$x_2       # Dummy var
  
  data_gen[[b]] <- data.frame(t = t_final, delta = delta_final, x = x_final, z = z_final)
} 

# Simulation process
data_address  <- "data_sim2.RData"
stan_address  <- "zmpg_model.stan"
output_address <- ""

load(data_address)
for(b in 1:B){
  cat("Fit:", b)
  # Data extraction
  data  <- data_gen[[b]]
  tt    <- data$t
  delta <- data$delta
  x     <- data$x
  z     <- data$z
  x1    <- model.matrix(~ x + z, data = data)
  x2    <- model.matrix(~ x + z, data = data)
  
  # Model fit
  set.seed(123+b)
  fit <- stan( file = stan_address,
                   data = list(n = length(t), Nbetas_omega = ncol(x1), Nbetas_mu = ncol(x2),
                               t = tt, delta = delta, x1 = x1, x2 = x2),
                   warmup = 300, iter = 1000,
                   chains = 3, seed = 123, cores = getOption("mc.cores", 3))
  
  sim  <- extract(fit)
  rmv  <- c("omega", "mu", "lp__")
  sim_final <- sim[!(names(sim) %in% rmv)]
  
  output_name <- paste0("sim2_", b, ".RData")
  output_address_b <- file.path(output_address, output_name)
  save(sim_final, file = output_address_b)
} 

# Obtaining mean, SD and 95% CP
N_ <- 20000
# Vector of fixed values
a0 <- c(phi=0.13, lambda=0.10, gamma=1.08, b00=1.5, b10=0.6, b20=3.5, b01=1.4, b11=0.1, b21=1.2)
pmean <- matrix(NA, ncol=length(a0), nrow=B)  
ci95   <- matrix(NA, ncol=2*length(a0), nrow=B)   
cp95   <- matrix(NA, ncol=length(a0), nrow=B)   

for(b in 1:B){
  load(paste0("sim2_", b, ".RData")) 
  pmean[b, ] <- c(
    round(mean(sim_final$phi), 3),
    round(mean(sim_final$lambda), 3),
    round(mean(sim_final$gamma), 3),
    round(mean(sim_final$beta1[,1]), 3),
    round(mean(sim_final$beta1[,2]), 3),
    round(mean(sim_final$beta1[,3]), 3),
    round(mean(sim_final$beta2[,1]), 3),
    round(mean(sim_final$beta2[,2]), 3),
    round(mean(sim_final$beta2[,3]), 3))
  
  ci95[b, ] <- c(
    round(quantile(sim_final$phi,       c(0.025, 0.975)), 3),
    round(quantile(sim_final$lambda,    c(0.025, 0.975)), 3),
    round(quantile(sim_final$gamma,      c(0.025, 0.975)), 3),
    round(quantile(sim_final$beta1[,1], c(0.025, 0.975)), 3),
    round(quantile(sim_final$beta1[,2], c(0.025, 0.975)), 3),
    round(quantile(sim_final$beta1[,3], c(0.025, 0.975)), 3),
    round(quantile(sim_final$beta2[,1], c(0.025, 0.975)), 3),
    round(quantile(sim_final$beta2[,2], c(0.025, 0.975)), 3),
    round(quantile(sim_final$beta2[,3], c(0.025, 0.975)), 3))
  
  for(j in 1:9){
    if(ci95[b, (2*j-1)] <= a0[j] && ci95[b, (2*j)] >= a0[j]){
      cp95[b, j] <- 1
    }else{
      cp95[b, j] <- 0
    }
  }
  # print(b)
}

output1 <- list(pmean = pmean, cp95 = cp95)
output_final <- data.frame(
  Parameters = c("phi","lambda", "gamma", "b00", "b10", "b20", "b01", "b11", "b21"),
  Mean = format(round(colMeans(pmean, na.rm = TRUE), 3), scientific = FALSE, nsmall = 3),  
  SD   = format(round(apply(pmean, 2, sd, na.rm = TRUE), 3), scientific = FALSE, nsmall = 3),  
  CP95 = format(round(colMeans(cp95, na.rm = TRUE), 3), scientific = FALSE, nsmall = 3)  
) 

# Classification
ome <- lisup <- list()
propor <- propor2 <- list()
clasi  <- rep(NA, B)
                               
for(b in 1:B){
  data <- data_gen[[b]]
  x1 <- data$x 
  x2 <- data$z 
  load(paste0("sim2_", b, ".RData"))
  phi <- (sim_final$phi)
  b0 <- (sim_final$beta1[,1]) # omega
  b1 <- (sim_final$beta1[,2])  
  b2 <- (sim_final$beta1[,3])
  b_0 <- (sim_final$beta2[,1]) # mu
  b_1 <- (sim_final$beta2[,2]) 
  b_2 <- (sim_final$beta2[,3]) 
  
  # Calculating omega and mu
  mu         <- exp(b_0+b_1*x1[b]+b_2*x2[b])
  lisup[[b]] <- (1 - exp(-mu/(1+mu*phi)))
  b_ome      <- (b0+b1*x1[b]+b2*x2[b])
  ome[[b]]   <- exp(b_ome)/(1+exp(b_ome))
  
  # Comparation
  compa      <-  ome[[b]] < lisup[[b]] 
  compa2     <-  ome[[b]] < 1
  propor[[b]]  <- mean(compa)
  propor2[[b]] <- mean(compa2)
  
  prop  <- propor[[b]]
  prop2 <- propor2[[b]]
  
  # Classification
  if(prop2 >= 0.95 & prop2 <= 1){
    if(prop >= 0.95 & prop <= 1){
      clasi[b] <- "Inflated"   
    } else if(prop >= 0.05 & prop < 0.95){
      clasi[b] <- "Tradicional" 
    } else if(prop >= 0.00 & prop < 0.05){
      clasi[b] <- "Deflated"   
    }
  }else {clasi[b] <- "Truncated"}
  
  print(b)
}

prop.table(table(clasi))
