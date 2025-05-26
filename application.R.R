library(lamW) 
library(survival)
library(xtable)
library(rstan)

# Data
data  <- read.table("lung_data.csv", header = TRUE, sep=",")
data1 <- na.omit(data)

# Survival information
tt    <- data1$time
delta <- data1$cens

# Covariates
x_1  <- as.numeric(scale(data1$age))
x_2  <- data1$gender
x_3  <- ifelse(data1$clinical_stage %in% c(1, 2), 0, 1)
x_4  <- data1$surgery
x_5  <- data1$radio
x_6  <- data1$chemo

x1 <- model.matrix(~ x_1+x_2+x_3+x_4+x_5+x_6, data = data1) #omega
x2 <- model.matrix(~ x_1+x_2+x_3+x_4+x_5+x_6, data = data1) #mu

set.seed(145)

# fitted model
m1 <- stan(file   = 'zmpg_model.stan', 
                 data   = list(n = length(t), Nbetas_omega = ncol(x1),  Nbetas_mu = ncol(x2),
                               t = tt, delta = delta, x1 = x1, x2 = x2),
                 warmup = 300, iter   = 1000,
                 init   = 1,   chains = 3,
                 seed   = 123, cores  = getOption("mc.cores",3))

print(m1)

# Diagnostics

pars <- c( "phi", "lambda", "gama","beta1", "beta2") 

plot(m1, plotfun = "trace", pars = pars, inc_warmup = TRUE)
plot(m1, plotfun = "hist", pars = pars, col = "lightblue",lty = 1)

# Classification of zero structures for each patient

post1 <- extract(m1)

phi <- (post1$phi)

b0 <- (post1$beta1[,1]) #omega
b1 <- (post1$beta1[,2]) 
b2 <- (post1$beta1[,3]) 
b3 <- (post1$beta1[,4])
b4 <- (post1$beta1[,5])
b5 <- (post1$beta1[,6])
b6 <- (post1$beta1[,7])

b_0 <- (post1$beta2[,1]) #mu
b_1 <- (post1$beta2[,2]) 
b_2 <- (post1$beta2[,3]) 
b_3 <- (post1$beta2[,4])
b_4 <- (post1$beta2[,5])
b_5 <- (post1$beta2[,6])
b_6 <- (post1$beta2[,7])

lisup <- list()
ome   <- list()

# Calculating omega and mu
for(i in 1:n){
  mu         <- exp(b_0+b_1*x1[i]+b_2*x2[i]+b_3*x3[i]+b_4*x4[i]+b_5*x5[i]+b_6*x6[i])
  b_ome      <- (b0+b1*x1[i]+b2*x2[i]+b3*x3[i]+b4*x4[i]+b5*x5[i]+b6*x6[i])
  lisup[[i]] <- (1 - exp(-mu/(1+mu*phi)))
  ome[[i]]   <- exp(b_ome)/(1+exp(b_ome))
}

propor <- propor2 <- rep(NA,n)

# comparating omega and threshold
for(i in 1:n){
  compa <-  ome[[i]] < lisup[[i]] 
  compa2 <-  ome[[i]] < 1
  propor[i] <- mean(compa)
  propor2[i] <- mean(compa2)
}

clasi<- rep(NA, n)

# classification of zero-modification
for(i in 1:n){
  prop <- propor[i]
  prop2 <- propor2[i]
  
  if(prop2 >= 0.95 & prop2 <= 1){
    if(prop >= 0.95 & prop <= 1){
      clasi[i] <- "Inflated"   
    } else if(prop >= 0.05 & prop < 0.95){
      clasi[i] <- "Tradicional" 
    } else if(prop >= 0.00 & prop < 0.05){
      clasi[i] <- "Deflated"   
    }
  }else {clasi[i] <- "Truncated"}
  
}

prop.table(table(clasi2))