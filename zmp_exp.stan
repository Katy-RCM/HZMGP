data{
  int n; 
  int Nbetas_omega;
  int Nbetas_mu;
  vector[n] t;
  vector[n] delta;
  matrix[n, Nbetas_omega] x1;
  matrix[n, Nbetas_mu] x2;
}

parameters{
  real theta2; // theta2 = log(lambda)
  vector[Nbetas_omega] beta1;
  vector[Nbetas_mu] beta2;
}

transformed parameters {
  
  real lambda = exp(theta2); 
  
  vector[n] omega; 

  // X' * beta1
  for (i in 1:n) {
    omega[i] = (exp(x1[i,]*beta1))/(1+exp(x1[i,]*beta1)); 
  }
  
  vector[n] mu; 

  // X' * beta2
  for (i in 1:n) {
    mu[i] = exp(x2[i,]*beta2); 
  }
  
}

model {
  vector[n] s0;
  vector[n] f0;
  vector[n] aux;
  vector[n] logst;
  vector[n] loght;
  
  // ZMP frailty model
  for(i in 1:n){
     // baseline functions
     s0[i]  = exp(-lambda * t[i]); // baseline functions 
     f0[i]  = exp(-lambda * t[i])*lambda; // baseline functions
     aux[i] = 1-exp(-mu[i]); //denom de la reparametrizacion
     
     //survival function
     logst[i] = log(1 - (omega[i] / aux[i]) + (omega[i] / aux[i]) * exp(-mu[i] * (1 - s0[i])));
     //hazard function
     loght[i] = log(omega[i]) + log(mu[i]) + log(f0[i]) - mu[i] * (1 - s0[i]) - log(aux[i]) - logst[i];
     
     target += delta[i] * loght[i] + logst[i];
  }

  // Prior distributions for parameters

  target += normal_lpdf(theta2 | 0, 10);  // lambda
  target += normal_lpdf(beta1 | 0, 10);   // Coefficients omega
  target += normal_lpdf(beta2 | 0, 10);   // Coefficients mu
}
