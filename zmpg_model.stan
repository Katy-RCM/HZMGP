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
  real theta1; // theta1 = log(phi)
  real theta2; // theta2 = log(lambda)
  real theta3; // theta3 = log(gama)
  vector[Nbetas_omega] beta1;
  vector[Nbetas_mu] beta2;
}

transformed parameters {
  
  real phi = exp(theta1);    
  real lambda = exp(theta2); 
  real gama = exp(theta3);   
  
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
  vector[n] h0;
  vector[n] b;
  vector[n] lamfc;
  vector[n] st;
  vector[n] ht;
  vector[n] p;
  
  
  // Likelihood using hazard function
  for (i in 1:n){
    
    s0[i]    = exp(-lambda * t[i]^gama);    // Baseline survival function
    h0[i]    = gama*lambda*(t[i]^(gama-1)); // Baseline hazard function
    b[i]     = mu[i]*phi/(1+(mu[i]*phi));
    lamfc[i] = lambert_w0(-b[i]*s0[i]*exp(-b[i]));
    p[i]     = omega[i]/(1-exp(-mu[i]/(1+(mu[i])*phi)));
    st[i]    = 1-p[i]+p[i]*exp((-1/phi)*(lamfc[i] + b[i]));  // Survival function
    ht[i]    = (-p[i]*h0[i]*exp((-1/phi)*(lamfc[i]+ b[i]))/(phi*st[i]))*(lamfc[i]/(1+lamfc[i])); // Hazard function
    target += delta[i] * log(ht[i]) + log(st[i]);  
  }
    
  // Prior distributions for parameters

  target += normal_lpdf(theta1 | 0, 10);  // phi
  target += normal_lpdf(theta2 | 0, 10);  // lambda
  target += normal_lpdf(theta3 | 0, 10);  // gama
  target += normal_lpdf(beta1 | 0, 10);   // Coefficients omega
  target += normal_lpdf(beta2 | 0, 10);   // Coefficients mu
}