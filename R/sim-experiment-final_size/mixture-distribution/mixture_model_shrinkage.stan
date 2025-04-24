## adapted from: https://osf.io/bcr4u

data {
	int N;
	vector[N] y;
	vector[N] x;
	
	int<lower=0> N_new;
  vector[N_new] x_new;
}

parameters {
	real alpha_mean;       
	vector[2] alpha_tau;
	real<lower=0> alpha_var;
	real beta_mean;        
	vector[2] beta_tau;
	real<lower=0> beta_var;
	real lambda_mean;      
	vector[2] lambda_tau;
	real<lower=0> lambda_var;
	real<lower=0, upper=1> p;         // mixing proportions
	real<lower=0> sigma_pos;
	real<lower=0> sigma_neg;
}

transformed parameters {
  real lambda_pos;      // box-cox hyper parameter
	real lambda_neg;
  real alpha_pos;       // intercepts
	real alpha_neg;
	real beta_pos;        // slopes
	real beta_neg;
  
  lambda_pos = lambda_mean + lambda_tau[1];
  lambda_neg = lambda_mean + lambda_tau[2];
  alpha_pos = alpha_mean + alpha_tau[1];
  alpha_neg = alpha_mean + alpha_tau[2];
  beta_pos = beta_mean + beta_tau[1];
  beta_neg = beta_mean + beta_tau[2];
  
}

model {
  lambda_mean ~ normal(0,2);
  lambda_var ~ cauchy(0, 1);//lambda_scale);
  // lambda_scale ~ cauchy(0, 1);
  alpha_mean ~ normal(0,2);
  alpha_var ~ cauchy(0, 1);
  beta_mean ~ normal(0,2);
  beta_var ~ cauchy(0, 1);
  for(j in 1:2){
    lambda_tau[j] ~ normal(0, lambda_var);
    alpha_tau[j] ~ normal(0, alpha_var);
    beta_tau[j] ~ normal(0, beta_var);
  }
  sigma_pos ~ cauchy(0,1);
  sigma_neg ~ cauchy(0,1);
  
  for(i in 1:N){
    if(y[i] >= 0){ // positive y values
      if(lambda_pos == 0){
        log(y[i]) ~ normal(alpha_pos + x[i]*beta_pos, sigma_pos);
        target += -log(y[i]) + log(p);
      }
      else {
        (y[i]^lambda_pos - 1)/lambda_pos ~ normal(alpha_pos + x[i]*beta_pos, sigma_pos); // transform y using box-cox
        target += (lambda_pos - 1)*log(y[i]) + log(p);
      }
    }
    else{  // negative y values
      if(lambda_neg == 0){
        log(-1*y[i]) ~ normal(alpha_neg + x[i]*beta_neg, sigma_neg);
        target += -log(-1*y[i]) + log(1-p);
      } 
      else {
        ((-1*y[i])^lambda_neg - 1)/lambda_neg ~ normal(alpha_neg + x[i]*beta_neg, sigma_neg); // transform y using box-cox
        target += (lambda_neg - 1)*log(-1*y[i]) + log(1-p);
      }
    }
  }
}

generated quantities { // for prediction intervals
  vector[N_new] y_new;
  vector[N_new] y_new_nominal;
  vector<lower=0,upper=1>[N_new] y_exclude;
  real<lower=0,upper=1> pos_flag; 
  
  for (n in 1:N_new){
    pos_flag = bernoulli_rng(p);
    y_exclude[n] = 0;
    if(pos_flag){
      y_new[n] = normal_rng(alpha_pos + x_new[n] * beta_pos, sigma_pos);
      if(lambda_pos == 0){
        y_new_nominal[n] = exp(y_new[n]);
      }
      else{
        if((lambda_pos * y_new[n] + 1) < 0){
          y_new_nominal[n] = 0;
          y_exclude[n] = 1;
        }
        else{
          y_new_nominal[n] = (lambda_pos * y_new[n] + 1)^(1/lambda_pos);
        }
      }
    }
    else{
      y_new[n] = normal_rng(alpha_neg + x_new[n] * beta_neg, sigma_neg);
      if(lambda_neg == 0){
        y_new_nominal[n] = -1*exp(y_new[n]);
      }
      else{
        if((lambda_neg * y_new[n] + 1) < 0){
          y_new_nominal[n] = 0; 
          y_exclude[n] = 1;
        }
        else{
          y_new_nominal[n] = -1*((lambda_neg * y_new[n] + 1)^(1/lambda_neg));
        }
      }
    }
  }
}
