## adapted from: https://osf.io/bcr4u

data {
	int N;
	vector[N] y;
	vector[N] x;
	
	int<lower=0> N_new;
  vector[N_new] x_new;
}

parameters {
	real alpha_pos;       // intercepts
	real alpha_neg;
	real beta_pos;        // slopes
	real beta_neg;
	real lambda_pos;      // box-cox hyper parameter
	real lambda_neg;
	// simplex[2] p;         // mixing proportions
	real<lower=0> sigma_pos;
	real<lower=0> sigma_neg;
}

// transformed parameters {
//   vector y_boxcox;
//   
//   if(y[i] >= 0){
//     if(lambda_pos == 0){
//       y_boxcox[i] = log(y[i])
//     }
//     else{
//       y_boxcox[i] = (y[i]^lambda_pos - 1)/lambda_pos 
//     }
//   }
//   
//   else{
//     if(lambda_pos == 0){
//       y_boxcox[i] = log(-1*y[i])
//     }
//     else{
//       y_boxcox[i] = ((-1*y[i])^lambda_pos - 1)/lambda_pos 
//     }
//   }
//   
// }

model {
  lambda_pos ~ normal(0,2);
  lambda_neg ~ normal(0,2);
  alpha_pos ~ normal(0,2);
  alpha_neg ~ normal(0,2);
  beta_pos ~ normal(0,2);
  beta_neg ~ normal(0,2);
  sigma_pos ~ cauchy(0,5);
  sigma_neg ~ cauchy(0,5);
  
  for(i in 1:N){
    if(y[i] >= 0){ // positive y values
      if(lambda_pos == 0){
        log(y[i]) ~ normal(alpha_pos + x[i]*beta_pos, sigma_pos);
        target += -log(y[i]);// + log(p);
      }
      else {
        (y[i]^lambda_pos - 1)/lambda_pos ~ normal(alpha_pos + x[i]*beta_pos, sigma_pos); // transform y using box-cox
        target += (lambda_pos - 1)*log(y[i]);// + log(p);
      }
    }
    else{  // negative y values
      if(lambda_neg == 0){
        log(-1*y[i]) ~ normal(alpha_neg + x[i]*beta_neg, sigma_neg);
        target += -log(-1*y[i]);// + log(1-p);
      } 
      else {
        ((-1*y[i])^lambda_neg - 1)/lambda_neg ~ normal(alpha_neg + x[i]*beta_neg, sigma_neg); // transform y using box-cox
        target += (lambda_neg - 1)*log(-1*y[i]);// + log(1-p);
      }
    }
  }
}

// for (n in 1:N) {
//   target += log_sum_exp(log(p) + normal_lpdf(y[n] | alpha_pos + x[i]*beta_pos, sigma_pos),
//                         log(1-p) + normal_lpdf(y[n] | alpha_neg + x[i]*beta_neg, sigma_neg));
// }


generated quantities { // for prediction intervals
  vector[N_new] y_new;
  vector[N_new] y_new_nominal;
  real<lower=0,upper=1> pos_flag; 
  
  for (n in 1:N_new){
    pos_flag = bernoulli_rng(0.86); // p
    if(pos_flag){
      y_new[n] = normal_rng(alpha_pos + x_new[n] * beta_pos, sigma_pos);
      if(lambda_pos == 0){
        y_new_nominal[n] = exp(y_new[n]);
      }
      else{
        if((lambda_pos * y_new[n] + 1) < 0){
          y_new_nominal[n] = 0; // double check this is the assumption we want to make
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
          y_new_nominal[n] = 0; // double check this is the assumption we want to make
        }
        else{
          y_new_nominal[n] = -1*((lambda_neg * y_new[n] + 1)^(1/lambda_neg));
        }
      }
    }
  }
}
