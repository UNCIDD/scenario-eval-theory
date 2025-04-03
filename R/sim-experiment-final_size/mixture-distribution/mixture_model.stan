## adapted from: https://osf.io/bcr4u

data {
	int N;
	vector[N] y;
	vector[N] x;
	
	int<lower=0> N_new;
  vector[N_new] x_new;
}

parameters {
	real alpha;
	real beta;
	real lambda;
	real<lower=0> sigma;
}

model {
	lambda ~ normal(0,2);
	alpha ~ normal(0,2);
	beta ~ normal(0,2);
	sigma ~ cauchy(0,5);

	if(lambda == 0){
		for(i in 1:N){
			log(y[i]) ~ normal(x[i]*beta,sigma);
			target += -log(y[i]);
		}
	} else {
		for(i in 1:N){
			(y[i]^lambda - 1)/lambda ~ normal(alpha + x[i]*beta,sigma);
			target += (lambda - 1)*log(y[i]);
		}
	}
}

generated quantities { // for prediction intervals
  vector[N_new] y_new;
  for (n in 1:N_new)
    y_new[n] = normal_rng(alpha + x_new[n] * beta, sigma);
}
