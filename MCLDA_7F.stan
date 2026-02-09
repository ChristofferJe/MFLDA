
data {
  int<lower=1> K;                                  // num signatures
  int<lower=1> C;                                  // num channels
  int<lower=1> M;                                  // num patients
  int<lower=1> N;                                  // total num observation pr. channel
  array[C] int<lower=1> T_c;                       // number of types in each channel
  array[N] int<lower=1> obs_c1;                    // observations for channel 1
  array[N] int<lower=1> obs_c2;                    // observations for channel 2
  array[N] int<lower=1> obs_c3;                    // observations for channel 3
  array[N] int<lower=1> obs_c4;                    // observations for channel 4
  array[N] int<lower=1> obs_c5;                    // observations for channel 5
  array[N] int<lower=1> obs_c6;                    // observations for channel 6
  array[N] int<lower=1> obs_c7;                    // observations for channel 7
  array[N] int<lower=1> pat;                       // patient ID for observation n
  vector<lower=0>[K] alpha;                        // signature distribution prior
  vector<lower=0>[T_c[1]] beta_c1;                 // prior of distribution over types for channel 1
  vector<lower=0>[T_c[2]] beta_c2;                 // prior of distribution over types for channel 2
  vector<lower=0>[T_c[3]] beta_c3;                 // prior of distribution over types for channel 3
  vector<lower=0>[T_c[4]] beta_c4;                 // prior of distribution over types for channel 4
  vector<lower=0>[T_c[5]] beta_c5;                 // prior of distribution over types for channel 5
  vector<lower=0>[T_c[6]] beta_c6;                 // prior of distribution over types for channel 6
  vector<lower=0>[T_c[7]] beta_c7;                 // prior of distribution over types for channel 7
}


parameters {
  array[M] simplex[K] theta;                       // topic dist for patient m
  array[K] simplex[T_c[1]] phi_c1;                  // word dist for topic k for channel 1
  array[K] simplex[T_c[2]] phi_c2;                  // word dist for topic k for channel 2
  array[K] simplex[T_c[3]] phi_c3;                  // word dist for topic k for channel 3
  array[K] simplex[T_c[4]] phi_c4;                  // word dist for topic k for channel 4
  array[K] simplex[T_c[5]] phi_c5;                  // word dist for topic k for channel 5
  array[K] simplex[T_c[6]] phi_c6;                  // word dist for topic k for channel 6
  array[K] simplex[T_c[7]] phi_c7;                  // word dist for topic k for channel 7
}

model {
  for (m in 1:M){
    theta[m] ~ dirichlet(alpha);        // prior 
  }

  for (k in 1:K){
    phi_c1[k] ~ dirichlet(beta_c1);
    phi_c2[k] ~ dirichlet(beta_c2);
    phi_c3[k] ~ dirichlet(beta_c3);
    phi_c4[k] ~ dirichlet(beta_c4);
    phi_c5[k] ~ dirichlet(beta_c5);
    phi_c6[k] ~ dirichlet(beta_c6);
    phi_c7[k] ~ dirichlet(beta_c7);
  }

  vector[N] log_likelihood;
  for (n in 1:N){
    array[K] real gamma;
    for (k in 1:K){
      gamma[k] = (log(theta[pat[n]][k]) 
                + log(phi_c1[k][obs_c1[n]])
                + log(phi_c2[k][obs_c2[n]])
                + log(phi_c3[k][obs_c3[n]])
                + log(phi_c4[k][obs_c4[n]])
                + log(phi_c5[k][obs_c5[n]])
                + log(phi_c6[k][obs_c6[n]])
                + log(phi_c7[k][obs_c7[n]]));
    }
    log_likelihood[n] = log_sum_exp(gamma); 
  }
  target += sum(log_likelihood);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N){
    array[K] real gamma;
    for (k in 1:K){
      gamma[k] = (log(theta[pat[n]][k]) 
                + log(phi_c1[k][obs_c1[n]])
                + log(phi_c2[k][obs_c2[n]])
                + log(phi_c3[k][obs_c3[n]])
                + log(phi_c4[k][obs_c4[n]])
                + log(phi_c5[k][obs_c5[n]])
                + log(phi_c6[k][obs_c6[n]])
                + log(phi_c7[k][obs_c7[n]]));
    }
    log_lik[n] = log_sum_exp(gamma); 
  } 
}
