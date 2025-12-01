functions {
    // Ctfun returns dCt, the delta Ct below Vinf:
    real CtVinffun(real t, real tinf, int lod, int Vinf, real tp, real Tp, real Tc, real Vp){
		if (t <= tinf)
		  return(Vinf); // Ct = 0
        else if (t >= tinf && t <= tp)
          return(Vp*(t-tinf)/Tp); // Viral load rises before peak: Tp = tp-tinf
        else // Viral load falls after peak: (t > tp)
           return(Vp+(lod-Vp)*(t-tp)/Tc); // Tc = tc-tp
    }
}

data{
  int<lower=0> N; // Number of concatenated data points
  int<lower=0> n_id; // Number of individuals 
  int lod; // Limit of detection of Ct (Ct=10)
  int Vinf; // Ct delta from lod at infection (Ct=0) 
  int<lower=0> id[N]; // Vector marking which datum belongs to which id
  real tSS[n_id]; // Vector of individual time of symptom onset 
  int nb_random; // Total number of parameters with a random effect
  int Ntest[n_id]; // Number of observations for each individual
  int max_Ntest; // Number of maximum of repeated test in the population
  matrix[n_id,max_Ntest] time; // array of dim n_id,max_Ntest containing time of observations 
  matrix[n_id, max_Ntest] y;  // array of dim n_id,max_Ntest containing Ct of observations
  matrix[n_id,max_Ntest] censor; // array of dim n_id,max_Ntest containing censor indication
}

transformed data{
    real log_0002 = log(0.0002); // log proba of false positive
    real log_9998 = log(0.9998); // log proba of true negative
    matrix[n_id,max_Ntest] use_obs_ll;  // 0 = lod; 1 = obs
    for(i in 1:n_id){
      for (n in 1:Ntest[i]){
        if(censor[i,n] == 1){
            use_obs_ll[i,n] = 0;
        } else{
            use_obs_ll[i,n] = 1;
        }
      }
    }
}

parameters{
  real<upper=50> mu_Vp; // fixed effect of Vp 
  real<lower=0> mu_Tp; // fixed effect of Tp
  real<lower=0> mu_Tc; // fixed effect of Tc
  real<lower=0, upper=14> mu_Tincub; // fixed effect of Tincub
  real<lower=0> sigma; // additive measurement error
  matrix[n_id,nb_random] eta_tilde; // Gaussian to center random effect to 0
  real<lower=0> eta_omega[nb_random]; // Variance of random effect for Vp, Tp, Tc, and Tincub
}

transformed parameters{
  real<lower=0, upper=14> Tincub[n_id]; // Incubation period cannot exceed 14 days
  matrix[n_id,nb_random] eta; // matrix of random effects for Vp, Tp, Tc, and Tincub
  real<upper=50> Vp[n_id]; // Vector of individual Vp
  real<lower=0> Tp[n_id]; // Vector of individual Tp
  real<lower=0> Tc[n_id]; // Vector of individual Tc
  row_vector[n_id] tp; // Vector of individual tp
  row_vector[n_id] tinf; // Vector of individual tinf
  matrix<upper=50>[n_id,max_Ntest] Ct; // array of dim n_id,N_max_Tests containing Ct pred of each observation
  real num_arg[n_id, max_Ntest]; // array containing individual LL contribution
  num_arg = rep_array(0, n_id, max_Ntest); // Initializing num_arg to 0 (defaut was NaN)
  Ct =  to_matrix(rep_array(0, n_id, max_Ntest)); // // Initializing Ct to 0 (defaut was NaN)
  for(j in 1:nb_random){
    eta[,j] = eta_tilde[,j]*eta_omega[j];
  }
  Vp = to_array_1d(mu_Vp * exp(eta[ ,1]));
  Tp = to_array_1d(mu_Tp * exp(eta[ ,2]));
  Tc = to_array_1d(mu_Tc * exp(eta[ ,3]));
  Tincub = to_array_1d(mu_Tincub * exp(eta[ ,4]));
  tinf = to_row_vector(tSS) - to_row_vector(Tincub);
  tp = to_row_vector(tinf) + to_row_vector(Tp);
    for(i in 1:n_id){  // Looping on each individual
      for (n in 1:Ntest[i]){  // Each observation of the individual
        Ct[i,n] = CtVinffun(time[i,n], tinf[i], lod, Vinf, tp[i], Tp[i], Tc[i], Vp[i]); // Prediction of the model
        if(time[i,n] < tinf[i]){ // Observations before the start of infection
            if(use_obs_ll[i,n] == 1){  // Data are uncensored
                num_arg[i,n] += log_0002; // proba of false positive
            } else{// Data are censored
                num_arg[i,n] += log_9998; // proba of true negative
            }
        } else{ // Observations after the start of infection
            if(use_obs_ll[i,n] == 1){ // Data are uncensored
                num_arg[i,n] += normal_lpdf(y[i,n] | Ct[i,n], sigma); // If infected 
            } else{ // Data are censored
                num_arg[i,n] += normal_lcdf(10 | Ct[i,n], sigma); // If infected 
            }
        }
    }
  }
}

model{
  // Priors //
  mu_Vp ~ normal(25, 2); // hierarchical mean (mu)
  mu_Tp ~ normal(6, 1); // hierarchical mean (mu)
  mu_Tc ~ normal(15, 2); // hierarchical mean (mu)
  mu_Tincub ~ normal(5, 1); // hierarchical mean (mu)
  sigma ~ normal(0, 3); // mesurement error
  to_vector(eta_tilde) ~ std_normal();
  to_vector(eta_omega) ~ normal(0.15, 0.2); // variance of random effect
  for(i in 1:n_id){ // Likelihood (looped on each observation of each individual)
      target += log_sum_exp(num_arg[i]);
  }
}
