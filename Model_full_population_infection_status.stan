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

  real always_negative[n_id]; // Vector of individual status (always negative =1, at least one non censored observation = 0)
  real number_negative_samples[n_id]; 
  real number_samples[n_id]; 
  real log_prop_negative[n_id];  
  real log_prop_positive[n_id];  

  // real PFP  ; // proba of false positive
  // real PTN ; // proba of true negative
  // real log_PFP  ;
  // real log_PTN ;
}

transformed data{
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
  real<lower=0, upper=1> Pinf; // Proportion of infected individuals in the population
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
  real num_arg[n_id, max_Ntest, 2]; // array containing individual LL contribution
  num_arg = rep_array(0, n_id, max_Ntest, 2); // Initializing num_arg to 0 (defaut was NaN)
  Ct =  to_matrix(rep_array(0, n_id, max_Ntest)); // // Initializing Ct to 0 (defaut was NaN)
  matrix[n_id,2] ll_Pinf; // matrix of two columns : Pinf and 1-Pinf for each individual
  ll_Pinf[,1] = to_vector(rep_array(log(Pinf), n_id)); // filling the matrix
  ll_Pinf[,2] = to_vector(rep_array(log1m(Pinf), n_id)); // filling the matrix
  for(j in 1:nb_random){
    eta[,j] = eta_tilde[,j]*eta_omega[j];
  }
  Vp = to_array_1d(mu_Vp * exp(eta[ ,1]));
  Tp = to_array_1d(mu_Tp * exp(eta[ ,2]));
  Tc = to_array_1d(mu_Tc * exp(eta[ ,3]));
  Tincub = to_array_1d(mu_Tincub * exp(eta[ ,4]));
  tinf = to_row_vector(tSS) - to_row_vector(Tincub);
  tp = to_row_vector(tinf) + to_row_vector(Tp);


  real log_lik[n_id]; //the log likelighood corresponding to each individual
  real P1 = 0; // some temporary variables
  real P2 = 0; 
    for(i in 1:n_id){  // Looping on each individual
    // Individuals with only negative tests
    if(always_negative[i]==1){
      P1 = 0;
      P2 = 0;
      for (n in 1:Ntest[i]){ 
        Ct[i,n] = CtVinffun(time[i,n], tinf[i], lod, Vinf, tp[i], Tp[i], Tc[i], Vp[i]);
        // if infected, negatives samples are all false negatives
        P1 += normal_lcdf(10 | Ct[i,n], sigma) ;  
      }
      // if uninfected, negatives samples are all true negatives
      P2 = log_prop_negative[i];
      
        log_lik[i] =  log(Pinf*exp(P1) + (1-Pinf)*exp(P2)); 
    }
    // Individuals with >= 1 test
    if(always_negative[i]==0){
      P1 = 0;
      P2 = 0;
      for (n in 1:Ntest[i]){  // Each observation of the individual
        Ct[i,n] = CtVinffun(time[i,n], tinf[i], lod, Vinf, tp[i], Tp[i], Tc[i], Vp[i]);
            
        // if infected
        if(use_obs_ll[i,n] == 1){ // Data is uncensored
            P1 += normal_lpdf(y[i,n] | Ct[i,n], sigma) ;
        } 
        else{ // Data is censored
            P1 += normal_lcdf(10 | Ct[i,n], sigma) ;  
        }
    }
      // if uninfected, positive samples are false positives
      P2 =  log_prop_positive[i];
      log_lik[i] =log(Pinf*exp(P1) + (1-Pinf)*exp(P2));
    }
  }
}

model{
  // Priors //
  mu_Vp ~ normal(25, 2); // hierarchical mean (mu)
  mu_Tp ~ normal(6, 1); // hierarchical mean (mu)
  mu_Tc ~ normal(15, 2); // hierarchical mean (mu)
  mu_Tincub ~ normal(5, 1); // hierarchical mean (mu)
  Pinf ~ beta(2,2); // probability to be infected => proportion of infected individuals in the population
  sigma ~ normal(0,3); // mesurement error
  to_vector(eta_tilde) ~ std_normal();
  to_vector(eta_omega) ~ normal(0.15,0.2); // variance of random effect 
  
  target += sum(log_lik);
}
