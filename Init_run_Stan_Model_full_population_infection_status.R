rm(list = ls())

set.seed(222) # for random initialisation

library(tidyverse)
library(rstan)

# Stan model --------------------------------------------------------------

stan_model <- stan_model("Model_full_population_infection_status.stan")

# Import data -------------------------------------------------------------

df_run <- read_csv("obs_S1_1.csv") 

df_run <- df_run |>
  rename(time = Time) |>
  rename(obs_diff = Diff_Ct) |>
  arrange(ID, time) |>
  as.data.frame()
  
# Df Unique ------------------------------------------------------------------

# Each individual has only 1 row (to extract information later)
df_unique <- df_run |>
  group_by(ID) |>
  slice(1) |>
  ungroup()

n = nrow(df_unique)
max_Ntest = max(df_unique$Ntest)
df_unique$ID_clean <- 1:n # to avoid non-consecutive IDs (Stan does not like it)

# Prepare data for fit ----------------------------------------------------

# Matrix of sampling time of observations
matrix_time <- df_run |>
  select(c(ID, time)) |>
  group_by(ID) |>
  mutate(zz = row_number()) |>
  pivot_wider(values_from = "time",
              names_from = "zz",
              names_prefix = "test_") |>
  ungroup() |>
  select(-ID) |>
  as.matrix()
matrix_time[is.na(matrix_time)] <- 99

# Matrix of observations (in diff from LOD)
matrix_obs_diff <- df_run |>
  select(c(ID, obs_diff)) |>
  group_by(ID) |>
  mutate(zz = row_number()) |>
  pivot_wider(values_from = "obs_diff",
              names_from = "zz",
              names_prefix = "test_") |>
  ungroup() |>
  select(-ID) |>
  as.matrix()
matrix_obs_diff[is.na(matrix_obs_diff)] <- 99

# Matrix of ward : is the observation is censored or not
matrix_censor <- df_run |>
  select(c(ID, censor)) |>
  group_by(ID) |>
  mutate(zz = row_number()) |>
  pivot_wider(values_from = "censor",
              names_from = "zz",
              names_prefix = "test_") |>
  ungroup() |>
  select(-ID) |>
  as.matrix()
matrix_censor[is.na(matrix_censor)] <- 99


ID_clean <- df_unique |>
  select(c(ID, ID_clean)) |>
  ungroup() |>
  as.data.frame()

df_run <- merge(ID_clean, df_run, by = "ID")

# Initialisation ----------------------------------------------------------

# Sample random initialisation value in the priors
init_fun <- function() {
  list(
    mu_Vp = rnorm(1, 25, 2),
    mu_Tp = rnorm(1, 6, 0.25),
    mu_Tc = rnorm(1, 15, 2),
    mu_Tincub = rnorm(1, 5, 0.25),
    Pinf = rbeta(1, 2, 2),
    sigma = truncnorm::rtruncnorm(1, a = 0, mean = 0, sd = 2),
    eta_tilde = matrix(rep(0, 4 * n), ncol = 4),
    eta_omega = truncnorm::rtruncnorm(4, a = 0, mean = 0.15, sd = 0.2)
  )
}

inits_list <- list(init_fun(), init_fun(), init_fun(), init_fun())

inits_to_print <- lapply(inits_list, function(x) x[1:6])
print(inits_to_print)

data_list = list(
  N = nrow(df_run), # Number of observation
  n_id = length(unique(df_run$ID_clean)), # Number of individuals
  lod = 10, # LOD threshold
  Vinf = 0, # VL at infection
  id = df_run$ID_clean, # list of IDs
  time = matrix_time, # time
  y = matrix_obs_diff, # vector of observations (diff compared to LOD)
  censor = matrix_censor, # is the data censored
  tSS = df_unique$tSS,
  nb_random = 4, # number of random effect
  max_Ntest = max_Ntest,  # number of maximum repeated data in the population
  Ntest = df_unique$Ntest
)

# Pre-processing the probability of not being infected (Eq4 in the article) -----------------------------------------------

PTN = 0.998 # arbitrary true negative probability
PFP = 1-PTN # false positive probability
log_PTN = log(PTN) # log of these proba
log_PFP = log(PFP)

data_list$PTN = PTN
data_list$PFP = PFP
data_list$log_PFP= log_PFP
data_list$log_PTN = log_PTN

data_list$always_negative = rep(0,data_list$n_id) 
data_list$number_negative_samples = rep(0,data_list$n_id)
data_list$number_samples = rep(0,data_list$n_id)
data_list$log_prop_negative = rep(0,data_list$n_id)
data_list$log_prop_positive = rep(0,data_list$n_id)

for(i in 1:data_list$n_id){ # loop in individuals
  data_list$number_negative_samples[i] = length(which(data_list$censor[i,] == 1))
  data_list$number_samples[i] = length(which(data_list$censor[i,] != 99))
  
  if(length(which(data_list$censor[i,]==0))== 0){
    data_list$always_negative[i] = 1
    data_list$log_prop_negative[i] = log(PTN^data_list$number_negative_samples[i]) # used in the likelihood in Rstan 
  }
  
  data_list$log_prop_positive[i] = log_PFP*(data_list$number_samples[i] - data_list$number_negative_samples[i]) + log_PTN*data_list$number_negative_samples[i];
  
}

# Run model ---------------------------------------------------------------

fit_startq <- Sys.time()
fit <- rstan::sampling(object = stan_model,
                       data = data_list, 
                       warmup = 200,
                       iter = 400,
                       control = list(adapt_delta=0.95, max_treedepth=15),
                       init = inits_list,
                       chains = 4,
                       cores = 4,
                       seed = 2020)
fit_endq <- Sys.time()

# print the time to run
print(paste0("Fit time: ",difftime(fit_endq, fit_startq, units="min")," mins"))

save(fit, file = "fit_Scenario1_full_population_1.Rsave")

