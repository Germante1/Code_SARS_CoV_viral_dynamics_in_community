rm(list = ls())

# Packages ---------------------------------------------------------------
library(tidyverse)

# Functions ---------------------------------------------------------------

# Viral dynamics function (Ct-like trajectory with Vinf formulation)
ct_vinf_fun <- function(t, tinf, lod, Vinf, tp, Tp, Tc, Vp) {
  if (t >= tinf && t <= tp) {
    return(Vp * (t - tinf) / Tp)
  } else if (t > tp) {
    return(Vp + (Vinf - lod - Vp) * (t - tp) / Tc)
  } else {
    return(0)
  }
}

# Simulate individual viral/temporal parameters for infected subjects
simulate_individual_parameters <- function(n, mu, eta_sd) {
  eta <- replicate(4, rnorm(n, 0, eta_sd)) |> as.data.frame()
  names(eta) <- c("eta_Vp", "eta_Tp", "eta_Tc", "eta_Tincub")
  
  param <- tibble(
    Vp = mu[1] * exp(eta$eta_Vp),
    Tp = mu[2] * exp(eta$eta_Tp),
    Tc = mu[3] * exp(eta$eta_Tc),
    Tincub = mu[4] * exp(eta$eta_Tincub)
  )
  
  tSS <- rep(0, n) # because time is relative to symptom onset
  
  tinf <- tSS - param$Tincub
  tp <- tinf + param$Tp
  tc <- tp + param$Tc
  
  tibble(ID = 1:n, tinf = tinf, tp = tp,
         tc = tc, tSS = tSS) |>
    bind_cols(param, eta)
}

# Simulate a full time grid for infected individuals
simulate_full_grid <- function(param_df, seq_time, lod, Vinf, sigma_obs) {
  n <- nrow(param_df)
  df <- expand_grid(ID = 1:n, t = seq_time) |>
    left_join(param_df, by = "ID")
  
  traj <- mapply(
    function(i, t) ct_vinf_fun(
      t = t, tinf = df$tinf[i], lod = lod,
      Vinf = Vinf, tp = df$tp[i], Tp = df$Tp[i],
      Tc = df$Tc[i], Vp = df$Vp[i]),
    seq_len(nrow(df)),
    df$t)
  
  obs <- traj + rnorm(nrow(df), 0, sigma_obs)
  
  df |>
    mutate(
      traj = 50 - traj,
      obs_raw = 50 - obs,
      obs = pmin(obs_raw, 40),
      censor = as.integer(obs >= 40)) |>
    select(ID, t, traj, obs, censor, everything())
}

# PARAMETERS USED IN ALL SCENARIOS ----------------------------------------

mu <- c(25, 6, 15, 5) # mean values for Vp, Tp, Tc, Tincub
lod <- 40 # limit of detection Ct
Vinf <- 50 # Ct at infection
seq_time <- -20:35 # window of observation
j_max <- 1 # number of datasets per scenario
sigma_obs <- 2 # observation noise SD

# SCENARIO 1  ----------------------------------------------------------

set.seed(111)

n_inf <- 60   # number of infected individuals
n_neg <- 60   # mirror negative population

for (j in 1:j_max) {
  
  param_inf <- simulate_individual_parameters(n_inf, mu, 0.15) |>
    mutate(infected = 1)
  
  df_inf <- param_inf |>
    group_by(ID) |>
    reframe(Time = seq(floor(tinf) - 1, ceiling(tc) + 1, 1)) |>
    left_join(param_inf, by = "ID") |>
    mutate(traj = pmap_dbl(
      list(Time, tinf, tp, Tp, Tc, Vp),
      ~ ct_vinf_fun(..1, ..2, lod, Vinf, ..3, ..4, ..5, ..6)
    ),
      traj = 50 - traj,
      obs = pmin(traj + rnorm(n(), 0, sigma_obs), 40),
      censor = as.integer(obs >= 40),
      Diff_Ct = 50 - obs) |>
    select(ID, Time, traj, obs, Diff_Ct, censor, tSS, infected)
  
  df_neg <- df_inf |>
    mutate(ID = ID + n_inf,
           infected = 0,
           traj = 50,
           obs = 40,
           censor = 1,
           Diff_Ct = 10) |> # 50 - 40
    select(ID, Time, traj, obs, Diff_Ct, censor, tSS, infected)
  
  df_obs <- bind_rows(df_inf, df_neg) |>
    group_by(ID) |>
    mutate(Ntest = n(),
           ID_test = row_number(),
           Tag_tinf = ifelse(sum(censor) != Ntest, 1, 0)) |>
    ungroup() |>
    arrange(ID, Time)|>
    select(ID, Ntest, ID_test, Time, traj, obs, Diff_Ct, censor, tSS, infected, Tag_tinf) 
  
  
  df_param <- bind_rows(
    param_inf |> select(ID, tSS, tinf, tp, tc, Vp, Tp, Tc, Tincub, infected),
    tibble(ID = (n_inf + 1):(n_inf + n_neg),
           tSS = 0, tinf = NA, tp = NA, tc = NA,
           Vp = NA, Tp = NA, Tc = NA, Tincub = NA,
           infected = 0)) |>
    arrange(ID)

  write.csv(df_param, paste0("BDD_article/param_S1_", j, ".csv"), row.names = FALSE)
  write.csv(df_obs, paste0("BDD_article/obs_S1_", j, ".csv"), row.names = FALSE)
}

# SCENARIO 2  -------------------------------------------------------------

set.seed(222)

n_inf <- 2000
n_neg <- 2000

probs_obs <- c(`1` = 0.75, `2` = 0.20, `3` = 0.04, `4` = 0.01)

for (j in 1:j_max) {
  
  param_inf <- simulate_individual_parameters(n_inf, mu, 0.15) |>
    mutate(infected = 1) |>
    mutate(tSS = 0)
  
  df_inf <- param_inf |>
    transmute(ID, infected, tinf, tp, Tp, Tc, Vp,
              t = map2(tinf, tc,
                       ~ sort(unique(round(runif(
                         n = sample(1:4, 1, prob = probs_obs),
                         min = .x - 2, max = .y + 2)))))) |>
    unnest(t, names_repair = "minimal") |>
    mutate(traj = pmap_dbl(
        list(t, tinf, tp, Tp, Tc, Vp),
        ~ ct_vinf_fun(..1, ..2, lod, Vinf, ..3, ..4, ..5, ..6)),
      traj = 50 - traj,
      obs = pmin(traj + rnorm(n(), 0, sigma_obs), 40),
      censor = as.integer(obs >= 40),
      Diff_Ct = 50 - obs) |>
    mutate(tSS = 0) |>
    rename(Time = t) |>
    select(ID, Time, traj, obs, Diff_Ct, censor, tSS, infected)
  
  df_neg <- df_inf |>
    mutate(ID = ID + n_inf,
           infected = 0,
           traj = 50,
           obs = 40,
           censor = 1,
           Diff_Ct = 10) |>
    select(ID, Time, traj, obs, Diff_Ct, censor, tSS, infected)
  
  df_obs <- bind_rows(df_inf, df_neg) |>
    group_by(ID) |>
    mutate(Ntest = n(),
           ID_test = row_number(),
           Tag_tinf = ifelse(sum(censor) != Ntest, 1, 0)) |>
    ungroup() |>
    arrange(ID, Time)|>
    select(ID, Ntest, ID_test, Time, traj, obs, Diff_Ct, censor, tSS, infected, Tag_tinf) 
  
  
  df_param <- bind_rows(
    param_inf |> select(ID, tSS, tinf, tp, tc, Vp, Tp, Tc, Tincub, infected),
    tibble(ID = (n_inf + 1):(n_inf + n_neg),
           tSS = 0, tinf = NA, tp = NA, tc = NA,
           Vp = NA, Tp = NA, Tc = NA, Tincub = NA,
           infected = 0)) |>
    arrange(ID)
  
  write.csv(df_param, paste0("BDD_article/param_S2_", j, ".csv"), row.names = FALSE)
  write.csv(df_obs, paste0("BDD_article/obs_S2_", j, ".csv"), row.names = FALSE)
}

# SCENARIO 3  -------------------------------------------------------------

set.seed(333)

library(data.table)
df_real <- as.data.frame(fread("BDD_article/community_PCR_tests_dataset.csv"))
detach("package:data.table", unload = TRUE)

df_real_pos <- df_real |>
  filter(infected == 1)

df_real_neg <- df_real |>
  filter(infected == 0)

id_pos <- unique(df_real_pos$ID)
id_neg <- unique(df_real_neg$ID)

# Directly sampling individuals from the community dataset 
# to randomly extract number of test and time of these tests
# Here we mimic as much as possible our datatset
func_sampling <- function(ID_list, n) {
  lapply(1:50, function(x) sample(ID_list, size = n, replace = FALSE))
}

n_inf = 2000
n_neg = 2000

# Sampling 2,O00 indiviuals from infected and always negative
sample_pos <- func_sampling(id_pos, n_inf)
sample_neg <- func_sampling(id_neg, n_neg)

for (j in 1:j_max) {
  
  temp_pos <- df_real_pos |>
    filter(ID %in% sample_pos[[j]]) |>
    select(ID, ID_test, time_since_tSS, infected) |>
    rename(ID_old = ID)
  
  map_pos <- temp_pos |>
    group_by(ID_old) |>
    slice(1) |>
    ungroup() |>
    mutate(ID = row_number()) |>
    select(ID_old, ID)
  
  temp_pos <- temp_pos |>
    inner_join(map_pos, by = "ID_old") |>
    rename(Time = time_since_tSS)
  
  temp_neg <- df_real_neg |>
    filter(ID %in% sample_neg[[j]]) |>
    select(ID, ID_test, time_since_tSS, infected) |>
    rename(ID_old = ID)
  
  map_neg <- temp_neg |>
    group_by(ID_old) |>
    slice(1) |>
    ungroup() |>
    mutate(ID = nrow(map_pos) + row_number()) |>
    select(ID_old, ID)
  
  temp_neg <- temp_neg |>
    inner_join(map_neg, by = "ID_old") |>
    rename(Time = time_since_tSS)
  
  param_inf <- simulate_individual_parameters(n_inf, mu, 0.15) |>
    mutate(infected = 1)
  
  df_grid <- simulate_full_grid(param_inf, seq_time, lod, Vinf, sigma_obs)
  
  df_inf <- df_grid |>
    semi_join(temp_pos, by = c("ID", "t" = "Time")) |>
    rename(Time = t) |>
    mutate(Diff_Ct = 50 - obs)
  
  df_neg <- temp_neg |>
    mutate(traj = 50, obs = 40, censor = 1,
           Diff_Ct = 10) |>
    select(ID, Time, traj, obs, Diff_Ct, censor, infected)
  
  df_obs <- bind_rows(df_inf, df_neg) |>
    group_by(ID) |> 
    mutate(Ntest = n(),
           ID_test = row_number(),
           Tag_tinf = ifelse(sum(censor) != Ntest, 1, 0)) |>
    ungroup() |>
    arrange(ID, Time) |>
    select(ID, Ntest, ID_test, Time, traj, obs, Diff_Ct, censor, tSS, infected, Tag_tinf) 
  
  df_param_inf <- param_inf |>
    select(ID, tSS, tinf, tp, tc, Vp, Tp, Tc, Tincub, infected)
  
  df_param_neg <- tibble(
    ID = n_inf+1:(n_inf+n_neg),
    tSS = 0,
    tinf = NA, tp = NA, tc = NA,
    Vp = NA, Tp = NA, Tc = NA, Tincub = NA,
    infected = 0)
  
  df_param <- bind_rows(df_param_inf, df_param_neg) |>
    arrange(ID)
  
  write.csv(df_param, paste0("BDD_article/param_S3_", j, ".csv"), row.names = FALSE)
  write.csv(df_obs, paste0("BDD_article/obs_S3_", j, ".csv"), row.names = FALSE)
}
