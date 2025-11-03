generate_simdata <- function(N, n_g0, seed) {
  
  set.seed(seed)
  
  # group 0 count and group 1 inferred
  n_g1 <- N - n_g0
  label <- c(rep(0L, n_g0), rep(1L, n_g1))  # 0/1 group label
  Z <- rbinom(N, 1, 0.5)
  
  # True treatment effects by label groups (targets)
  HR_g0_true <- 0.9
  HR_g1_true <- 0.7
  beta_g0   <- log(HR_g0_true)
  beta_g1   <- log(HR_g1_true)
  
  # Baseline hazard and censoring
  lambda0 <- 0.1
  censor_time <- ifelse(rbinom(N, 1, 0.9) == 1, runif(N, min = 12, max = 24), runif(N, min = 0, max = 12))
  
  # Individual event rate respecting label-specific treatment HR
  linpred <- ifelse(label == 1L, beta_g1 * Z, beta_g0 * Z)
  rate    <- lambda0 * exp(linpred)
  T_event <- rexp(N, rate = rate)
  
  time   <- pmin(T_event, censor_time)
  event <- as.integer(T_event <= censor_time)
  
  df_full <- data.frame(time = time, event = event, Z = Z, label = label)
  df_full <- df_full[sample(nrow(df_full)), ]
  rownames(df_full) <- NULL
  df_full
}