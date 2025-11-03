library(survival)
library(dplyr)
library(here)


source(here("MAPLE", "target.R"))
source(here("MAPLE", "loss.R"))
source(here("MAPLE", "SA.R"))
source(here("MAPLE", "simdata.R"))




# ===============================
# 1) Read data
# ===============================
N <- 400
n_g0 <- 200
df_full <- generate_simdata(N = N, n_g0 = n_g0, seed = 0) # low = 0, high = 1
source(here("MAPLE", "constraints_target.R"))



# ===============================
# 2) Format the lock columns
# ===============================
# change df_full column group to locked_label
df_full$locked_label <- NA
df_full$locked <- NA
df_full$max <- NA
# stratify by label and Z to get maximum time
for (lbl in unique(df_full$label)) {
  for (z in unique(df_full$Z)) {
    idx_subset <- which(df_full$label == lbl & df_full$Z == z)
    times_subset <- df_full$time[idx_subset]
    max_time <- max(times_subset)
    idx_max_time <- idx_subset[which(times_subset == max_time)]
    df_full$locked_label[idx_max_time] <- lbl
    df_full$locked[idx_max_time] <- TRUE
    df_full$max[idx_max_time] <- TRUE
  }
}


# ===============================
# ===============================
# 3) Run simulated annealing
# ===============================
res_sa <- simulated_annealing_swap_multi(
  df_hidden      = df_full,
  target_stats   = target_stats,
  constraints    = constraints,
  seeds = 2,
  loss_fun = compute_rel_worst_loss,
  target_fn = compute_summary_stat_simdata,
  T0 = 1.0,
  cooling = 0.998,
  min_iter = 10000,
  max_iter = 20000,
  stagnation_stop_iter = 3000,
  verbose = TRUE)


save_res <- list(
  best_loss = res_sa$best_loss,
  best_labels = res_sa$best_labels,
  best_seed = res_sa$best_seed
)

