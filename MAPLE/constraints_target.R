library(dplyr)

# ===============================
# Specify constraints for labels g in {0,1}
# ===============================
constraints <- list(
  n_0_z1_ev = df_full %>% filter(label == 0 & Z == 1 & event == 1) %>% nrow(),
  n_0_z1    = df_full %>% filter(label == 0 & Z == 1) %>% nrow(),
  n_0_z0_ev = df_full %>% filter(label == 0 & Z == 0 & event == 1) %>% nrow(),
  n_0_z0    = df_full %>% filter(label == 0 & Z == 0) %>% nrow(),
  n_1_z1_ev = df_full %>% filter(label == 1 & Z == 1 & event == 1) %>% nrow(),
  n_1_z1    = df_full %>% filter(label == 1 & Z == 1) %>% nrow(),
  n_1_z0_ev = df_full %>% filter(label == 1 & Z == 0 & event == 1) %>% nrow(),
  n_1_z0    = df_full %>% filter(label == 1 & Z == 0) %>% nrow()
)



# ===============================
# Specify target statistics
# ===============================
target_stats <- compute_summary_stat_simdata(
  data = df_full,
  label = df_full$label
)
