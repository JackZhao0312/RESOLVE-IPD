safe_cox_hr_ci_generic <- function(df_subset) {
  if (nrow(df_subset) == 0L || length(unique(df_subset$Z)) < 2L) {
    return(c(HR = NA_real_, LCL = NA_real_, UCL = NA_real_))
  }
  fit <- tryCatch(
    survival::coxph(survival::Surv(time, event) ~ Z, data = df_subset),
    error = function(e) NULL
  )
  if (is.null(fit)) return(c(HR = NA_real_, LCL = NA_real_, UCL = NA_real_))
  beta <- stats::coef(fit)
  if (length(beta) == 0L || anyNA(beta)) return(c(HR = NA_real_, LCL = NA_real_, UCL = NA_real_))
  hr <- as.numeric(exp(beta[1]))
  ci <- tryCatch(stats::confint(fit), error = function(e) NULL)
  if (is.null(ci) || anyNA(ci[1,])) return(c(HR = hr, LCL = NA_real_, UCL = NA_real_))
  c(HR = hr, LCL = exp(ci[1,1]), UCL = exp(ci[1,2]))
}

safe_km_median_ci_generic <- function(df_subset) {
  if (nrow(df_subset) == 0L) return(c(median = NA_real_, LCL = NA_real_, UCL = NA_real_))
  fit <- survival::survfit(survival::Surv(time, event) ~ 1, conf.type = "log-log", data = df_subset)
  tab <- tryCatch(summary(fit)$table, error = function(e) NULL)
  med <- NA_real_; lcl <- NA_real_; ucl <- NA_real_
  if (!is.null(tab)) {
    if (is.null(dim(tab))) {
      med <- suppressWarnings(as.numeric(tab["median"]))
      lcl <- suppressWarnings(as.numeric(tab["0.95LCL"]))
      ucl <- suppressWarnings(as.numeric(tab["0.95UCL"]))
    } else {
      cn <- colnames(tab)
      med_col <- which(cn == "median")
      lcl_col <- which(cn %in% c("0.95LCL", "LCL"))
      ucl_col <- which(cn %in% c("0.95UCL", "UCL"))
      if (length(med_col)) med <- suppressWarnings(as.numeric(tab[1, med_col]))
      if (length(lcl_col)) lcl <- suppressWarnings(as.numeric(tab[1, lcl_col[1]]))
      if (length(ucl_col)) ucl <- suppressWarnings(as.numeric(tab[1, ucl_col[1]]))
    }
  }
  c(median = med, LCL = lcl, UCL = ucl)
}

compute_summary_stat_simdata <- function(data, label) {
  # data has columns: time, event, Z; label is integer vector (0/1)
  stopifnot(nrow(data) == length(label))
  df <- data
  df$label <- as.integer(label)
  labs <- sort(unique(df$label))
  labs <- labs[labs %in% c(0L, 1L)]
  out <- list(
    target_HR_0 = NA_real_, target_HR_0_lcl = NA_real_, target_HR_0_ucl = NA_real_,
    target_HR_1 = NA_real_, target_HR_1_lcl = NA_real_, target_HR_1_ucl = NA_real_,
    target_0_z0_median = NA_real_, target_0_z0_median_lcl = NA_real_, target_0_z0_median_ucl = NA_real_,
    target_0_z1_median = NA_real_, target_0_z1_median_lcl = NA_real_, target_0_z1_median_ucl = NA_real_,
    target_1_z0_median = NA_real_, target_1_z0_median_lcl = NA_real_, target_1_z0_median_ucl = NA_real_,
    target_1_z1_median = NA_real_, target_1_z1_median_lcl = NA_real_, target_1_z1_median_ucl = NA_real_
  )
  for (g in labs) {
    ss <- df[df$label == g, , drop = FALSE]
    hr <- safe_cox_hr_ci_generic(ss)
    out[[paste0("target_HR_", g)]] <- round(hr[["HR"]],2)
    out[[paste0("target_HR_", g, "_lcl")]] <- round(hr[["LCL"]],2)
    out[[paste0("target_HR_", g, "_ucl")]] <- round(hr[["UCL"]],2)
    for (z in c(0L, 1L)) {
      ssz <- ss[ss$Z == z, , drop = FALSE]
      med <- safe_km_median_ci_generic(ssz)
      base <- paste0("target_", g, "_z", z, "_median")
      out[[base]] <- round(med[["median"]],1)
      out[[paste0(base, "_lcl")]] <- round(med[["LCL"]],1)
      out[[paste0(base, "_ucl")]] <- round(med[["UCL"]],1)
    }
  }
  out
}