.build_strata <- function(df_hidden, constraints) {
  # Always stratify by (Z, event)
  interaction(paste0("Z", df_hidden$Z), paste0("E", df_hidden$event), drop = TRUE)
}

# Compute the per-(Z, group) maximum allowed time from rows flagged as max==TRUE.
# If a (Z, group) has no max row, treat its allowed max time as +Inf (no cap).
.compute_tmax_table <- function(df_hidden, groups = 0:1) {
  Zvals <- sort(unique(df_hidden$Z))
  mk <- function(z, g) paste0("z", z, "_g", g)
  tmax_tbl <- as.list(rep(Inf, length(Zvals) * length(groups)))
  names(tmax_tbl) <- as.vector(outer(Zvals, groups, Vectorize(function(z,g) mk(z,g))))

  has_max <- !is.null(df_hidden$max) && any(df_hidden$max %in% TRUE)
  if (!has_max) return(tmax_tbl)
  # Use locked_label on max rows to identify subgroup
  idx_max <- which((df_hidden$max %in% TRUE) & !is.na(df_hidden$locked_label))
  if (!length(idx_max)) return(tmax_tbl)
  for (z in Zvals) {
    for (g in groups) {
      ii <- idx_max[df_hidden$Z[idx_max] == z & as.integer(df_hidden$locked_label[idx_max]) == g]
      if (length(ii)) {
        tmax_tbl[[mk(z,g)]] <- suppressWarnings(max(df_hidden$time[ii], na.rm = TRUE))
      }
    }
  }
  tmax_tbl
}

.make_initial_assignment_generic <- function(df_hidden, constraints) {
  # Concise initializer: assume per-(Z, event, group) counts provided; sample directly.
  N <- nrow(df_hidden)
  groups <- 0:1
  Zvals <- sort(unique(df_hidden$Z))
  label <- rep(NA_integer_, N)
  # Respect locked rows (locked==TRUE) and max==TRUE rows: pre-fill labels and adjust pools/counts
  locked_col <- if (!is.null(df_hidden$locked)) df_hidden$locked else rep(FALSE, N)
  max_col    <- if (!is.null(df_hidden$max)) df_hidden$max else rep(FALSE, N)
  locked_lab <- if (!is.null(df_hidden$locked_label)) df_hidden$locked_label else rep(NA, N)
  is_locked_any <- (locked_col %in% TRUE) | (max_col %in% TRUE)
  locked_idx <- which(is_locked_any & !is.na(locked_lab))
  if (length(locked_idx)) {
    label[locked_idx] <- as.integer(locked_lab[locked_idx])
  }
  # Precompute max-time caps per (Z,group)
  tmax_tbl <- .compute_tmax_table(df_hidden, groups)
  key_zg <- function(z,g) paste0("z", z, "_g", g)

  for (z in Zvals) {
    idxZ  <- which(df_hidden$Z == z)
    idxE1 <- idxZ[df_hidden$event[idxZ] == 1L]
    idxE0 <- setdiff(idxZ, idxE1)
    # Remove locked rows from pools; their labels are already set
    locked_here <- which(is_locked_any[idxZ] & !is.na(locked_lab[idxZ]))
    if (length(locked_here)) {
      locked_idxZ <- idxZ[locked_here]
      idxE1 <- setdiff(idxE1, locked_idxZ)
      idxE0 <- setdiff(idxE0, locked_idxZ)
    }
    e1_counts <- setNames(integer(length(groups)), as.character(groups))
    e0_counts <- setNames(integer(length(groups)), as.character(groups))
    for (g in groups) {
      totZ <- as.integer(constraints[[sprintf("n_%d_z%d", g, z)]]);
      evZ  <- as.integer(constraints[[sprintf("n_%d_z%d_ev", g, z)]]);
      # Adjust counts by subtracting locked rows already assigned to group g in this Z and event stratum
      if (is.na(totZ)) totZ <- 0L
      if (is.na(evZ)) evZ <- 0L
      # locked (incl. max) events for this (z,g)
      idxZ_all <- which(df_hidden$Z == z)
      locked_z <- which(is_locked_any[idxZ_all] & !is.na(locked_lab[idxZ_all]))
      if (length(locked_z)) {
        idx_locked_z <- idxZ_all[locked_z]
        locked_lab_z <- as.integer(locked_lab[idx_locked_z])
        locked_ev1 <- idx_locked_z[df_hidden$event[idx_locked_z] == 1L & locked_lab_z == g]
        locked_ev0 <- idx_locked_z[df_hidden$event[idx_locked_z] == 0L & locked_lab_z == g]
        evZ <- max(0L, evZ - length(locked_ev1))
        totZ <- max(0L, totZ - (length(locked_ev1) + length(locked_ev0)))
      }
      e1_counts[as.character(g)] <- max(0L, evZ)
      e0_counts[as.character(g)] <- max(0L, totZ - evZ)
    }
    poolE1 <- idxE1; poolE0 <- idxE0
    # Capacity check per group given max-time caps (exclude locked already)
    for (g in groups) {
      tcap <- tmax_tbl[[key_zg(z,g)]]
      eligE1_g <- poolE1[df_hidden$time[poolE1] <= tcap]
      eligE0_g <- poolE0[df_hidden$time[poolE0] <= tcap]
      if (length(eligE1_g) < e1_counts[as.character(g)]) {
        stop(sprintf("Initialization impossible: not enough eligible E=1 rows for z=%d,g=%d under max-time caps", z, g))
      }
      if (length(eligE0_g) < e0_counts[as.character(g)]) {
        stop(sprintf("Initialization impossible: not enough eligible E=0 rows for z=%d,g=%d under max-time caps", z, g))
      }
    }
    for (gname in sample(names(e1_counts))) {
      ng1 <- as.integer(e1_counts[[gname]])
      if (ng1 > 0) {
        g <- as.integer(gname)
        tcap <- tmax_tbl[[key_zg(z,g)]]
        elig <- poolE1[df_hidden$time[poolE1] <= tcap]
        sel <- if (length(elig)) sample(elig, ng1, replace = FALSE) else integer(0)
        if (length(sel)) { label[sel] <- as.integer(gname); poolE1 <- setdiff(poolE1, sel) }
      }
    }
    for (gname in sample(names(e0_counts))) {
      ng0 <- as.integer(e0_counts[[gname]])
      if (ng0 > 0) {
        g <- as.integer(gname)
        tcap <- tmax_tbl[[key_zg(z,g)]]
        elig <- poolE0[df_hidden$time[poolE0] <= tcap]
        sel <- if (length(elig)) sample(elig, ng0, replace = FALSE) else integer(0)
        if (length(sel)) { label[sel] <- as.integer(gname); poolE0 <- setdiff(poolE0, sel) }
      }
    }
    # Any leftover pools would indicate capacity/accounting issues; fail safe
    if (length(poolE1) || length(poolE0)) {
      stop(sprintf("Initialization failed to assign all rows for z=%d; check constraints vs. max-time caps", z))
    }
  }
  # Check constraints and locked-label consistency for the initial assignment
  if(.labels_satisfy_constraints(label, df_hidden, constraints)){
    return(label)
  } else{
    stop("Failed to generate initial assignment satisfying constraints")
  }
  
  label
}

.propose_swap_preserving_constraints <- function(label_vec, df_hidden, constraints) {
  strata <- .build_strata(df_hidden, constraints)
  idx_by_stratum <- split(seq_len(length(label_vec)), strata)
  # Create combined lock mask: locked==TRUE OR max==TRUE with non-NA locked_label
  N <- nrow(df_hidden)
  locked_col <- if (!is.null(df_hidden$locked)) df_hidden$locked else rep(FALSE, N)
  max_col    <- if (!is.null(df_hidden$max)) df_hidden$max else rep(FALSE, N)
  locked_lab <- if (!is.null(df_hidden$locked_label)) df_hidden$locked_label else rep(NA, N)
  is_locked_any <- ((locked_col %in% TRUE) | (max_col %in% TRUE)) & !is.na(locked_lab)
  # Precompute tmax caps
  groups <- sort(unique(as.integer(label_vec)))
  tmax_tbl <- .compute_tmax_table(df_hidden, groups)
  key_zg <- function(z,g) paste0("z", z, "_g", g)
  # choose a stratum with at least two different labels present (excluding locked)
  cand <- Filter(function(v){
    vv <- v
    vv <- vv[!is_locked_any[vv]]
    length(vv) >= 2L && length(unique(label_vec[vv])) >= 2L
  }, idx_by_stratum)
  if (!length(cand)) return(NULL)
  # Attempt to find a permissible swap that respects tmax caps
  attempt <- 0L
  while (attempt < 100L) {
    ids <- cand[[sample.int(length(cand), 1)]]
    ids <- ids[!is_locked_any[ids]]
    if (length(ids) < 2L) { attempt <- attempt + 1L; next }
    labs <- label_vec[ids]
    ulabs <- unique(labs)
    if (length(ulabs) < 2L) { attempt <- attempt + 1L; next }
    chosen_labs <- sample(ulabs, 2)
    i <- sample(ids[labs == chosen_labs[1]], 1)
    j <- sample(ids[labs == chosen_labs[2]], 1)
    # Check max-time feasibility for swapped labels
    zi <- df_hidden$Z[i]; zj <- df_hidden$Z[j]
    li <- label_vec[i]; lj <- label_vec[j]
    # After swap: i -> lj, j -> li
    t_i <- df_hidden$time[i]; t_j <- df_hidden$time[j]
    cap_i <- tmax_tbl[[key_zg(zi, lj)]]
    cap_j <- tmax_tbl[[key_zg(zj, li)]]
    if ((t_i <= cap_i) && (t_j <= cap_j)) return(c(i, j))
    attempt <- attempt + 1L
  }
  NULL
}

.eval_loss_generic <- function(label_vec, df_hidden, target_stats, loss_fun, target_fn, eps = 1e-8, na_penalty = 1) {
  # Compute estimated summary stats using the provided target_fn, then compute loss
  est_stats <- target_fn(df_hidden, as.integer(label_vec))
  loss_fun(est_stats, target_stats, eps = eps, na_penalty = na_penalty)
}

# Validate that a proposed labeling meets the provided constraints and honors locked rows
.labels_satisfy_constraints <- function(label_vec, df_hidden, constraints) {
  # Check locked rows first
  N <- nrow(df_hidden)
  locked_col <- if (!is.null(df_hidden$locked)) df_hidden$locked else rep(FALSE, N)
  max_col    <- if (!is.null(df_hidden$max)) df_hidden$max else rep(FALSE, N)
  locked_lab <- if (!is.null(df_hidden$locked_label)) df_hidden$locked_label else rep(NA, N)
  is_locked_any <- ((locked_col %in% TRUE) | (max_col %in% TRUE)) & !is.na(locked_lab)
  li <- which(is_locked_any)
  if (length(li)) {
    if (any(as.integer(label_vec[li]) != as.integer(locked_lab[li]))) return(FALSE)
  }
  # Enforce max-time caps: for each row, assigned label must have time <= tmax(Z, label)
  groups <- sort(unique(as.integer(label_vec)))
  tmax_tbl <- .compute_tmax_table(df_hidden, groups)
  key_zg <- function(z,g) paste0("z", z, "_g", g)
  # Rows with labels not in groups get skipped (NA/invalid)
  idx_all <- which(!is.na(label_vec))
  if (length(idx_all)) {
    z_all <- df_hidden$Z[idx_all]
    g_all <- as.integer(label_vec[idx_all])
    t_all <- df_hidden$time[idx_all]
    caps <- mapply(function(z,g) tmax_tbl[[key_zg(z,g)]], z_all, g_all)
    # Treat missing caps as +Inf by design
    if (any(t_all > caps)) return(FALSE)
  }
  # Check only the provided constraint keys
  if (is.null(constraints) || !length(constraints)) return(TRUE)
  Z <- df_hidden$Z
  EV <- df_hidden$event
  for (nm in names(constraints)) {
    target <- constraints[[nm]]
    if (is.null(target)) next
    m <- regexec("^n_([01])_z([01])(?:_ev)?$", nm)
    rr <- regmatches(nm, m)[[1]]
    if (length(rr) == 0) next
    g <- as.integer(rr[2])
    z <- as.integer(rr[3])
    has_ev <- grepl("_ev$", nm)
    idx <- which(as.integer(label_vec) == g & Z == z)
    if (has_ev) idx <- idx[EV[idx] == 1L]
    val <- length(idx)
    if (!identical(val, as.integer(target))) return(FALSE)
  }
  TRUE
}

 

simulated_annealing_swap <- function(
  df_hidden,
  target_stats,
  constraints,
  loss_fun,
  target_fn,
  T0 = 1.0,
  cooling = 0.995,
  min_iter = 1000L,
  max_iter = 10000L,
  random_seed = NULL,
  verbose = TRUE,
  stagnation_stop_iter = 2000L
) {
  if (!is.null(random_seed)) set.seed(random_seed)

  # Feasible initialization
  cur_label <- .make_initial_assignment_generic(df_hidden, constraints)
  # initial labels before diversification are not exposed; keep only post-diversification

  # Diversification phase removed per request
  init_labels <- as.integer(cur_label)

  cur_loss  <- .eval_loss_generic(cur_label, df_hidden, target_stats, loss_fun, target_fn)
  best_lab  <- cur_label
  best_loss <- cur_loss

  T <- T0
  accept_count <- 0L
  last_improve_it <- 0L

  for (it in seq_len(max_iter)) {
    # Only attempt pair swaps within (Z,event) strata
  pair <- .propose_swap_preserving_constraints(cur_label, df_hidden, constraints)
    if (is.null(pair)) {
      T <- T * cooling
      if (verbose && it %% 1000L == 0L) {
        message(sprintf("[it=%d] (no move) T=%.4g cur=%.4g best=%.4g", it, T, cur_loss, best_loss))
      }
      next
    }

    prop_lab <- cur_label
    i <- pair[1]; j <- pair[2]
    tmp <- prop_lab[i]; prop_lab[i] <- prop_lab[j]; prop_lab[j] <- tmp
    # Re-enforce locked rows on proposal
    N <- nrow(df_hidden)
    locked_col <- if (!is.null(df_hidden$locked)) df_hidden$locked else rep(FALSE, N)
    max_col    <- if (!is.null(df_hidden$max)) df_hidden$max else rep(FALSE, N)
    locked_lab <- if (!is.null(df_hidden$locked_label)) df_hidden$locked_label else rep(NA, N)
    is_locked_any <- ((locked_col %in% TRUE) | (max_col %in% TRUE)) & !is.na(locked_lab)
    li <- which(is_locked_any)
    if (length(li)) prop_lab[li] <- as.integer(locked_lab[li])

  prop_loss <- .eval_loss_generic(prop_lab, df_hidden, target_stats, loss_fun, target_fn)
    d <- prop_loss - cur_loss
    accept <- (d <= 0) || (runif(1) < exp(-d / max(T, 1e-12)))
    if (accept) {
      cur_label <- prop_lab
      cur_loss  <- prop_loss
      accept_count <- accept_count + 1L
      if (cur_loss < best_loss) {
        best_loss <- cur_loss
        best_lab  <- cur_label
        last_improve_it <- it
      }
    }

    # Early exit if perfect fit
    if (cur_loss <= 0 || best_loss <= 0) {
      if (verbose) message(sprintf("[it=%d] loss reached 0 -> stopping early", it))
      break
    }

    T <- T * cooling
    if (it %% 100L == 0L) {
      if (verbose && it %% 1000L == 0L) {
        message(sprintf("[it=%d] T=%.4g cur=%.4g best=%.4g acc=%d", it, T, cur_loss, best_loss, accept_count))
        accept_count <- 0L
      }
    }
    # Stop on stagnation: start counting ONLY after min_iter
    if (it >= min_iter) {
      baseline_it <- max(last_improve_it, min_iter)
      if ((it - baseline_it) >= stagnation_stop_iter) {
        if (verbose) message(sprintf("[it=%d] no improvement for %d iterations after min_iter -> stopping", it, stagnation_stop_iter))
        break
      }
    }
  }

  list(
    best_loss = best_loss,
    best_labels  = as.integer(best_lab),
    init_labels = init_labels
  )
}


simulated_annealing_swap_multi <- function(
  df_hidden,
  target_stats,
  constraints,
  seeds,
  loss_fun,
  target_fn,
  T0,
  cooling,
  min_iter,
  max_iter,
  verbose,
  stagnation_stop_iter
) {

  n <- length(seeds)
  results <- vector("list", n)
  losses  <- rep(NA_real_, n)

  if (verbose) {
    message(sprintf("Running SA across %d seeds...", n))
  }

  for (k in seq_len(n)) {
    s <- seeds[k]
    if (verbose) message(sprintf("[seed=%d] starting", s))
    res <- simulated_annealing_swap(
        df_hidden       = df_hidden,
        target_stats    = target_stats,
        constraints     = constraints,
        loss_fun        = loss_fun,
        target_fn       = target_fn,
        T0              = T0,
        cooling         = cooling,
        min_iter        = min_iter,
        max_iter        = max_iter,
        random_seed     = s,
        verbose         = verbose,
        stagnation_stop_iter = stagnation_stop_iter
      )
    results[[k]] <- res
    if (!is.null(res)) losses[k] <- res$best_loss
    if (verbose) {
      if (is.na(losses[k])) {
        message(sprintf("[seed=%d] best_loss: NA (failed)", s))
      } else {
        message(sprintf("[seed=%d] best_loss: %.6g", s, losses[k]))
      }
    }
  }

  good <- which(!is.na(losses))
  if (length(good) == 0L) stop("All SA runs failed for the provided seeds.")

  best_idx  <- good[which.min(losses[good])]
  best_seed <- seeds[best_idx]
  best_res  <- results[[best_idx]]

  best_res$best_seed <- best_seed
  best_res$seed_summary <- data.frame(
    seed = seeds[good], best_loss = losses[good], row.names = NULL
  )
  # Expose only the initial labeling for the best seed at the multi level
  best_res$best_init_label <- results[[best_idx]]$init_labels
  # Remove other initial-label fields for conciseness
  best_res$init_labels <- NULL
  best_res$init_labels_raw <- NULL
  best_res$init_data <- NULL

  if (verbose) message(sprintf("Best across seeds: seed=%d best_loss=%.6g", best_seed, best_res$best_loss))
  best_res
}
