compute_rel_worst_loss <- function(
	est_stats,
	target_stats,
	eps = 1e-8,
	na_penalty = 1
){
	# Only compare keys present and non-NULL in target_stats
	cand <- intersect(names(est_stats), names(target_stats))
	cand <- cand[!vapply(cand, function(k) is.null(target_stats[[k]]), logical(1))]

	rel_losses <- vapply(cand, function(nm) {
		obs_val <- est_stats[[nm]]
		target_val <- target_stats[[nm]]
		if (!is.finite(obs_val) || !is.finite(target_val)) return(as.numeric(na_penalty))
		denom <- max(abs(target_val), eps)
		loss_val <- abs(obs_val - target_val) / denom
		if (!is.finite(loss_val)) as.numeric(na_penalty) else loss_val
	}, numeric(1L))

	if (anyNA(rel_losses)) rel_losses[is.na(rel_losses)] <- as.numeric(na_penalty)
	max(rel_losses)
}

