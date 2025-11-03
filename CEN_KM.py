import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from typing import List, Dict, Tuple, Optional, Union, Any, Iterable

def _km_survival_at(times: np.ndarray, events: np.ndarray, t0: float) -> float:
    """ 
    Calculates the Kaplan–Meier survival estimate
    events: 1 = event (death), 0 = right-censored
    times: array of observed times (events or censoring)
    events: array of 1=event (death), 0=censor
    t0: the time at which we want S(t0)
    """

    mask_event = (events == 1) # select only event (death) times
    event_times = np.unique(times[mask_event]) # unique event times
    event_times = event_times[event_times <= t0] # keep the events that occur before or at t0
    if event_times.size == 0:
        return 1.0 # if no event, survival = 1.0
    S = 1.0
    # Loop through each event time (u) up to t0
    for u in np.sort(event_times):
        r_u = np.sum(times >= u) # number at risk at u
        d_u = np.sum((times == u) & (events == 1)) # number of death at u
        if r_u > 0:
            S *= (1.0 - d_u / r_u) # update KM survival
    return float(S)

def _try_add_censors_in_bin(
    base_t: np.ndarray,
    base_e: np.ndarray,
    p0: int,
    n: int,
    S_target: float,
    t_prev: float,
    t_i: float,
    bin_cens_times: np.ndarray,
    tol: float,
    cap: int,
    rng,
    debug: bool = True
) -> Tuple[np.ndarray, np.ndarray, int, float, int]:
    """
    Try to reduce |KM - S_target| by adding extra censors in [t_prev, t_i).
    If the cap is reached and the tolerance still not satisfied, 
    return the original arrays (no censors added).
    Returns (new_t, new_e, new_p, est_S, extras_added).
    """
    # Backup original
    t_orig = base_t.copy()
    e_orig = base_e.copy()
    p_orig = int(p0)
    est_orig = _km_survival_at(t_orig, e_orig, t_i)

    # Working copies
    t_work = t_orig.copy()
    e_work = e_orig.copy()
    p_work = p_orig
    est_curr = est_orig
    extras = 0
    can_sample = bin_cens_times.size > 0

    while (abs(est_curr - S_target) > tol) and (extras < cap) and (p_work < n):
        if not can_sample:
            break

        c_time = rng.choice(bin_cens_times)
        t_work[p_work] = c_time
        e_work[p_work] = 0
        p_work += 1
        extras += 1
        est_curr = _km_survival_at(t_work, e_work, t_i)

        if debug:
            print(f'Added censor at {c_time:.4f}, new est_S={est_curr:.6f}, error = {est_curr - S_target:.6f}')

    # If still outside tolerance after reaching cap, revert to original
    if (abs(est_curr - S_target) > tol) and (extras >= cap):
        if debug:
            print(f"Reached censor cap ({cap}) but |est_S - target|={abs(est_curr - S_target):.6f} > tol={tol:.6f}, reverting all.")
        return t_orig, e_orig, p_orig, est_orig, 0

    return t_work, e_work, p_work, est_curr, extras

# commit pointer already handled in branches/fallbacks; proceed to next bin
def _at_risk_count(t_arr: np.ndarray, tau: float) -> int:
    # at-risk just prior to tau: times >= tau
    return int(np.sum(t_arr >= tau))

def _nextbefore(x: float) -> float:
    # the largest float strictly less than x (tiny epsilon step)
    return float(np.nextafter(x, -np.inf))

def _align_single_boundary(
    IPD_t: np.ndarray,
    IPD_e: np.ndarray,
    tau: float,
    target_at_risk: int,
    prev_tau: float,
    p: int,
    eps: float = 1e-9,
    rng: Optional[np.random.Generator] = None,
    interval_length: float = 0.1,
    debug: bool = False
):
    """
    Adjust censor times (only) to match target at-risk at time `tau`.

    New behavior:
    - Takes `p`: the next available free slot (index) to write NEW censor rows.
    - If delta > 0 (need to INCREASE at-risk at tau): 'sample' from previously available
      censor marks (any existing censors with time < tau) and *add* new censor rows at time=tau,
      starting at `p` (original donors remain unchanged). Updates `p`.
    - If delta < 0 (need to DECREASE at-risk at tau): remove excess by moving *duplicated*
      censor marks at `tau` (exact matches first) to just before `tau` (bounded by `prev_tau`);
      if still short, also move the earliest censors with time > tau.

    Returns
    -------
    IPD_t, IPD_e, p, n_adjusted
      - n_adjusted signed:
          +n : increased at-risk (added new censors at tau)
          -n : decreased at-risk (moved censors from >= tau to before tau)
    """
    if rng is None:
        rng = np.random.default_rng()

    current = _at_risk_count(IPD_t, tau)
    delta = target_at_risk - current

    if delta == 0:
        return IPD_t, IPD_e, 0

    cens_idx = np.where(IPD_e == 0)[0]

    if delta > 0:
        # Need to DECREASE at-risk at tau by removing existing censors
        need = delta
        moved = 0

        # First find all censors between [prev_tau, tau)
        donors = cens_idx[(tau-interval_length < IPD_t[cens_idx]) & (IPD_t[cens_idx] < tau)]

        # Remove the censors from this interval by moving them to after tau
        if donors.size != 0:
            IPD_t[donors[:need]] = [IPD_t[-1] for _ in range(min(need, donors.size))]
            moved += need

    if delta < 0:
        # Need to increase at-risk at tau by ADDING new censor rows at tau.
        # We "sample" donors from existing censors with time < tau (for bookkeeping),
        # but we do NOT remove the donors; we simply create new rows at tau.
        donors = cens_idx[(tau-interval_length < IPD_t[cens_idx]) & (IPD_t[cens_idx] < tau)]

        if donors.size != 0:
            if debug:
                print(len(donors), "donors available to sample from")
            k = delta
            pick = donors if donors.size == k else rng.choice(donors, size=k, replace=True)
            # Write new censor rows at tau into free slots starting from p
            end = p + k
            IPD_t[p:end] = IPD_t[pick]
            IPD_e[p:end] = 0
            p = end
        else:
            if debug:
                print("No donors available; adding new event at tau without sampling")
            k = delta
            end = p + k
            IPD_t[p:end] = rng.uniform(low=max(tau-interval_length+1e-6, prev_tau-1e-6), size=k)
            IPD_e[p:end] = 1
            p = end

        return IPD_t, IPD_e, -k
    
def get_ipd(
    n: int,
    t: Iterable[float],
    S: Iterable[float],
    cens_t: Optional[Iterable[float]] = None,
    match_tol: float = 5e-3,               # tolerance for |est_S - S_target|
    max_extra_censors_per_bin: int = 100,  # guardrail per bin
    random_state: Optional[int] = None,
    risk_table: Optional[pd.DataFrame] = None,  # for reproducible sampling,
    debug: bool = False
) -> pd.DataFrame:
    """
    Reconstruct individual patient data (IPD) from digitized KM data, allowing
    multiple patients to share the same censor time.

    NEW: When deciding whether to keep the last placed death or revert it, we
    branch two candidates ("after" vs "before"), try adding extra censors within
    the bin for both branches (up to a cap), then pick the branch with the lower
    absolute error to the target survival at the drop time.
    """
    rng = np.random.RandomState(random_state)

    # ---- Preprocessing ----
    t = np.asarray(list(t), dtype=float)  # event (drop) times
    S = np.asarray(list(S), dtype=float)  # survival proportions

    if cens_t is None:
        cens = np.array([], dtype=float)
        has_cens = False
    else:
        cens = np.asarray(list(cens_t), dtype=float)
        cens = cens[~np.isnan(cens)]
        has_cens = cens.size > 0

    # Sorts
    t = np.sort(t) if t.size > 0 else np.array([], dtype=float)
    S = np.sort(S)[::-1] if S.size > 0 else np.array([], dtype=float)
    if has_cens:
        cens = np.sort(cens)

    # Clamp invalids
    if t.size > 0:
        t = np.where(t < 0, 0.0, t)
    if S.size > 0:
        S = np.where(S > 1, 1.0, S)

    # Initialize KM curve with origin
    t = np.concatenate(([0.0], t))
    S = np.concatenate(([1.0], S))

    # Determine max follow-up
    maxFU_candidates = [np.max(t)] if t.size > 0 else []
    if has_cens and cens.size > 0:
        maxFU_candidates.append(np.max(cens))
    maxFU = float(np.max(maxFU_candidates)) if maxFU_candidates else 0.0

    # Prepare risk times if provided
    if risk_table is not None and len(risk_table) > 0:
        rt = risk_table.copy()
        # enforce monotone increasing times
        rt = rt.sort_values("time").reset_index(drop=True)
        risk_times = rt["time"].to_numpy(dtype=float)
        risk_targets = rt["n_at_risk"].to_numpy(dtype=int)
    
        interval_length = risk_times[1] - risk_times[0]

    risk_times = None
    risk_targets = None
    risk_idx = 0

    # Allocate shells
    IPD_t = np.full(n, maxFU, dtype=float)   # default: censored at maxFU
    IPD_e = np.zeros(n, dtype=int)           # default: censored
    p = 0                                    # next free slot

    last_avail_censor_bin = np.array([cens[0]], dtype=float)
    # Main loop over bins between consecutive drops
    for i in range(1, len(t)):
        if debug:
            print(f"Processing bin {i}/{len(t)-1}, p={p}/{n}", 'Target S:', S[i] if i < len(S) else (S[-1] if len(S) > 0 else 1.0))

        t_prev, t_i = t[i - 1], t[i]
        S_target = S[i] if i < len(S) else (S[-1] if len(S) > 0 else 1.0)

        # --- Risk-table alignment: fix all risk times τ with τ <= t_i and not yet aligned ---
        if risk_times is not None:
            while risk_idx < len(risk_times) and risk_times[risk_idx] <= t_i:
                tau = risk_times[risk_idx]
                target = risk_targets[risk_idx]
                prev_tau = 0.0 if risk_idx == 0 else float(risk_times[risk_idx - 1])
                # print(tau, prev_tau)
                output =_align_single_boundary(IPD_t, IPD_e, tau, target, prev_tau,p,interval_length=interval_length,rng=rng, debug=debug)
                IPD_t, IPD_e, n_adjusted = output
                p += n_adjusted
                risk_idx += 1
                if debug:
                    print('After alignment:',_at_risk_count(IPD_t, risk_times[risk_idx-1]) if risk_times is not None and risk_idx>0 else "N/A")

        # Start Processing
        if debug:
            print(f"Processing bin {i}/{len(t)-1} with p={p}/{n} at time {t[i]:.4f}", 'Target S:', S_target)
        # ---- 1) Place known censors in this bin ----
        bin_mask = np.array([], dtype=bool)
        bin_cens_times = np.array([], dtype=float)
        if has_cens:
            bin_mask = (cens >= t_prev) & (cens < t_i)
            bin_cens_times = cens[bin_mask]  # may be empty, may contain duplicates

            if bin_cens_times.size > 0 and p < n:
                n_cens_place = int(min(bin_cens_times.size, n - p))
                if n_cens_place > 0:
                    csel = bin_cens_times[:n_cens_place]
                    IPD_t[p:p + n_cens_place] = csel
                    IPD_e[p:p + n_cens_place] = 0
                    p += n_cens_place
            
            # Update the rolling "last available" bin if this bin has censors
            if bin_cens_times.size > 0:
                last_avail_censor_bin = bin_cens_times

        # ---- 2) Add deaths at the drop to match target survival ----
        n_died = 0
        est_S = _km_survival_at(IPD_t, IPD_e, t_i)
        diff = est_S - S_target

        death_indices_this_bin = []   # absolute indices in IPD arrays

        # Iteratively add deaths until we cross or meet the target
        while (diff > 0.0) and (p + n_died < n):
            idx = p + n_died
            IPD_t[idx] = t_i
            IPD_e[idx] = 1
            death_indices_this_bin.append(idx)
            n_died += 1
            est_S = _km_survival_at(IPD_t, IPD_e, t_i)
            diff = est_S - S_target
            if debug:
                print(f'  Placed death at {t_i:.4f}, new est_S={est_S:.6f}, error = {diff:.6f}')

        # ---- 2.5) Branch-and-compare with extra censors (AFTER: m+1, BEFORE: m, WAY_BEFORE: m-1) ----
        # We now decide among three candidates:
        #   AFTER      -> keep all deaths placed in this bin (n_died = m+1)
        #   BEFORE     -> revert the last placed death          (n_died = m)
        #   WAY_BEFORE -> revert the last TWO placed deaths     (n_died = m-1), if n_died >= 2
        if n_died > 0:
            # ---------- AFTER (m+1) ----------
            IPD_t_after = IPD_t.copy()
            IPD_e_after = IPD_e.copy()
            p_after = p + n_died
            est_after = _km_survival_at(IPD_t_after, IPD_e_after, t_i)
            IPD_t_after2, IPD_e_after2, p_after2, est_after2, extras_after = _try_add_censors_in_bin(
                IPD_t_after, IPD_e_after, p_after, n, S_target, t_prev, t_i,
                bin_cens_times, match_tol, max_extra_censors_per_bin, rng, debug=debug
            )
            err_after = abs(est_after2 - S_target)

            # ---------- BEFORE (m) ----------
            last_idx = death_indices_this_bin[-1]
            IPD_t_before = IPD_t.copy()
            IPD_e_before = IPD_e.copy()
            # revert the last death to a tail censor at maxFU
            IPD_t_before[last_idx] = maxFU
            IPD_e_before[last_idx] = 0
            p_before = p + (n_died - 1)
            est_before = _km_survival_at(IPD_t_before, IPD_e_before, t_i)
            IPD_t_before2, IPD_e_before2, p_before2, est_before2, extras_before = _try_add_censors_in_bin(
                IPD_t_before, IPD_e_before, p_before, n, S_target, t_prev, t_i,
                bin_cens_times, match_tol, max_extra_censors_per_bin, rng, debug=debug
            )
            err_before = abs(est_before2 - S_target)

            # ---------- WAY_BEFORE (m-1), if available ----------
            have_way_before = (n_died >= 2)
            if have_way_before:
                last2_idx = death_indices_this_bin[-2]
                IPD_t_way = IPD_t_before.copy()
                IPD_e_way = IPD_e_before.copy()
                # also revert the second-to-last death
                IPD_t_way[last2_idx] = maxFU
                IPD_e_way[last2_idx] = 0
                p_way = p + (n_died - 2)
                est_way = _km_survival_at(IPD_t_way, IPD_e_way, t_i)
                IPD_t_way2, IPD_e_way2, p_way2, est_way2, extras_way = _try_add_censors_in_bin(
                    IPD_t_way, IPD_e_way, p_way, n, S_target, t_prev, t_i,
                    bin_cens_times, match_tol, max_extra_censors_per_bin, rng, debug=debug
                )
                err_way = abs(est_way2 - S_target)

            # ---------- Choose the best branch ----------
            # Build comparison table
            candidates = [
                ("after",  err_after,  IPD_t_after2,  IPD_e_after2,  p_after2,  est_after2,  0),  # pop 0
                ("before", err_before, IPD_t_before2, IPD_e_before2, p_before2, est_before2, 1),  # pop 1
            ]
            if have_way_before:
                candidates.append(
                    ("way_before", err_way, IPD_t_way2, IPD_e_way2, p_way2, est_way2, 2)        # pop 2
                )

            # Prefer the candidate with the smallest absolute error;
            # tie-breaker: more deaths first (after > before > way_before), you can change this if you prefer.
            # So we'll sort by (error, tie_rank) where AFTER has lowest tie_rank.
            tie_rank = {"after": 0, "before": 1, "way_before": 2}
            best = min(candidates, key=lambda c: (c[1], tie_rank[c[0]]))
            best_name, best_err, best_t, best_e, best_p, best_est, pops = best

            if debug:
                msg = [f"  AFTER:  extras={extras_after},  est_S={est_after2:.6f}, err={err_after:.6f}",
                       f"  BEFORE: extras={extras_before}, est_S={est_before2:.6f}, err={err_before:.6f}"]
                if have_way_before:
                    msg.append(f"  WAY_BEFORE: extras={extras_way}, est_S={est_way2:.6f}, err={err_way:.6f}")
                print("\n".join(msg))
                print(f"  >> Chosen branch: {best_name} (pop {pops}) with err={best_err:.6f}")

            # Apply the chosen branch
            IPD_t = best_t
            IPD_e = best_e
            p = best_p
            est_curr = best_est

            # Update death bookkeeping to reflect how many last deaths were reverted
            for _ in range(pops):
                if death_indices_this_bin:
                    death_indices_this_bin.pop()
                    n_died -= 1

            reverted_last_death = (best_name != "after")

        else:
            # No deaths placed; just attempt extra censors (single branch)
            IPD_t, IPD_e, p, est_curr, _ = _try_add_censors_in_bin(
                IPD_t, IPD_e, p, n, S_target, t_prev, t_i,
                bin_cens_times, match_tol, max_extra_censors_per_bin, rng, debug=debug
            )
            reverted_last_death = False  # nothing was placed
                
    # ---- Tail censors at/after the last drop time ----
    if has_cens and p < n:
        last_drop = np.max(t) if t.size > 0 else 0.0
        tail_mask = cens >= last_drop
        tail_times = cens[tail_mask]
        if tail_times.size > 0:
            n_tail = int(min(tail_times.size, n - p))
            IPD_t[p:p + n_tail] = tail_times[:n_tail]
            IPD_e[p:p + n_tail] = 0
            p += n_tail

    return pd.DataFrame({"time": IPD_t.astype(float), "event": IPD_e.astype(int)})