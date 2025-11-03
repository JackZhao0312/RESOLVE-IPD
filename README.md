# RESOLVE-IPD

This repository provides practical utilities and examples to turn digitized KM survival data into per-patient IPD suitable for downstream analysis, simulation, and reproducible research. It is geared toward scenarios where published KM curves and (optionally) risk tables are available, but raw patient-level data are not.

- Language: Python / R
- License: MIT

## Highlights

- IPD reconstruction from digitized KM survival points:
  - Accepts event (drop) times and survival probabilities from the KM curve.
  - Optionally accepts observed censor times.
- Flexible censor handling:
  - Allows multiple individuals to share the same censor time.
  - Adds additional censors within bins to better match target survival within a tolerance.
- Branch-and-compare strategy:
  - Evaluates multiple candidate death placements and nearby censor additions to minimize error to the target KM survival at each drop time.
- Optional risk-table alignment (experimental):
  - Aligns the number-at-risk at specified time points by relocating/adding censors around risk table boundaries.
- Reproducibility:
  - Deterministic results with `random_state`.

## Repository layout

- [CEN_KM.py](https://github.com/JackZhao0312/RESOLVE-IPD/blob/main/CEN_KM.py) — Core implementation, including:
  - `_km_survival_at(...)`: KM survival at a given time.
  - `_try_add_censors_in_bin(...)`: Iteratively add censors in a time bin to reduce error to target survival.
  - `_align_single_boundary(...)`: Experimental at-risk alignment around risk-table times.
  - `get_ipd(...)`: Main API to reconstruct IPD as a DataFrame with `time` and `event` columns.
- [CEN-KM Example.ipynb](https://github.com/JackZhao0312/RESOLVE-IPD/blob/main/CEN-KM%20Example.ipynb) — End-to-end example notebook for the CEN-KM workflow.
- [VEC-KM Example.ipynb](https://github.com/JackZhao0312/RESOLVE-IPD/blob/main/VEC-KM%20Example.ipynb) — Example notebook for an alternative KM setting/workflow.
- [MAPLE](https://github.com/JackZhao0312/RESOLVE-IPD/tree/main/MAPLE) — Plausible subgroup label recoverying algorithm.
- [LICENSE](https://github.com/JackZhao0312/RESOLVE-IPD/blob/main/LICENSE) — MIT License.

## Installation

This project uses standard scientific Python dependencies.

- Python 3.9+ recommended
- Dependencies:
  - numpy
  - pandas
  - lifelines
  - jupyter (for running notebooks)
  - matplotlib or plotnine (optional, for plotting)
  - PyMuPDF (for loading and digitizing PDF)

## Quick start (Python API)

Below is a minimal example using the core `get_ipd` function in [CEN_KM.py](https://github.com/JackZhao0312/RESOLVE-IPD/blob/main/CEN_KM.py).

```python
import numpy as np
import pandas as pd
from CEN_KM import get_ipd
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt

# Example digitized KM data (replace with your own)
# t: event/drop times where the KM curve steps down
# S: corresponding survival proportions after each drop time
t = [2.0, 4.0, 6.0, 8.0, 12.0]             # event/drop times
S = [0.92, 0.85, 0.78, 0.70, 0.60]         # survival levels
cens_t = [3.1, 5.4, 7.9, 10.2, 11.8, 11.8] # optional censor times (may contain duplicates)

# Optional risk table (experimental support)
# Provide n_at_risk at specific times to guide alignment
risk_df = pd.DataFrame({
    "time": [0.0, 6.0, 12.0],
    "n_at_risk": [100, 72, 45],
})

# Target cohort size
n = 100

ipd = get_ipd(
    n=n,
    t=t,
    S=S,
    cens_t=cens_t,
    match_tol=5e-3,
    max_extra_censors_per_bin=100,
    random_state=42,
    risk_table=risk_df,   # optional; experimental
    debug=False
)

print(ipd.head())
# Columns: time (float), event (int; 1=death, 0=censored)

# Validate by refitting KM on reconstructed IPD
kmf = KaplanMeierFitter()
kmf.fit(ipd["time"], event_observed=ipd["event"], label="Reconstructed KM")
kmf.plot(ci_show=False)
plt.title("Reconstructed KM")
plt.show()
```

## Examples

- Run the notebooks:
  - [CEN-KM Example.ipynb](https://github.com/JackZhao0312/RESOLVE-IPD/blob/main/CEN-KM%20Example.ipynb)
  - [VEC-KM Example.ipynb](https://github.com/JackZhao0312/RESOLVE-IPD/blob/main/VEC-KM%20Example.ipynb)


## Reproducibility

- Always set `random_state` when comparing across runs or reporting results.
- Record the exact inputs (`t`, `S`, `cens_t`, risk table) used for reconstruction.

## Contributing

Contributions are welcome!
- Open an issue for bugs, questions, or feature requests.
- Submit pull requests with clear descriptions and tests or examples where applicable.

## License

This project is licensed under the [MIT License](https://github.com/JackZhao0312/RESOLVE-IPD/blob/main/LICENSE).
