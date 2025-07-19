### `header_highdim.RData`
RData file containing regression design matrix, response variable, and hyper-parameters for Bayesian linear regression model.

---

### `distou_revise.cpp`
C++ source file implementing the Discrete Over-Relaxation algorithm.

---

### `sparse_linear_algos.cpp`
C++ script implementing various discrete samplers, with dependencies on `distou_revise.cpp`.

---

### `repeated_runs.R`
Runs repeated MCMC chains for each sampler to replicate and evaluate results. Results are saved as `.RData` files, with filenames indicating the sampler used.

---

###`PIP_calculation.R`
R script for calculating metrics of PIP.

---

### `plots`
Notebook that generates all figures and tables for analysis and reporting.