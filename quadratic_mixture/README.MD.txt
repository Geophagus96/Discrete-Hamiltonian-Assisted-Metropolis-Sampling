### `probparalleltest.R`
R script for generating a **quadratic mixture distribution**.  
Also constructs probability tables for the exact distribution and saves them as `probtable.RData`.

---

### `distou_revise.cpp`
C++ source file implementing the Discrete Over-Relaxation algorithm.

---

### `algos.R`
R script implementing various discrete samplers, with dependencies on `distou_revise.cpp`.

---

### `repeated_runs.R`
Runs repeated MCMC chains for each sampler to replicate and evaluate results. Results are saved as `.RData` files, with filenames indicating the sampler used.

---

### `plot_all`
Notebook that generates all figures and tables for analysis and reporting.