## Code Organization

This repository contains multiple implementations of Bayesian inference for the bivariate Weibull failure model of advanced ceramics. The code is organized by language and inference method to keep the workflow clear and avoid unnecessary clutter.

### MATLAB: ABC–Metropolis–Hastings (MH)

The `matlab/mh/` folder contains the MATLAB implementation of the **Approximate Bayesian Computation – Metropolis–Hastings (ABC–MH)** workflow.

This part of the code is used for:

- defining the bivariate Weibull model parameters,
- generating proposal values in the MH sampler,
- evaluating the discrepancy / distance criterion in the ABC framework,
- running posterior sampling,
- and post-processing parameter estimates.

This implementation is the main **MATLAB-based likelihood-free Bayesian workflow** in the repository.

---

### R + Stan: Hamiltonian Monte Carlo (HMC)

The **Hamiltonian Monte Carlo (HMC)** implementation is split into two parts:

- `r/hmc_rstan/` contains the **R scripts**
- `stan/` contains the **Stan model file**

The R scripts are responsible for:

- preparing the data,
- calling the Stan model through **RStan**,
- running posterior sampling,
- and post-processing / visualizing the results.

The Stan file defines the probabilistic model structure and enables efficient gradient-based Bayesian inference using HMC.

This part of the repository is the main **gradient-based Bayesian workflow** used for comparison with the MATLAB ABC–MH implementation.

---

### Synthetic bivariate data generation

The `matlab/data_generation/` folder contains scripts for generating synthetic bivariate data, including:

- **strength data**
- **time-to-failure data**

These scripts are useful for validating the inference workflow, testing the samplers, and reproducing benchmark examples.

---

### External dependency for random bivariate data generation

For generating random bivariate data, this repository uses functionality based on **Chebfun**.

Please install Chebfun separately before running the corresponding data-generation scripts.

- GitHub: https://github.com/chebfun/chebfun
- MATLAB File Exchange: https://www.mathworks.com/matlabcentral/fileexchange/47023-chebfun-current-version

Chebfun is **not included** in this repository.

---

### Suggested repository layout

```text
.
├── matlab/
│   ├── mh/                  # MATLAB implementation of ABC–MH
│   └── data_generation/     # Synthetic bivariate data generation
│
├── r/
│   └── hmc_rstan/           # R scripts for HMC workflow
│
├── stan/
│   └── *.stan              # Stan model for HMC
│
├── data/
│   ├── raw/
│   └── processed/
│
├── figures/                # Flowcharts, benchmark plots, comparison figures
└── docs/                   # Optional notes, methodology, references
