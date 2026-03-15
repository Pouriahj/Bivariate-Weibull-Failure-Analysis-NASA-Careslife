# Bivariate Weibull Failure Analysis of Advanced Ceramics (NASA CARES/LIFE Benchmark)

This repository contains implementations of Bayesian inference for a **5-parameter bivariate Weibull distribution** applied to failure analysis of advanced ceramics, calibrated to experimental data and benchmarked against **NASA’s CARES/LIFE** Batdorf-based reliability predictions.

The repository includes:

- a **MATLAB implementation** of **Approximate Bayesian Computation (ABC)** using **Metropolis–Hastings (MH)**,
- and an **R/Stan implementation** of **Hamiltonian Monte Carlo (HMC)**.

The code accompanies the article:

> P. Hajizadeh, M. Khosravi, M. Ravandi,  
> *Failure analysis of advanced ceramics using bivariate Weibull distribution and Bayesian estimation*,  
> Engineering with Computers (2025).  
> DOI: https://doi.org/10.1007/s00366-025-02143-x

The goal of this repository is to provide a **reproducible implementation** of the statistical model and Bayesian algorithms, and to make it easy to compare these results with NASA’s **CARES/LIFE** physics-based Batdorf model.

---

## Code Organization

This repository contains multiple implementations of Bayesian inference for the bivariate Weibull failure model of advanced ceramics. The code is organized by language and inference method to keep the workflow clear and avoid unnecessary clutter.

### MATLAB: ABC–Metropolis–Hastings (MH)

The `matlab/mh/` folder contains the MATLAB implementation of the **ABC–MH** workflow.

This part of the code is used for:

- defining the bivariate Weibull model,
- generating proposal values in the MH sampler,
- evaluating the discrepancy/distance criterion in the ABC framework,
- running posterior sampling

This implementation is the main **MATLAB-based likelihood-free Bayesian workflow** in the repository.

### R + Stan: Hamiltonian Monte Carlo (HMC)

The **Hamiltonian Monte Carlo (HMC)** implementation is split into two parts:

- `r/hmc_rstan/` contains the **R scripts**
- `stan/` contains the **Stan model file**

The R scripts are responsible for:

- preparing the data,
- calling the Stan model through **RStan**,
- running posterior sampling,

The Stan file defines the probabilistic model structure and enables efficient gradient-based Bayesian inference using HMC.

This part of the repository is the main **gradient-based Bayesian workflow** used for comparison with the MATLAB ABC–MH implementation.

### External dependency for random bivariate data generation

For generating random bivariate data, this repository uses functionality based on **Chebfun**.

Please install Chebfun separately before running the corresponding data-generation scripts.

- GitHub: https://github.com/chebfun/chebfun
- MATLAB: https://www.mathworks.com/matlabcentral/fileexchange/47023-chebfun-current-version

---

## 1. Bivariate Weibull model

We use a 5-parameter bivariate Weibull distribution for two failure variables $x$ and $y$:

- shape parameters: $\beta_1, \beta_2 > 0$
- scale parameters: $\theta_1, \theta_2 > 0$
- dependence parameter: $\delta > 0$

The cumulative distribution function (CDF) is:

$$
\Large
F(x, y \mid \theta_1, \theta_2, \beta_1, \beta_2, \delta)
= 1 - \exp\left(
  -\left(
  \left(\frac{x}{\theta_1}\right)^{\beta_1/\delta}
  +
  \left(\frac{y}{\theta_2}\right)^{\beta_2/\delta}
  \right)^{\delta}
  \right)
$$

where $x > 0$, $y > 0$.

The joint probability density function $f(x,y \mid \cdot)$ is obtained by differentiating the CDF with respect to $x$ and $y$ (see Eq. (2) in the paper) and is used inside the ABC algorithms.

We collect the parameters as

$$
\Large
\theta = (\theta_1, \theta_2, \beta_1, \beta_2, \delta).
$$

---

## 2. Priors

We distinguish between the priors used in the two Bayesian algorithms.

### 2.1 Priors for ABC–MCMC (Metropolis–Hastings)

For the ABC–MCMC sampler we use **independent exponential priors** for all five positive parameters:

$$
\Large
\theta_1 \sim \mathrm{Exp}(\lambda_1), \quad
\theta_2 \sim \mathrm{Exp}(\lambda_2), \quad
\beta_1  \sim \mathrm{Exp}(\lambda_3), \quad
\beta_2  \sim \mathrm{Exp}(\lambda_4), \quad
\delta   \sim \mathrm{Exp}(\lambda_5).
$$

The rate parameters $\lambda_1,\dots,\lambda_5$ are chosen to be weakly informative and can be modified directly in the ABC–MCMC MATLAB script.

### 2.2 Priors for ABC–HMC (Stan)

For the ABC–HMC sampler implemented in **Stan**, we use **independent uniform priors** over physically reasonable ranges:

$$
\Large
\theta_1 \sim \mathrm{Unif}(a_1, b_1), \quad
\theta_2 \sim \mathrm{Unif}(a_2, b_2), \quad
\beta_1  \sim \mathrm{Unif}(a_3, b_3), \quad
\beta_2  \sim \mathrm{Unif}(a_4, b_4), \quad
\delta   \sim \mathrm{Unif}(a_5, b_5).
$$

The bounds $(a_j, b_j)$ are selected based on engineering knowledge of the material and test geometry and are explicitly coded in the Stan model file (e.g. `stan/bivariate_weibull_abc_hmc.stan`).

---

## 3. ABC distance and tolerance

Given a parameter vector $\theta$, we generate synthetic bivariate data
$(x', y')$ from the bivariate Weibull model and compare it with the observed
dataset $(x, y)$.

We write the observed data as

$$
\Large
x =
\begin{pmatrix}
x_1 \\
x_2 \\
\vdots \\
x_n
\end{pmatrix},
\qquad
y =
\begin{pmatrix}
y_1 \\
y_2 \\
\vdots \\
y_n
\end{pmatrix},
$$

and the generated data as

$$
\Large
x' =
\begin{pmatrix}
x'_1 \\
x'_2 \\
\vdots \\
x'_n
\end{pmatrix},
\qquad
y' =
\begin{pmatrix}
y'_1 \\
y'_2 \\
\vdots \\
y'_n
\end{pmatrix},
$$

where $(x', y')$ are generated under the candidate parameter vector $\theta$.

The ABC discrepancy is defined through the norm of the difference between the
observed and generated paired datasets:

$$
\Large
d(\theta)
=\left\|
\begin{pmatrix}
x_1 & y_1 \\
x_2 & y_2 \\
\vdots & \vdots \\
x_n & y_n
\end{pmatrix}
-\begin{pmatrix}
x'_1 & y'_1 \\
x'_2 & y'_2 \\
\vdots & \vdots \\
x'_n & y'_n
\end{pmatrix}
\right\|.
$$

A candidate parameter is accepted by the ABC criterion if

$$
d(\theta) \le \varepsilon,
$$

where $\varepsilon$ is the ABC tolerance. Smaller values of $\varepsilon$
enforce closer agreement between the simulated and observed bivariate datasets,
but they also reduce the acceptance rate and increase computational cost.

---

## 4. ABC–MCMC (Metropolis–Hastings)

ABC–MCMC runs a Metropolis–Hastings chain on

$$
\Large
\theta = (\theta_1, \theta_2, \beta_1, \beta_2, \delta).
$$

At each iteration:

1. **Propose new parameters**  
   Draw a proposal $\theta^{\*}$ using component-wise positive proposals (e.g. Gamma random walks) centered around the current $\theta$.

2. **Compute Metropolis–Hastings ratio**  
   Combine the prior and bivariate Weibull likelihood:

$$
\Large
\alpha_{\mathrm{MH}} =
\min\left(
  1,\,
  \frac{\pi(\theta^{\*})}{\pi(\theta)}
  \prod_{i=1}^{N}
  \frac{
    f\left(x_i^{\text{obs}}, y_i^{\text{obs}} \mid \theta^{\*}\right)
  }{
    f\left(x_i^{\text{obs}}, y_i^{\text{obs}} \mid \theta\right)
  }
  \frac{
    q(\theta \mid \theta^{\*})
  }{
    q(\theta^{\*} \mid \theta)
  }
\right)
$$

where $q$ is the proposal density.

3. **ABC filter**  
   Simulate synthetic data under $\theta^{\*}$, compute $d(\theta^{\*})$, and reject immediately if $d(\theta^{\*}) > \varepsilon$.

4. **Accept / reject**  
   If $d(\theta^{\*}) \le \varepsilon$, accept $\theta^{\*}$ with probability $\alpha_{\mathrm{MH}}$; otherwise keep $\theta$.

After burn-in, the chain provides ABC posterior samples for $\theta$.

---

## 5. ABC–HMC (Hamiltonian Monte Carlo)

ABC–HMC uses Hamiltonian dynamics to explore the posterior distribution of $\theta$ more efficiently than standard random-walk MCMC methods.

In this study, the HMC workflow is implemented in **Stan** using **independent
uniform priors** for all five parameters. The method introduces an auxiliary
momentum variable $m$ with Gaussian distribution

$$
\Large
m \sim \mathcal{N}(0, M),
$$

where $M$ is the mass matrix.

The HMC algorithm proceeds by alternating between two steps:

1. **Randomize the momentum** by sampling
   $m \sim \mathcal{N}(0, M)$.

2. **Simulate the Hamiltonian trajectory** of the system $(\theta, m)$ by
   numerically integrating Hamilton’s equations:

$$
\Large
\frac{d\theta}{dt}=
\frac{\partial H}{\partial m}=
M^{-1} m
$$

and

$$
\Large
\frac{dm}{dt}=
-\frac{\partial H}{\partial \theta}=
\frac{\nabla p(\theta \mid x, y)}{p(\theta \mid x, y)}.
$$

The numerical integration is carried out using the **leapfrog integrator**,
which depends on two tuning parameters:

- the number of leapfrog steps $L$
- the step size $\varepsilon$

Starting from the current state $(\theta^t, m^t)$, the leapfrog updates are

$$
\Large
m^{t+\varepsilon/2}
\leftarrow
m^t
-\frac{\varepsilon}{2}
\frac{\partial U(\theta^t)}{\partial \theta}
$$

$$
\Large
\theta^{t+\varepsilon}
\leftarrow
\theta^t-
\varepsilon M^{-1} m^{t+\varepsilon/2}
$$

$$
\Large
m^{t+\varepsilon}
\leftarrow
m^{t+\varepsilon/2}
-\frac{\varepsilon}{2}
\frac{\partial U(\theta^{t+\varepsilon})}{\partial \theta},
$$

repeated for $t = 0 : L-1$, leading to a proposed state
$(\theta^\*, m^\*)$.

After the simulation, the proposed state is accepted with probability

$$
\Large
\alpha
=\min\left(
1,\;
\frac{p(\theta^\* \mid x, y)}{p(\theta \mid x, y)}
\right),
$$

following the Metropolis–Hastings acceptance rule.

The implementation in this repository
uses the **No-U-Turn Sampler (NUTS)** through **Stan**.

---

## Suggested repository layout

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
│   └── *.stan               # Stan model for HMC
│
├── data/
│   ├── raw/
│   └── processed/
│
├── figures/                 # Flowcharts, benchmark plots, comparison figures
└── docs/                    # Optional notes, methodology, references
