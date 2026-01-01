# Bivariate Weibull Bayesian Estimation with ABC (MCMC & HMC)

This repository contains MATLAB implementations of **Approximate Bayesian Computation (ABC)**
using **Metropolis–Hastings MCMC** and **Hamiltonian Monte Carlo (HMC)** for a
**5-parameter bivariate Weibull distribution** applied to failure analysis of advanced ceramics.

The code accompanies the article:

> P. Hajizadeh, M. Khosravi, M. Ravandi,  
> *Failure analysis of advanced ceramics using bivariate Weibull distribution and Bayesian estimation*,  
> Engineering with Computers (2025).  
> DOI: https://doi.org/10.1007/s00366-025-02143-x

The goal of this repo is to provide a **reproducible implementation** of the statistical
model and ABC algorithms and to make it easy to compare these results with
NASA’s **CARES/LIFE** physics-based Batdorf model.

---

## 1. Bivariate Weibull model

We use a 5-parameter bivariate Weibull distribution
for two failure variables $x$ and $y$:

- shape parameters: $\beta_1, \beta_2 > 0$  
- scale parameters: $\theta_1, \theta_2 > 0$  
- dependence parameter: $\delta > 0$

The cumulative distribution function (CDF) is:

$$
\Large
F(x, y \mid \theta_1, \theta_2, \beta_1, \beta_2, \delta)= 1 - \exp\left(-\left(\left(\frac{x}{\theta_1}\right)^{\beta_1/\delta}+\left(\frac{y}{\theta_2}\right)^{\beta_2/\delta}\right)^{\delta}\right)
$$

where $x > 0$, $y > 0$.

The joint probability density function $f(x,y \mid \cdot)$ is obtained by
differentiating the CDF with respect to $x$ and $y$ (see Eq. (2) in the paper)
and is used inside the ABC algorithms.

We collect the parameters as

$$
\Large
\theta = (\theta_1, \theta_2, \beta_1, \beta_2, \delta).
$$

---

## 2. Priors

We distinguish between the priors used in the two ABC algorithms.

### 2.1 Priors for ABC–MCMC (Metropolis–Hastings)

For the ABC–MCMC sampler we use **independent exponential priors** for all
five positive parameters

$$
\Large
\theta_1 \sim \mathrm{Exp}(\lambda_1), \quad
\theta_2 \sim \mathrm{Exp}(\lambda_2), \quad
\beta_1  \sim \mathrm{Exp}(\lambda_3), \quad
\beta_2  \sim \mathrm{Exp}(\lambda_4), \quad
\delta   \sim \mathrm{Exp}(\lambda_5).
$$

The rate parameters $\lambda_1,\dots,\lambda_5$ are chosen to be weakly
informative and can be modified directly in the ABC–MCMC MATLAB script.

### 2.2 Priors for ABC–HMC (Stan)

For the ABC–HMC sampler implemented in **Stan**, we use **independent uniform
priors** over physically reasonable ranges:

$$
\Large
\theta_1 \sim \mathrm{Unif}(a_1, b_1), \quad
\theta_2 \sim \mathrm{Unif}(a_2, b_2), \quad
\beta_1  \sim \mathrm{Unif}(a_3, b_3), \quad
\beta_2  \sim \mathrm{Unif}(a_4, b_4), \quad
\delta   \sim \mathrm{Unif}(a_5, b_5).
$$

The bounds $(a_j, b_j)$ are selected based on engineering knowledge of the
material and test geometry and are explicitly coded in the Stan model file
(e.g. `stan/bivariate_weibull_abc_hmc.stan`).

---

## 3. ABC distance and tolerance

Given parameters $\theta$, we simulate synthetic failure pairs
$(x_i^{\text{sim}}, y_i^{\text{sim}})$ under the bivariate Weibull model and
compare them with the observed failure data $(x_i^{\text{obs}}, y_i^{\text{obs}})$.

We use a simple Euclidean distance in the 2D failure space (or on suitable
summaries, depending on the script):

$$
\Large
d(\theta)=
\left\|
\begin{pmatrix}
x^{\text{sim}} \\
y^{\text{sim}}
\end{pmatrix}-
\begin{pmatrix}
x^{\text{obs}} \\
y^{\text{obs}}
\end{pmatrix}
\right\|_2
$$

A parameter proposal is **accepted by ABC** if

$$
\Large
d(\theta) \le \varepsilon,
$$

where $\varepsilon$ is an ABC tolerance. In practice, the scripts allow
you to tune $\varepsilon$ to balance accuracy vs. computational effort.

---

## 4. ABC–MCMC (Metropolis–Hastings)

ABC–MCMC runs a Metropolis–Hastings chain on

$$
\Large
\theta = (\theta_1, \theta_2, \beta_1, \beta_2, \delta).
$$

At each iteration:

1. **Propose new parameters**  
   Draw a proposal $\theta^{\*}$ using component-wise positive proposals
   (e.g. Gamma random walks) centered around the current $\theta$.

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
   Simulate synthetic data under $\theta^{\*}$, compute $d(\theta^{\*})$,
   and reject immediately if $d(\theta^{\*}) > \varepsilon$.

4. **Accept / reject**  
   If $d(\theta^{\*}) \le \varepsilon$, accept $\theta^{\*}$ with probability
   $\alpha_{\mathrm{MH}}$; otherwise keep $\theta$.

After burn-in, the chain provides ABC posterior samples for $\theta$.
