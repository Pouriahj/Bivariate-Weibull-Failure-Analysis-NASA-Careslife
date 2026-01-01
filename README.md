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

The cumulative distribution function (CDF) is


$$
\Large
F(x, y; \theta_{1}, \theta_{2}, \beta_{1}, \beta_{2}, \delta)
=
1 - \exp\left\{
-
\left[
\left( \frac{x}{\theta_1} \right)^{\beta_1 / \delta}
+
\left( \frac{y}{\theta_2} \right)^{\beta_2 / \delta}
\right]^{\delta}
\right\}
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

All parameters are **strictly positive**, so we assign independent
exponential priors

$$
\Large
\theta_1 \sim \mathrm{Exp}(\lambda_1), \quad
\theta_2 \sim \mathrm{Exp}(\lambda_2), \quad
\beta_1  \sim \mathrm{Exp}(\lambda_3), \quad
\beta_2  \sim \mathrm{Exp}(\lambda_4), \quad
\delta   \sim \mathrm{Exp}(\lambda_5).
$$

The hyperparameters $\lambda_1, \ldots, \lambda_5$ are chosen to be weakly
informative and can be edited directly in the MATLAB scripts.

---

## 3. ABC distance and tolerance

Given parameters $\theta$, we simulate synthetic failure pairs
$(x_i^{\text{sim}}, y_i^{\text{sim}})$ under the bivariate Weibull model and
compare them with the observed failure data $(x_i^{\text{obs}}, y_i^{\text{obs}})$.

We use a simple Euclidean distance in the 2D failure space (or on suitable
summaries, depending on the script):

$$
\Large
d(\theta)
=
\left\|
\begin{pmatrix}
x^{\text{sim}} \\
y^{\text{sim}}
\end{pmatrix}
-
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
   \alpha_{\mathrm{MH}}
   =
   \min\Biggl(
     1,\;
     \frac{\pi(\theta^{\*})}{\pi(\theta)}
     \prod_{i=1}^{N}
     \frac{
       f\!\left(x_i^{\text{obs}}, y_i^{\text{obs}} \mid \theta^{\*}\right)
     }{
       f\!\left(x_i^{\text{obs}}, y_i^{\text{obs}} \mid \theta\right)
     }
     \frac{
       q(\theta \mid \theta^{\*})
     }{
       q(\theta^{\*} \mid \theta)
     }
   \Biggr)
   $$

   where $q$ is the proposal density.

3. **ABC filter**  
   Simulate synthetic data under $\theta^{\*}$, compute $d(\theta^{\*})$,
   and reject immediately if $d(\theta^{\*}) > \varepsilon$.

4. **Accept / reject**  
   If $d(\theta^{\*}) \le \varepsilon$, accept $\theta^{\*}$ with probability
   $\alpha_{\mathrm{MH}}$; otherwise keep $\theta$.

After burn-in, the chain provides ABC posterior samples for $\theta$.

**Flowchart (ABC–MCMC)**

```markdown
![Flowchart of the ABC–MCMC algorithm](figures/abc_mcmc_bivariate.png)
