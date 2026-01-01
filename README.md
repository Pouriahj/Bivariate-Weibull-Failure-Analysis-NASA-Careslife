# Bivariate-Weibull-Approximate-Bayesian-Computation-ABC-

This repository accompanies the article

P. Hajizadeh, M. Khosravi & M. Ravandi
Failure analysis of advanced ceramics using bivariate Weibull distribution and Bayesian estimation, Engineering with Computers (2025).

The paper proposes a 5‑parameter bivariate Weibull distribution to model the joint behaviour of failure strength and time to failure (TTF) of advanced ceramics.  Approximate Bayesian Computation (ABC) is used to infer the parameters when the likelihood is intractable.  Two Markov chain Monte Carlo samplers are explored:
	•	ABC–MCMC (Metropolis–Hastings) – parameters are proposed from Gamma distributions and accepted using the Metropolis–Hastings ratio subject to an ABC distance threshold.
	•	ABC–HMC (Hamiltonian Monte Carlo) – an auxiliary momentum is introduced, Hamiltonian dynamics are simulated via a leap‑frog integrator and proposals are accepted with a standard HMC acceptance probability.  The paper implements HMC via the No‑U‑Turn Sampler (NUTS).

The goal of this repository is to provide Matlab and Python code for reproducing the ABC–MCMC and ABC–HMC analyses, along with documentation and instructions.

⸻

1  Model: Bivariate Weibull distribution

A bivariate Weibull random vector $(X,Y)$ with scale parameters $\theta_1,\theta_2$, shape parameters $\beta_1,\beta_2$ and association parameter $\delta\ge 0$ has joint cumulative distribution function (CDF)

$$
\Large
F(x,y\mid \theta_1,\theta_2,\beta_1,\beta_2,\delta)
= 1 - \exp\Biggl(
	•	\Bigl[\Bigl(\tfrac{x}{\theta_1}\Bigr)^{\beta_1 \delta}
+ \Bigl(\tfrac{y}{\theta_2}\Bigr)^{\beta_2 \delta}\Bigr]^{\delta}
\Biggr).
$$

The probability density function (PDF) can be written:

$$
\begin{aligned}
\Large
f(x,y\mid \theta) &= \prod_{i=1}^{n}
\frac{\beta_1}{\theta_1}
\Bigl(\tfrac{x_i}{\theta_1}\Bigr)^{\beta_1/\delta-1}
\frac{\beta_2}{\theta_2}
\Bigl(\tfrac{y_i}{\theta_2}\Bigr)^{\beta_2/\delta-1}
\Bigl[\bigl(\tfrac{x_i}{\theta_1}\bigr)^{\beta_1 \delta}+\bigl(\tfrac{y_i}{\theta_2}\bigr)^{\beta_2 \delta}\Bigr]^{\delta-2}

\exp\Bigl{-\Bigl[\bigl(\tfrac{x_i}{\theta_1}\bigr)^{\beta_1 \delta}+\bigl(\tfrac{y_i}{\theta_2}\bigr)^{\beta_2 \delta}\Bigr]^\delta\Bigr}!
\end{aligned}
$$

where $\theta=(\theta_1,\theta_2,\beta_1,\beta_2,\delta)^T$ is the parameter vector and $(x_i,y_i)$ are observed data points.  The goal is to infer $\theta$ given paired data $(x,y)$.

Priors

Different priors are used depending on the sampler:
	•	ABC–MCMC: Independent exponential priors are assigned to each parameter $\theta_i$, i.e. $\theta_i\sim \mathrm{Exp}(\lambda_i)$.
	•	ABC–HMC: Weakly informative uniform priors $\theta_i \sim \mathrm{Uniform}(a_i, b_i)$ are used.

These priors encode positivity and broad parameter ranges.

Posterior

Bayesian inference uses Bayes’ rule to combine the likelihood and prior:

$$
\Large
p(\theta\mid x,y) \propto f(x,y\mid\theta),\pi(\theta).
$$

Direct sampling from the posterior is difficult, so Approximate Bayesian Computation (ABC) is adopted.  ABC replaces the likelihood with a discrepancy between simulated and observed data.

⸻

2  ABC distance and tolerance

Given a candidate parameter set $\theta$, the ABC procedure simulates data $\bigl(\tilde{x},\tilde{y}\bigr)$ from the bivariate Weibull model and compares it to the observed data $(x,y)$ using a norm:

$$
\Large
d(\theta) = | (\tilde{x},\tilde{y}) - (x,y) |_2.
$$

A proposal is accepted if $d(\theta) \le \varepsilon$, where $\varepsilon$ is a user‑defined tolerance.  Smaller tolerances yield more accurate but less efficient samplers.  When the tolerance is a single value (ABC–MCMC) or a decreasing sequence (ABC–HMC with NUTS), the accepted parameters approximate the posterior distribution.

⸻

3  ABC–MCMC (Metropolis–Hastings)

The ABC–MCMC sampler generates a Markov chain whose stationary distribution approximates the ABC posterior.  For each parameter component $\theta_j$:
	1.	Proposal: Draw $\theta_j^*$ from a Gamma distribution centred at the current value; this choice preserves positivity and matches the scale of the bivariate Weibull posterior.
	2.	Compute Metropolis–Hastings ratio:
$$
\Large
\alpha_{\mathrm{MH}} = \min\Biggl(1, \frac{\pi(\theta^)}{\pi(\theta)}
\prod_{i=1}^{n}
\frac{f(x_i,y_i\mid \theta^)}{f(x_i,y_i\mid \theta)}
\times \frac{q(\theta_j\mid\theta_j^)}{q(\theta_j^\mid\theta_j)}
\Biggr),
$$
where $q(\cdot)$ is the proposal density.
	3.	ABC filter: Simulate data under $\theta^$ and compute $d(\theta^)$; accept the proposal with probability $\alpha_{\mathrm{MH}}$ only if $d(\theta^*) \le \varepsilon$.
	4.	Iterate through all parameter components.  After a burn‑in period the chain samples from the ABC posterior.

Flowchart (ABC–MCMC)

Add your own flowchart illustrating this algorithm, similar to Fig. 4 in the paper. For example, include boxes for draw proposal, simulate data, compute distance, Metropolis–Hastings acceptance, and update parameter.

⸻

4  ABC–HMC (Hamiltonian Monte Carlo)

Hamiltonian Monte Carlo augments the parameter vector $\theta$ with a momentum vector $m$ of the same dimension and defines a Hamiltonian

$$
\Large
H(\theta, m) = U(\theta) + K(m) = -\log p(\theta\mid x,y) + \tfrac{1}{2}, m^T M^{-1} m,
$$

where $M$ is a positive‑definite mass matrix.  The sampler proceeds as follows:
	1.	Initialize $(\theta,m)$; draw $m \sim \mathcal{N}(0, M)$.
	2.	Leapfrog integration: Integrate Hamilton’s equations

$$
\Large
\frac{d\theta}{dt} = M^{-1} m, \quad
\frac{dm}{dt} = \nabla_{\theta} \log p(\theta\mid x,y)
$$

for $L$ steps with step size $\varepsilon$ to propose $(\theta^*, m^*)$.
3. ABC filter: Simulate synthetic data under $\theta^$ and compute $d(\theta^)$; accept the proposed state with probability

$$
\Large
\alpha_{\mathrm{HMC}} = \min\biggl(1, \frac{p(\theta^\mid x,y),\mathcal{N}(m^;0,M)}{p(\theta\mid x,y),\mathcal{N}(m;0,M)}\biggr),
$$

subject to $d(\theta^*) \le \varepsilon$.

The paper implements HMC using the No‑U‑Turn Sampler (NUTS) algorithm, which automatically adapts the number of leap‑frog steps and step sizes for efficient exploration.

Flowchart (ABC–HMC)

Create a flowchart showing the introduction of momentum variables, leapfrog integration steps, ABC distance check, and acceptance step.  For example, start from sample momentum, simulate Hamiltonian trajectory, compute distance and acceptance, update.

⸻

5  Repository structure (suggested)

To mirror the structure of the previous repository, organise your project as follows:
