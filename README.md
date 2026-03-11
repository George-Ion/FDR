# Computational Statistics & Multiple Testing

This repository contains a suite of algorithms and simulations focused on modern computational statistics, specifically addressing **False Discovery Rate (FDR)**, **Maximum Likelihood Estimation**, and **Stochastic Sampling**.

## Core Statistical Frameworks

### 1. False Discovery Rate (FDR) Control
When testing $m$ hypotheses, the **Benjamini-Hochberg (BH)** procedure is implemented to control the FDR. The algorithm identifies the largest $k$ such that:

$$P_{(k)} \le \frac{k}{m} \alpha$$

All hypotheses $H_{(i)}$ for $i = 1, \dots, k$ are rejected, ensuring that $E[V / R] \le \alpha$.


### 2. The EM Algorithm (`EM_algorithm.R`)
Used for parameter estimation in models with latent variables (e.g., Gaussian Mixture Models). The algorithm alternates between two steps:

* **E-Step**: Compute the expected value of the log-likelihood:
    $$Q(\theta | \theta^{(t)}) = E_{Z|X,\theta^{(t)}} [\log L(\theta; X, Z)]$$
* **M-Step**: Update parameters by maximizing $Q$:
    $$\theta^{(t+1)} = \text{arg max}_{\theta} Q(\theta | \theta^{(t)})$$


### 3. Metropolis-Hastings MCMC (`metropolis_hastings.R`)
To sample from a target distribution $P(x)$ where the normalization constant is unknown, we use a proposal distribution $q(x' | x)$. The acceptance probability $\alpha$ is defined as:

$$\alpha(x, x') = \min \left( 1, \frac{P(x') q(x|x')}{P(x) q(x'|x)} \right)$$


---

## Implementation Matrix

| Method | Script | Application | Computational Constraint |
| :--- | :--- | :--- | :--- |
| **FDR** | `fdr_analysis.R` | Multiple Testing | Dependency between tests |
| **EM** | `EM_algorithm.R` | Latent Variables | Convergence to local optima |
| **MCMC** | `metropolis_hastings.R` | Bayesian Inference | High autocorrelation |
| **Bootstrap** | `bootstrap.R` | Uncertainty | Sample size $n$ and iterations $B$ |

---

##  Usage
Each script is written in **R** and includes a simulated dataset to demonstrate convergence. To run the simulations:

```R
source("metropolis_hastings.R")
# Check trace plots for convergence
