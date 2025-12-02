# A Bayesian survival model induced by hurdle zero-modified power series discrete frailty with dispersion: an application in lung cancer

This repository contains the code accompanying the paper:

**“A Bayesian Survival Model Induced by Hurdle Zero-Modified Power Series Discrete Frailty with Dispersion: An Application in Lung Cancer”**  
*Molina, K.C., Martínez-Minaya, J., Alvares, D., & Tomazella, V.D. (2025).*  
https://arxiv.org/abs/2505.23568
It provides reproducible material for:

- A simulation study evaluating the HZMGP-based frailty model.
- An application to a large lung cancer dataset from São Paulo.

---

## Repository Structure

```
.
├── zmpg_model.stan          # Stan implementation of the HZMGP survival model
├── simulation.R             # Simulation study (Section 5)
├── application.R            # Lung cancer analysis (Section 6)
└── README.md                # This file
```

---

## Model Summary

The repository implements the Bayesian survival model described in the manuscript:

- Frailty term: **Hurdle Zero-Modified Generalized Poisson (HZMGP)** with dispersion.
- Individual parameters:
  - Zero-modification probability \( \\omega_i \\ )
  - Mean \( \\mu_i \\ )
- Weibull baseline hazard.
- Bayesian inference via Hamiltonian Monte Carlo (NUTS).

This structure enables:

- Flexible frailty types (zero-inflated, zero-deflated, zero-truncated, or standard).
- Posterior-based classification of individuals following the decision rule in Section 4.4.
- Identification of long-term survivors through the discrete frailty component.

---

## Simulation Study

`simulation.R` reproduces Section 5:

- Two simulation scenarios with different patterns of zero-modification and dispersion.
- Sample sizes: 1,000; 10,000; 20,000.
- Outputs posterior means, standard deviations, coverage probabilities, and zero-modification classifications.

---

## Application: Lung Cancer Survival Analysis

`application.R` reproduces Section 6:

- Data: 30,900 lung cancer patients from São Paulo.
- Covariates: age, gender, clinical stage, surgery, radiotherapy, chemotherapy.
- Outputs posterior summaries and classification into ZIGP, GP, ZDGP, or ZTGP.
- Includes decision-rule visualizations for individual classification.

---

## Citation

If you use this repository, please cite:

> Molina, K.C., Martínez-Minaya, J., Alvares, D., & Tomazella, V.D. (2025).  
> *A Bayesian Survival Model Induced by Hurdle Zero-Modified Power Series Discrete Frailty with Dispersion: An Application in Lung Cancer*. arXiv:2505.23568.

