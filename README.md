# Bayesian-Mixture-Models

This repository contains R implementations of various Bayesian mixture models, along with experimental scripts for evaluating their performance in different data scenarios. The models are implemented using NIMBLE and JAGS, facilitating flexible Bayesian inference.

## Repository Structure
```bash
.
├── code/
│   ├── Bayes_lddp_var.R        # Covariate-dependent variance LDDP model
│   ├── Bayes_lddp.R            # Standard LDDP model
│   ├── Bayes_mix_normal.R      # Gaussian Mixture Model (GMM)
│   ├── Bayes_mix_poisson.R     # Poisson Mixture Model (PMM)
│   ├── Bayes_mix_student_t.R   # Student-t Mixture Model
│
├── experiments/
│   ├── ex1/          # Experiment 1: Priors of Overfitted Mixtures
│   ├── ex2/          # Experiment 2: Over-dispersed Data
│   ├── ex3/          # Experiment 3: Under-dispersed Data
│   ├── ex4/          # Experiment 4: Heavy-Tailed Data
│   ├── ex_DP/        # Experiment on Dirichlet Process Mixtures
│   ├── ex_covariate/ # Experiment on Covariate-Dependent Mixtures
│
└── README.md  # Project overview
```

## Setup
Before running the models, ensure that NIMBLE and RJAGS are installed on your system.

