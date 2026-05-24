# DCB PYP Model

This repository contains an R implementation of a copy-number alteration (CNA) timing model for PCAWG data. The main script is `DCB_PYP_Final2.R`.

The project starts from the TickTack idea of estimating when copy-number gains occurred from mutation burdens and variant allele frequencies (VAFs), then extends the temporal clustering step with a Pitman-Yor Process (PYP). The goal is to infer evolutionary "clocks" without fixing the number of clusters in advance, while separating large clonal CNA waves from low-frequency noisy events.

## Main Model

The implemented pipeline has five conceptual stages:

1. Load PCAWG mutation, copy-number, and metadata tables.
2. Estimate segment-specific CNA timing from VAF and copy-number state.
3. Cluster timing estimates using a Pitman-Yor Process mixture.
4. Validate the inferred clusters with retrospective likelihood and posterior predictive checks.
5. Compare the PYP clustering against the original TickTack package.

The model focuses on patient `DO51964` in the saved validation object, but the script automatically searches among candidate samples with many amplifications.

## Timing Estimation

For each amplified copy-number segment, the code collects mutations falling inside the segment and requires at least `MIN_MUT_PER_SEGMENT = 5` mutations.

For a segment with total copy number `N`, tumor purity `pi`, and `j` mutated copies, the expected VAF is:

```text
p_j = (pi * j) / (N * pi + 2 * (1 - pi))
```

The observed mean VAF is converted into the estimated proportion of pre-gain mutations:

```text
theta_2 = (VAF_obs - p_1) / (p_2 - p_1)
theta_1 = 1 - theta_2
```

The segment timing estimate `tau` is then obtained with copy-number-specific formulas:

```text
CNLOH:      tau = (2 * theta_2) / (2 * theta_2 + theta_1)
Trisomy:   tau = (3 * theta_2) / (theta_1 + 2 * theta_2)
Tetrasomy: tau = (2 * theta_2) / (2 * theta_2 + theta_1)
HighAmp:   tau = (theta_2 * total_cn) / (theta_1 + theta_2 * total_cn)
```

This gives each CNA segment a pseudo-time value in `[0, 1]`, where lower values indicate earlier events and higher values indicate later events.

## Pitman-Yor Process Extension

The core extension is in `fit_pitman_yor_clustering()`.

Instead of using a fixed finite mixture, the model uses a truncated stick-breaking Pitman-Yor Process mixture:

```text
G ~ PYP(alpha, d, G0)
G0 = Uniform(0, 1)
tau_i | cluster k ~ Normal(mu_k, sigma^2)
```

The implementation uses:

```text
ALPHA_PYP    = 2.0
DISCOUNT_PYP = 0.5
SIGMA_LIK    = 0.10
K_TRUNC      = 15
```

The discount parameter lets the model produce uneven, power-law-like cluster sizes. This is useful for CNA timing because major clonal bursts may contain hundreds of events, while neutral or subclonal noise may appear as isolated small clusters.

The sampler updates:

- cluster assignments `z_i`
- cluster centers `mu_k`
- stick-breaking weights `v_k`

The Gaussian likelihood is motivated by the second part of the slide deck: segment-level timing estimates are treated as noisy continuous observations, with the Central Limit Theorem supporting an approximate Normal working model when each segment has enough mutations.

## Validation

The script validates the PYP result in two ways.

First, `compute_retrospective_likelihood()` reconstructs the expected mutation mixture for each segment and scores observed VAFs under a binomial model using an approximate depth of `100`.

Second, `posterior_predictive_vaf()` simulates predicted VAF distributions from the inferred timing model and compares observed versus predicted VAFs using Kolmogorov-Smirnov tests per cluster.

These validation steps are designed to check whether the inferred temporal clusters explain the observed VAF structure, not only whether the clusters look separated in timing space.

## TickTack Comparison

The script also converts the inferred segment and mutation data into the format expected by the `tickTack` R package:

```r
convert_to_ticktack()
analyze_patient_ticktack()
compare_methods()
```

The comparison reports:

- number of clusters
- cluster sizes
- within-cluster variance
- inter-cluster separation
- singleton clusters
- TickTack AIC/BIC model selection when available

This makes the PYP model a direct extension of the TickTack framework rather than a completely separate analysis.

## Result Highlight

The saved validation result `DO51964_Validation_Results.rds` reports 7 active PYP clusters for patient `DO51964`:

| Cluster | CNAs | Mean tau | Interpretation |
| --- | ---: | ---: | --- |
| 1 | 106 | 0.044 | Early phase establishment |
| 2 | 193 | 0.372 | Major clonal expansion |
| 3 | 1 | 0.529 | Subclonality / noise |
| 4 | 217 | 0.665 | Sustained progression phase |
| 5 | 1 | 1.000 | Isolated terminal event |
| 6 | 1 | 0.964 | Isolated terminal event |
| 7 | 294 | 0.908 | Late "Hopeful Monster" burst |

The important biological interpretation is that the PYP separates the dominant CNA waves from singleton events. Cluster 7 is the largest late amplification burst, while clusters 3, 5, and 6 are isolated low-frequency events that would otherwise risk contaminating the major temporal clocks.

## Files

```text
DCB_PYP_Final2.R                  Main model and analysis pipeline
DO51964_Validation_Results.rds    Saved validation result for DO51964
DCB_SLIDES_Favagrossa.pdf         Slide deck describing TickTack and the PYP extension
wasserstein_heatmap.png           Supporting result figure
Data/                             PCAWG input data
```

## Requirements

The main script uses R and the following packages:

```r
data.table
dplyr
tickTack
```

The local R environment used during cleanup had `data.table` and `dplyr` installed, but `tickTack` was missing. Install TickTack before running the complete comparison pipeline.

## Running

From the repository root:

```bash
Rscript DCB_PYP_Final2.R
```

The script expects the PCAWG files to be present in `Data/`. Output is saved as:

```text
<patient_id>_CompleteResults.rds
```

## Current Assumptions and Limitations

- Tumor purity is fixed at `PURITY = 1.0`.
- Timing is represented as abstract pseudo-time, not calendar time.
- The PYP sampler is a practical truncated implementation with `K_TRUNC = 15`.
- The retrospective likelihood uses an approximate fixed sequencing depth of `100`.
- Hyperparameters `alpha`, `discount`, and `sigma` may need sensitivity analysis across tumor types.
- Singleton clusters are useful for noise isolation, but broader use would require standardized post-processing rules.

## Thesis Contribution

The central contribution is replacing a fixed temporal clustering approach with a non-parametric Pitman-Yor Process mixture. This allows the model to infer the number and size of evolutionary clocks from the data, capture large CNA bursts, and isolate rare noisy events. In the DO51964 analysis, this produces a high-resolution reconstruction with 7 temporal clusters and a clear late amplification burst.
