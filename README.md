# Replication: Synthetic Confounder Sensitivity Analysis

## Overview

This repository contains the code for replicating the empirical
application in Section 5 and 6 of the paper. The application estimates
the effect of age-10 cognitive skills on labour-market outcomes at
age 46 using data from the **1970 British Cohort Study (BCS70)**,
and applies the synthetic confounder sensitivity analysis to assess
robustness to omitted variable bias.

## Data access

The BCS70 is managed by the **Centre for Longitudinal Studies (CLS)**
at UCL and distributed through the **UK Data Service (UKDS)**. The
data are freely available for academic use but require registration.

1. Register at [https://ukdataservice.ac.uk](https://ukdataservice.ac.uk).
2. Download the following BCS70 studies. Each study corresponds to a
   survey wave or set of linked records:

| Study Number | Title / Content | Wave / Age |
|:-------------|:----------------|:-----------|
| SN 2666 | BCS70: Birth and 22-Month Subsample | Birth (1970) |
| SN 2699 | BCS70: Five-Year Follow-Up | Age 5 (1975) |
| SN 3723 | BCS70: Ten-Year Follow-Up | Age 10 (1980) |
| SN 8547 | BCS70: Age 46 Survey (2016) | Age 46 (2016) |

Additional linked files may be required depending on which
covariate blocks are used. The control variables span seven
covariate sets (set1–set7) covering birth and perinatal
characteristics, age-5 family background, age-5 and age-10 health,
parental education, maternal employment, accommodation, parental
job satisfaction, medical and sensory indicators, physical
development, and teacher-rated classroom behaviour.

3. Place the downloaded Stata files in a directory of your choice
   (referred to below as `$BCSRAW`).

## Variables used

### Treatment

| Variable | Description | Source |
|:---------|:------------|:-------|
| `pc1` | First principal component of five age-10 cognitive test scores (Fun Maths Test; BAS word definitions, number recall, similarities, matrices), each age-standardised before PCA | Constructed in `datapreparation_2.do` from SN 3723 |

### Outcomes (age 46, Study SN 8547)

| Variable | Description |
|:---------|:------------|
| `employedFT` | Binary: employed full-time at age 46 |
| `isManagerSupervisor` | Binary: holds a managerial or supervisory role at age 46 |

Both outcomes are drawn from the 2016 wave only (age 46). The
underlying panel file `panel_labour_outcomes.dta` contains multiple
waves per individual; the replication code restricts to `wave == 2016`
to avoid treating repeated observations as independent.

### Controls

The full control set is:

```
global cogcontrols perserverance_2_sd $set1 $set2 $set3 $set4 $set5 $set6 $set7
```

where the covariate blocks are:

| Block | Content | Source wave |
|:------|:--------|:------------|
| `perserverance_2_sd` | Teacher-rated perseverance (item j139), standardised | Age 10 |
| `set1` | Mother's age at delivery (linear + quadratic), region of birth (10 dummies), maternal smoking during pregnancy, abnormal gestation, sex, birth weight | Birth |
| `set2` | Mother/father stayed in education, parity, previous stillbirth/preterm/miscarriage, breastfed | Age 5 |
| `set3` | Inpatient stay, hearing/sight problems, child health rating (age 5); parental education dummies (mother/father: non-level, secondary, degree); maternal employment (age 5); non-white; no siblings | Age 5 |
| `set4` | Child health (age 10), inpatient stay (age 10), father/mother unemployment/stay-at-home (age 10), benefit receipt, income band dummies (7 bands) | Age 10 |
| `set5` | Accommodation type dummies (7 types), number of rooms | Age 10 |
| `set6` | Mother's/father's job satisfaction, eyesight/hearing/speech problems | Age 10 |
| `set7` | Perseverance (age 10, teacher questionnaire), classroom behaviour principal components (pcB, pcD, pcE, pcF, pcG, pcH) | Age 10 |

Note: `set1` is redefined in the replication script to expand
`i.c_region_birth` and `i.c_sexbirth` into explicit dummy variables,
because the sensitivity package's internal FWL residualisation
cannot consume Stata's `i.` factor-variable syntax.

### Benchmark variables

The sensitivity analysis benchmarks against observed covariates to
calibrate the plausible strength of unobserved confounders:

| Variable | Description |
|:---------|:------------|
| `sex_2` | Female indicator |
| `c_motherStayhome10` | Mother stay-at-home (age 10) |
| `c_fatherUnemployed10` | Father unemployed (age 10) |
| `pcB` | Classroom behaviour principal component B |
| `c_heightAt10_sd` | Height at age 10 (standardised) |
| `c_weightAt10_sd` | Weight at age 10 (standardised) |

## Data preparation

The analysis dataset is built by three preparation scripts, run in
the following order:

```stata
do "$CODE/datapreparation_3.do"
do "$CODE/datapreparation_1.do"
do "$CODE/datapreparation_2.do"
```

These scripts:

1. **`datapreparation_3.do`**: Loads and cleans the birth and
   perinatal survey data. Creates the `isInPerinatalWave` indicator
   and the birth-related control variables.

2. **`datapreparation_1.do`**: Loads and merges the age-5 and age-10
   follow-up data. Creates the family background, health, parental
   education, accommodation, and teacher questionnaire variables.
   Defines the covariate-block globals `$set2`–`$set7`.

3. **`datapreparation_2.do`**: Constructs the cognitive skills
   treatment variable `pc1` via PCA on the five age-standardised
   age-10 test scores. Constructs the perseverance measures.

After running these scripts, the sample is restricted:

```stata
keep if isInPerinatalWave == 1 & pc1 != . ///
     & (perserverance_1_sd != . | perserverance_2_sd != .)
```

The age-46 labour-market outcomes are then merged:

```stata
preserve
    use "$BCSRAW/panel_labour_outcomes.dta", clear
    keep if wave == 2016
    tempfile labour_age46
    save `labour_age46'
restore
merge 1:1 bcsid using `labour_age46'
drop if _merge == 2
drop _merge
```

The final analysis sample contains approximately **5,344**
observations for `employedFT` and **3,615** for
`isManagerSupervisor` (the difference reflects item-level
missingness in the outcome variables).

## Running the sensitivity analysis

The main replication script is `sensitivity_v18_validate_full_with_sims-2.do`.
It contains both the sensitivity package code (programs `sens_abc`,
`sens_build_z`, `sens_fit`, `sens_bisect`, `sensitivity`, etc.) and
the application code at the bottom.

Before running, edit the file paths at the top of the application
section to point to your local directories:

```stata
* Edit these paths:
do "$CODE/datapreparation_3.do"         // -> your path to datapreparation scripts
do "$CODE/datapreparation_1.do"
do "$CODE/datapreparation_2.do"
use "$BCSRAW/panel_labour_outcomes.dta"  // -> your path to BCS70 data
log using "$OUTPUT/bcs70_sens.smcl"      // -> your desired output path
```

The script runs:

1. **Sensitivity grid**: For each outcome, the sensitivity command
   evaluates the adjusted AME on an 11×11 grid of
   (ρ_D, ρ_Y) ∈ {−0.5, −0.4, ..., 0.5}².

2. **Bisection**: Computes explain-away (ρ_D\*) and sign-change
   (ρ_D⁺) thresholds at each ρ_Y level.

3. **Benchmarking**: Evaluates the sensitivity of the AME when the
   confounding strength matches that of each benchmark variable.

4. **Profile plots**: Generates sensitivity profiles at ρ_Y = 0.1
   and ρ_Y = 0.3 for each outcome.

Outputs are saved to `bcs70_sens.smcl`.

## Running the simulations

The simulation studies (Section [X] of the paper) are also contained
in `sensitivity_v18_validate_full_with_sims-2.do`. They do not
require BCS70 data — all data are generated within the script. To
run the simulations only:

```stata
set seed 12345
sim_v3, nsim(10000) nobs(5000)
```

Output is saved to `simulations.smcl`. The simulation generates
10,000 replications of 5,000 observations each across six DGP cells,
comparing the synthetic confounder method against CH-OLS and
CH-adapted.

## Software requirements

- **Stata 16** or later (the script sets `version 16`)
- No additional Stata packages are required; the sensitivity analysis
  programs are self-contained within the .do file

## Important caveats

1. The analysis uses **unweighted logit** without clustering. The
   BCS70 application in the paper is an illustration of the
   sensitivity method, not a full analysis of the cohort data.
   Production-quality estimates would incorporate attrition weights
   and clustered standard errors.

2. Only **binary outcomes** are analysed. The sensitivity package
   supports logit, probit, cloglog, poisson, nbreg, ologit, and
   mlogit.  

3. The `sensitivity_profile` command (called for the profile plots)
   is used only for the plots in the appendix (mlogit). You may
   comment out the four `sensitivity_profile` lines before the
   `foreach` loop; the main grid, bisection, and benchmark results
   are unaffected.

## File manifest

| File | Description |
|:-----|:------------|
| `sensitivity_v18_validate_full_with_sims-2.do` | Main replication script (sensitivity package + simulations + BCS70 application) |
| `datapreparation_1.do` | Data preparation: age-5/age-10 variables and covariate blocks |
| `datapreparation_2.do` | Data preparation: PCA for cognitive skills, perseverance |
| `datapreparation_3.do` | Data preparation: birth/perinatal variables |
| `simulations.smcl` | Simulation output log |
| `bcs70_sens.smcl` | BCS70 application output log |
| `README.md` | This file |

## Citation

If you use this code, please cite:

> Fé, E. (2026). "Synthetic Confounder Sensitivity Analysis for
> Nonlinear Models." Available at SSRN: [Working paper 7395379](http://dx.doi.org/10.2139/ssrn.7395379)
