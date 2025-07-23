# SPA<sub>GRM</sub> introduction

SPA<sub>GRM</sub> is a scalable and accurate analysis framework to control for sample relatedness in large-scale genome-wide association studies (GWAS). In the paper [SPA<sub>GRM</sub>: effectively controlling for sample relatedness in large-scale genome-wide association studies of longitudinal traits](https://www.nature.com/articles/s41467-025-56669-1), we applied SPA<sub>GRM</sub> to analyze 79 longitudinal traits extracted from UK Biobank primary care data.

**SPA<sub>GRM</sub> is now implemented in the [GRAB package](https://wenjianbi.github.io/grab.github.io/). Please click here to download.**

**Detailed documentation about how to use SPA<sub>GRM</sub> is available at [SPA<sub>GRM</sub> documentation](https://hexupku.github.io/SPAGRM.github.io/).**

**Summary statistics of 79 longitudinal traits extracted from UK Biobank primary care data is available at [here](https://zenodo.org/records/14633793).**

# About this repository

This repository contains: 
1) old version of SPA<sub>GRM</sub>;
2) manhattan and QQ plots of GWAS results for 79 longitudinal traits;
3) materials for reproducing the experiments, including real data analyses and simulation studies.

# SPA<sub>GRM</sub> reproducibility

Scripts to reproduce the experiments performed for the SPA<sub>GRM</sub> manuscript:

**A scalable and accurate analysis framework to control for sample relatedness in large-scale genome-wide association studies and its application to 79 longitudinal traits in UK Biobank (to be updated)**

## Real data analysis
### 1. extract longitudinal traits from UK Biobank primary care data

In this paper, we extracted 79 longitudinal traits from UKB primary care data. The UK Biobank (UKB) primary care data (Category ID: 3001) is derived from electronic health records (EHRs) maintained by General Practitioners (GPs) from multiple data providers in England, Wales, and Scotland. As of the latest release in September 2019, approximating 230,000 UKB participants have been linked to their corresponding primary care data. This dataset includes clinical event records (Field ID: 42040) spanning over 30 years, rich in information of diagnoses, history, symptoms, lab results, and procedures. Two controlled clinical terminologies, Read version 2 (Read v2) and Clinical Terms Version 3 (CTV3) are used to record these primary clinical events. Read v2 and CTV3 clinical terms to extract 79 longitudinal traits are available at real_data\1.extract_pheno\phecode-2023-05-10XH.R, and gp_clinical data should be downloaded from UKB primary care data (Field ID: 42040).

```
In real_data\1.extract_pheno:
phecode-2023-05-10XH.R           # defines the Read v2 and CTV3 code terms for 79 longitudinal traits.
real_data_extra-2023-05-01XH.R   # preprocesses the gp_clinical table.
```

### 2. preprocess longitudinal traits

We used the following code to extract longitudinal traits and perform quality control (QC). QC is very important for longitudinal trait analyses. We excluded records containing implausible values and unit errors. We also filtered out outliers at the population and individual level.

```
In real_data\2.preprocess:
XXX_preprocess.R   # extracts longitudinal traits from gp_clinical table and preprocesses each longitudinal trait.
```

### 3. fit the null model via WiSER package

We use the WiSER, a julia package to fit the linear mixed model for each longitudinal trait. User can also fit other models like generalized estimation equations (GEE) as diaplayed in our [online tutorials](https://hexupku.github.io/SPAGRM.github.io/).

```
In real_data\3.null_fitter:
XXX_null_fitter.jl   # fits the null model and obtains model residuals for each longitudinal trait.
```

### 4. GWAS analysis
- SPA<sub>GRM</sub> analysis 
```
In real_data\4.1.SPAGRM:
XXX_SPAGRM.R   # conducts the SPAGRM analysis with two parts: 1). pre-calculate of the joint distribution of genotypes and 2). conduct genome-wide association studies for each SNP.
```
- TrajGWAS analysis
```
In real_data\4.2.TrajGWAS:
XXX_TrajGWAS.jl   # conducts the TrajGWAS analysis.
```

### 5. PRS adjustment (optional)

SPA<sub>GRM</sub> can further gain statistical power through incorporating polygenic scores (PGSs) as covariates with fixed effects. It's an optional term. We employ the idea to implement a two-stage strategy, SPA<sub>GRM</sub>-PGS. In stage 1, we conduct the first round of GWAS via SPA<sub>GRM</sub> and then calculate the Leave One Chromosome Out (LOCO)-PGS based on the summary statistics. In stage 2, the LOCO-PGS is included as an additional covariate for a second round of GWAS via SPA<sub>GRM</sub>. In supplementary note, we used three longitudinal traits to demonstrate this.

```
In real_data\5.PRS_adjustment (optional):
1.select_index_SNPs_from_results_of_SPAGRM.R   # selects independent SNPs that pass the given p value threshold based on the GWAS results of SPAGRM.
2.estimate_effect_sizes_XXX.jl                 # estimate the effect sizes of index SNPs via WiSER package.
3.construct_LOCO-PRS.R                         # constructs the Leave One Chromosome Out Polygenic Risk Scores (LOCO-PRS).
4.adjust_PRS_in_null_model_fitting_XXX.jl      # adds the PRS into the null model fitting.
5.second_round_of_GWAS_via_SPAGRM.R            # conducts the second round of GWAS via SPAGRM.
```

## Simulation study

### 1. genotype simulation

To mimic the genotype distribution in real data, we simulated genotype data using real genotype data of White British subjects in UK Biobank by performing gene-dropping simulations. We simulated 100,000 common variants (MAF > 0.05) and rare variants (MAF < 0.05 and MAC > 20) from genotype calls (field ID: 22418) and sequencing data (field ID: 23155), respectively. The following R script may be lengthy, but the general idea is using the real genotype file as the input of function `GRAB.SimuGMatFromGenoFile()` in GRAB package.

```
In simulation\1.simulate_genotype:
SimulatedGenotypeUsingWES-2023-02-16XH.R   # contains R script to simulate genotypes.
```

### 2. phenotype simulation

We simulated longitudinal traits following TrajGWAS model. The number of measurements was simulated equally distributed ranging from 6 to 15. Three covariates in the mean and within-subject (WS) variability formulas were simulated as: the first one is time-invariant following a Bernoulli distribution with a probability of 0.5; the second one is time-invariant variable following the standard normal distribution, and the third one is time-varying with each measurement following an independent standard normal distribution. One time-varying covariate was simulated following an independent standard normal distribution as random slope. We additionally added random effects that followed a multivariate normal distribution related to GRM to mimic sample relatedness. We chose similar parameters as in TrajGWAS. In addition to the above, we also used the inverse-gamma distribution to generate WS variability instead of the log-normal model as in Trajgwas. We also used the GEE model with various working correlation structures to simulate longitudinal traits.

```
In simulation\2.simulate_phenotype:
longitudinal_pheno_simu_function-2023-02-16XH.R   # contains R script to simulate longitudinal traits as described above.
longitudinal_pheno_simu_function-2024-04-16XH.R   # contains R script to simulate longitudinal traits using inverse-gamma distribution to generate WS variability.
longitudinal_pheno_simu_function-2024-08-13XH.R   # contains R script to simulate longitudinal traits using GEE model with various working correlation structures.
```

### 3. empirical type I error rates

In each scenario, we conducted 1e9 type I error simulations using genotypes and phenotypes simulated above. 

```
In simulation\3.typeIerror, contains R scripts to conduct 1e9 type I error simulations.
```

### 4. empirical power

In each scenario, we conducted 1e3 power simulations using genotypes and phenotypes simulated above. 

```
In simulation\4.power, contains R scripts to conduct 1e3 empirical power simulations.
```
