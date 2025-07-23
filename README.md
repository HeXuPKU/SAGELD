# SAGELD introduction

SPA<sub>GRM</sub> is a genome-wide gene-environmental interaction method of longitudinal traits while controlling for sample relatedness in large-scale datasets. 

**SAGELD is now implemented in the [GRAB package](https://wenjianbi.github.io/grab.github.io/). Please click here to download.**

**Detailed documentation about how to use SAGELD is available at [SAGELD documentation](https://hexupku.github.io/SPAGRM.github.io/docs/Step%201c%20Lonigtudinal%20trand.html).**

**Summary statistics of genome-wide gene-environmental interaction analyses in UK Biobank is available at [here](https://zenodo.org/records/16345203).**


# SAGELD reproducibility

This repository contains: 
1) In-house implementations of several existing longitudinal data analysis methods such as GALLOP and SCEBE;
2) Scripts to reproduce the experiments performed for the SAGELD manuscript, including real data analyses and simulation studies.

## Real data analysis
### preprocess longitudinal traits

In this paper, we extracted 11 longitudinal traits from UKBB primary care data and combined measurements from baseline assessments.
1) To extracted longitudinal traits from gp clinical data, users can refer to the code repository of [here](https://www.nature.com/articles/s41467-025-56669-1).
2) To extracted drug information from gp scripts data, users can refer to the code repository of [here](https://www.nature.com/articles/s41467-025-58152-3).
3) To extracted drug information from baseline assessments, users can refer to the code repository of [here](https://www.nature.com/articles/s41467-019-09572-5).

After preparing these data, we add medication information to each observation by
```
ProcessPhenotypes/1.Process_Drug_Information:
antidiabetics.R           # define treatment intervals and add an binary indicator to decide whether a specific blood glucose measurement is in the treatment interval.
antihypertensives.R       # define treatment intervals and add an binary indicator to decide whether a specific blood pressure measurement is in the treatment interval.
statins.R                 # define treatment intervals and add an binary indicator to decide whether a specific blood lipid measurement is in the treatment interval.
```

We then combined measurements from UKBB primary care data and baseline assessments by
```
ProcessPhenotypes/2.Combine_Primary_Care_and_Baseline_Assessments:
add_baseline_assessments.R # Combine the data from the two sources in chronological order and an binary indicator to distinguish the source of the data.
```

For GxBMI analyses, we match BMI values to four cardiovascular traits by
```
ProcessPhenotypes/3.Match_BMI_to_Cardiovascular_Traits:
add_imputed_BMI.R # Match BMI values to those measured no later than the trait observation date and no earlier than one year prior. 
```


### Gxage analyses

We applied SAGELD to perform Gxage analyses for 11 longitudinal traits extracted from UKBB primary care data or from combining baseline assessments via
```
Gxage_analyses:
SAGELD_XXX_trajectory.R
```

### GxBMI analyses

We applied SAGELD to perform Gxage analyses for four cardiovascular traits via
```
Gxage_analyses:
SAGELD_XXX_bmi.R         # we truncated BMI values within trait-specific ranges prior to analysis to ensure an approximately linear relationship.
```

## Simulation studies

### In-house implementations of several existing longitudinal data analysis methods, including

```
GALLOP_Implementation.R  # In-house implementation of GALLOP algorithm
SCEBE_Implementation.R   # In-house implementation of SCEBE algorithm
```

### 1. genotype simulation

To mimic the genotype distribution in real data, we simulated genotype data using real genotype data of White British subjects in UK Biobank by performing gene-dropping simulations. We simulated 100,000 common variants (MAF > 0.01) from genotype calls (field ID: 22418) and sequencing data (field ID: 23155), respectively. Users can refer to the code repository of [here](https://github.com/HeXuPKU/SPAGRM/tree/main/simulation/1.simulate_genotype) for more details.

```
SimulateGenotype.R       # contains R script to simulate genotypes.
```

### 2. phenotype simulation

We simulated longitudinal traits with varying numbers of observations and structures of between-subject (BS) and within-subject (WS) variabilities.

```
SimulatePhenotype.R      # contains R script to simulate longitudinal traits.
```

### 3. type I error rates

We conducted 1e9 replications in scenario 1 and 1e3 replications in scenario 2 using genotypes and phenotypes simulated above. 

```
typeIerror/typeIerror_quan1.R      # contains R scripts to conduct 1e9 type I error simulations for continuous environmental exposures in scenario 1 (betaG = 0 & betaGxE = 0).
typeIerror/typeIerror_binary1.R    # contains R scripts to conduct 1e9 type I error simulations for discrete environmental exposures in scenario 1 (betaG = 0 & betaGxE = 0).
typeIerror/typeIerror_quan2.R      # contains R scripts to conduct 1e9 type I error simulations for continuous environmental exposures in scenario 2 (betaG ≠ 0 & betaGxE = 0).
typeIerror/typeIerror_binary2.R    # contains R scripts to conduct 1e9 type I error simulations for continuous environmental exposures in scenario 2 (betaG ≠ 0 & betaGxE = 0).
```

### 4. statistical power

We conducted 1e3 replications in scenario 3 and 4 using genotypes and phenotypes simulated above. 

```
typeIerror/typeIerror_quan3.R      # contains R scripts to conduct 1e9 type I error simulations for continuous environmental exposures in scenario 3 (betaG = 0 & betaGxE ≠ 0).
typeIerror/typeIerror_binary3.R    # contains R scripts to conduct 1e9 type I error simulations for discrete environmental exposures in scenario 3 (betaG = 0 & betaGxE ≠ 0).
typeIerror/typeIerror_quan4.R      # contains R scripts to conduct 1e9 type I error simulations for continuous environmental exposures in scenario 4 (betaG ≠ 0 & betaGxE ≠ 0).
typeIerror/typeIerror_binary4.R     # contains R scripts to conduct 1e9 type I error simulations for continuous environmental exposures in scenario 4 (betaG ≠ 0 & betaGxE ≠ 0).
```
