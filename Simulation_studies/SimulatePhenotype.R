
library(dplyr)
library(tidyr)
#### function to generate longitudinal phenotype with GxE interaction.
pheno_simu = function(nSub,                           # number of independent individuals.
                      nFam = 0,                       # number of families.
                      FamMode = "4-members",          # Structure of families, default is 4.
                      # mi = 6:15,                    # Vector of observations for each individual, out of use.
                      medianObs = 4,                  # median of number of observations for each individual.
                      Env = "quan",                   # Environment is quantitative/binary variable.
                      MedFreq = 0.3,                  # Mimic medication frequency (Env = "bianry").
                      MedPrpo = 0.5,                  # Mimic proportion of people on medication.
                      betaG = 0,                      # Genetic effect for mean level, a number.
                      betaGxE = 0,                    # Genetic effect for GxE interaction, a number.
                      Geno = NULL,                    # Genotype vector, default equal NULL.
                      bVectomain = NULL,              # random effect for main effect, match the family structure.
                      bVectoGxE = NULL,               # random effect for interaction effect, match the family structure.
                      rho = 0.8,                      # weight parameter to weight random intercept and slope effects.
                      hete = "medium")                # heterogeneity parameter to determine between-subject variability and within-subject variability.
{
  if(nFam == 0){
    N = nSub
    SubjID = UnrelatSubjID = paste0("U-", 1:nSub)
  }else{
    if(FamMode == "4-members"){
      N = nSub + 4 * nFam
      if(nSub == 0)
      {
        SubjID = paste0("F4-", rep(1:nFam, each=4), "-", 1:4)
        UnrelatSubjID = paste0("F4-", rep(1:nFam, each=2), "-", 1:2)
      }else
      {
        SubjID = c(paste0("F4-", rep(1:nFam, each=4), "-", 1:4), paste0("U-", 1:nSub))
        UnrelatSubjID = c(paste0("F4-", rep(1:nFam, each=2), "-", 1:2), paste0("U-", 1:nSub))
      }
    }else if(FamMode == "10-members"){
      N = nSub + 10 * nFam
      if(nSub == 0)
      {
        SubjID = paste0("F10-", rep(1:nFam, each=10), "-", 1:10)
        UnrelatSubjID = paste0("F10-", rep(1:nFam, each=4), "-", 1:4)
      }else
      {
        SubjID = c(paste0("F10-", rep(1:nFam, each=10), "-", 1:10), paste0("U-", 1:nSub))
        UnrelatSubjID = c(paste0("F10-", rep(1:nFam, each=4), "-", 1:4), paste0("U-", 1:nSub))
      }
    }else if(FamMode == "20-members"){
      N = nSub + 20 * nFam
      if(nSub == 0)
      {
        SubjID = paste0("F20-", rep(1:nFam, each=20), "-", 1:20)
        UnrelatSubjID = paste0("F20-", rep(1:nFam, each=8), "-", 1:8)
      }else
      {
        SubjID = c(paste0("F20-", rep(1:nFam, each=20), "-", 1:20), paste0("U-", 1:nSub))
        UnrelatSubjID = c(paste0("F20-", rep(1:nFam, each=8), "-", 1:8), paste0("U-", 1:nSub))
      }
    }else{
      stop("Please check if the FamMode is in '4-members', '10-members' or '20-members'!")
    }
  }
  
  if(is.null(Geno)){
    Geno = rep(0, N)
  }else if(length(Geno) != N){
    stop("Length of geno should correspond to number of individuals.")
  }
  
  if(is.null(bVectomain)){
    bVectomain = rep(0, N)
  }else if(length(bVectomain) != N){
    stop("Length of random vector should correspond to number of individuals.")
  }
  if(is.null(bVectoGxE)){
    bVectoGxE = rep(0, N)
  }else if(length(bVectoGxE) != N){
    stop("Length of random vector should correspond to number of individuals.")
  }
  
  # mi = 6:15 is out of use, now I use gamma distribution gamma(1.5, 0.3) to mimic LDL in real scenario.
  # m = sample(mi, N, replace = TRUE); M = sum(m) # number of measurements per subject and total measurements.
  
  if(medianObs == 2)
  {
    m = ceiling(5 * rgamma(n = N, shape = 1.5, scale = 0.3)); M = sum(m)
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    # 1.000   1.000   2.000   2.764   4.000  18.000 
  }else if(medianObs == 4)
  {
    m = ceiling(10 * rgamma(n = N, shape = 1.5, scale = 0.3)); M = sum(m)
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    # 1.000   2.000   4.000   4.975   7.000  35.000 
  }else if(medianObs == 8)
  {
    m = ceiling(20 * rgamma(n = N, shape = 1.5, scale = 0.3)); M = sum(m)
    # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    # 1.000   4.000   8.000   9.518  13.000  74.000
  }else
  {
    stop("Current this function don't support this option.")
  }
  
  # generate covariates and corresponding coefficients.
  xone = rbinom(N, 1, 0.5); xone = rep(xone, times = m) # time-invariant binary variable, mimic sex.
  xtwo = rnorm(N, 0, 1); xtwo = rep(xtwo, times = m)    # time-invariant standard normal variable, mimic height.
  xthree  = rnorm(M, 0, 1)                              # time-varying standard normal variable, mimic weight.
  Geno = rep(Geno, times = m)                           # hard-called genotypes.
  
  if(Env == "quan")
  {
    Env = rnorm(M, 0, 1)    # time-varying standard normal variable, mimic age.
  }else if(Env == "binary")
  {
    # Env = rbinom(M, 1, MedFreq)  # time-varying binary variable, mimic medication.
    
    isonMed = sample(c(TRUE, FALSE), N, replace = TRUE, prob = c(MedPrpo, 1 - MedPrpo))
    
    Env = lapply(1:N, function(i){
      if(isonMed[i]){
        rbinom(m[i], 1, MedFreq)
      }else{
        rep(0, m[i])
      }
    }) %>% unlist(use.names = FALSE)
    
  }else
  {
    stop("Env should be 'quan' or 'binary'.\n")
  }
  
  # corresponding coefficients.
  beta = c(1, 0.5, 0.5, 0.5, 0.5)
  
  # generate random effects.
  randomeffect_main = rep(rnorm(N), times = m);
  randomeffect_main = sqrt(rho) * randomeffect_main
  
  randomeffect_GxE = rep(rnorm(N), times = m) * Env;
  randomeffect_GxE = sqrt(1 - rho) * randomeffect_GxE
  
  # generate polygenic effects.
  polygeniceffect_main = rep(bVectomain, times = m);
  polygeniceffect_main = sqrt(rho) * polygeniceffect_main
  
  polygeniceffect_GxE = rep(bVectoGxE, times = m) * Env;
  polygeniceffect_GxE = sqrt(1 - rho) * polygeniceffect_GxE
  
  error_term = rnorm(M, 0, 1) # Within subject variability

  # long_pheno = beta[1] + beta[2] * xone + beta[3] * xtwo + beta[4] * xthree + beta[5] * Env +
  #   betaG * Geno + betaGxE * Geno * Env + randomeffect_main + randomeffect_GxE + polygeniceffect_main + polygeniceffect_GxE + error_term
  
  if(hete == "medium")
  {
    # between-subject variability : within-subject variability = 1.6:1
    long_pheno = beta[1] + beta[2] * xone + beta[3] * xtwo + beta[4] * xthree + beta[5] * Env +
      betaG * Geno + betaGxE * Geno * Env + randomeffect_main + randomeffect_GxE + polygeniceffect_main + polygeniceffect_GxE + error_term
    
  }else if(hete == "high")
  {
    # between-subject variability : within-subject variability = 2.2:0.4
    long_pheno = beta[1] + beta[2] * xone + beta[3] * xtwo + beta[4] * xthree + beta[5] * Env +
      betaG * Geno + betaGxE * Geno * Env + sqrt(1.75) * randomeffect_main + randomeffect_GxE + polygeniceffect_main + polygeniceffect_GxE + sqrt(0.4) * error_term
    
  }else if(hete == "low")
  {
    # between-subject variability : within-subject variability = 1:1.6
    long_pheno = beta[1] + beta[2] * xone + beta[3] * xtwo + beta[4] * xthree + beta[5] * Env +
      betaG * Geno + betaGxE * Geno * Env + sqrt(0.25) * randomeffect_main + randomeffect_GxE + polygeniceffect_main + polygeniceffect_GxE + sqrt(1.6) * error_term
    
  }else
  {
    stop("Currently, we only support hete = high, medium, or low.\n")
  }
  
  covdata = data.frame(xone = xone, xtwo = xtwo, xthree = xthree, Env = Env, pheno = long_pheno) %>%
    mutate(SubjID = rep(SubjID, times = m)) %>% mutate(UNRELATED = SubjID %in% UnrelatSubjID)
  
  return(covdata)
}
