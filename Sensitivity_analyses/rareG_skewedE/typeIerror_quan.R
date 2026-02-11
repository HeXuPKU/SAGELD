
args = commandArgs(TRUE)
print(args)
print(sessionInfo())
n.cpu = as.numeric(args)
print(n.cpu)

library(dplyr)
library(tidyr)
library(lme4)
library(optimx)
library(afex)
library(GRAB)
source("your_working_directory/SimulatePhenotype_skewed_Env.R")

Envtype = "quan"

tables = data.frame(nrep = rep(1:100, times = 3),
                    Skewness = rep(c("Normal", "Moderate", "High"), each = 100))

nrep = tables$nrep[n.cpu]; Skewness = tables$Skewness[n.cpu]

nSub = 25000; nFam = 2500
randMat = data.table::fread("your_working_directory/randMat10Members.txt")
bvectomain = c(as.numeric(t(randMat[sample(1e6, nFam),])), rnorm(nSub))
bvectoGxE = c(as.numeric(t(randMat[sample(1e6, nFam),])), rnorm(nSub))

longpheno = pheno_simu(nSub = nSub, nFam = nFam, FamMode = "10-members", 
                       medianObs = 4, Env = "quan", Skewness = Skewness, betaG = 0, betaGxE = 0,  
                       Geno = NULL, bVectomain = bvectomain, bVectoGxE = bvectoGxE, rho = 0.8)

nullmodel = lmer(pheno ~ xone + xtwo + xthree + Env + (Env|SubjID), data = longpheno, 
                 control = lmerControl(optimizer ='optimx', optCtrl = list(method = 'nlminb')))


SAGELD_NULL_Model = SAGELD.NullModel(NullModel = nullmodel,
                                     PlinkFile = "your_working_directory/nSub_25000_nFam_2500",
                                     SparseGRMFile = "your_working_directory/SparseGRM_0.05.txt",
                                     PairwiseIBDFile = "your_working_directory/PairwiseIBD.txt",
                                     control = list(ControlOutlier = FALSE))

GRAB.Marker(objNull = SAGELD_NULL_Model,
            GenoFile = "your_working_directory/nSub_25000_nFam_2500.bed",
            OutputFile = paste0("your_working_directory/common_", Envtype, "_Skewness_", Skewness, "/SAGELD.", nrep, ".txt"),
            control = list(SPA_Cutoff = 2,
                           zeta = 0,
                           tol = 1e-4,
                           impute_method = "mean",
                           missing_cutoff = 0.05,
                           min_maf_marker = 2e-4,
                           min_mac_marker = 20,
                           nMarkersEachChunk = 10000))

GRAB.Marker(objNull = SAGELD_NULL_Model,
            GenoFile = "your_working_directory/nSub_25000_nFam_2500Rare.bed",
            OutputFile = paste0("your_working_directory/rare_", Envtype, "_Skewness_", Skewness, "/SAGELD.", nrep, ".txt"),
            control = list(SPA_Cutoff = 2,
                           zeta = 0,
                           tol = 1e-4,
                           impute_method = "mean",
                           missing_cutoff = 0.05,
                           min_maf_marker = 2e-4,
                           min_mac_marker = 20,
                           nMarkersEachChunk = 10000))
