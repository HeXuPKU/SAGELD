
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
source("your_working_directory/SimulatePhenotype_time_invariant.R")

scr = 2; Envtype = "quan"

tables = data.frame(nrep = 1:100)

nrep = tables$nrep[n.cpu]

GenoMat = GRAB::GRAB.ReadGeno(GenoFile = "your_working_directory/causalSNPs.bed",
                              control = list(AllMarkers = TRUE,
                                             imputeMethod = "mean"))

maf = colMeans(GenoMat$GenoMat)/2
GenoMatSD=t((t(GenoMat$GenoMat) - 2*maf)/sqrt(2*maf*(1-maf)))

Geno = GenoMatSD %*% rep(0.1, 10) %>% as.numeric()

nSub = 25000; nFam = 2500
randMat = data.table::fread("your_working_directory/randMat10Members.txt")
bvectomain = c(as.numeric(t(randMat[sample(1e6, nFam),])), rnorm(nSub))
bvectoGxE = c(as.numeric(t(randMat[sample(1e6, nFam),])), rnorm(nSub))

longpheno = pheno_simu(nSub = nSub, nFam = nFam, FamMode = "10-members", 
                       medianObs = 4, Env = "quan", betaG = 1, betaGxE = 0,  
                       Geno = Geno, bVectomain = bvectomain, bVectoGxE = bvectoGxE, rho = 0.8)

unrelated_subject = c(paste0("F10-", rep(1:2500, each = 4), "-", 1:4), paste0("U-", 1:25000))

nullmodel_full = lmer(pheno ~ xone + xtwo + xthree + Env + (1|SubjID), data = longpheno,
                      control = lmerControl(optimizer ='optimx', optCtrl = list(method = 'nlminb')))

# or alternative
# nullmodel_full = lmer(pheno ~ xone + xtwo + xthree + Env + (Env|SubjID), data = longpheno,
#                       control = lmerControl(optimizer ='optimx', optCtrl = list(method = 'nlminb')))


SAGELD_NULL_Model = SAGELD.NullModel(NullModel = nullmodel_full,
                                     PlinkFile = "your_working_directory/nSub_25000_nFam_2500",
                                     SparseGRMFile = "your_working_directory/SparseGRM_0.05.txt",
                                     PairwiseIBDFile = "your_working_directory/PairwiseIBD.txt",
                                     control = list(ControlOutlier = FALSE))

GRAB.Marker(objNull = SAGELD_NULL_Model,
            GenoFile = "your_working_directory/causalSNPs.bed",
            OutputFile = paste0("your_working_directory/", Envtype, "/SAGELD.", nrep, ".txt"),
            control = list(SPA_Cutoff = 2,
                           zeta = 0,
                           tol = 1e-4,
                           impute_method = "mean",
                           missing_cutoff = 0.05,
                           min_maf_marker = 2e-4,
                           min_mac_marker = 20,
                           nMarkersEachChunk = 10000))
