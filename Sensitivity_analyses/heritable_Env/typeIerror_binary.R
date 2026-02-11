
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
source("your_working_directory/SimulatePhenotype.R")

scr = 2; Envtype = "binary"

tables = data.frame(nrep = rep(1:100, times = 3),
                    causalnum = rep(c(500, 1000, 2000), each = 100))

nrep = tables$nrep[n.cpu]; causalnum = tables$causalnum[n.cpu]

GenoMat = GRAB::GRAB.ReadGeno(GenoFile = paste0("your_working_directory/random", causalnum, "SNPs.bed"),
                              control = list(AllMarkers = TRUE,
                                             imputeMethod = "mean"))

maf = colMeans(GenoMat$GenoMat)/2
GenoMatSD=t((t(GenoMat$GenoMat) - 2*maf)/sqrt(2*maf*(1-maf)))

# We are conducting type I error simulation (betaG != 0; betaG2E != 0; betaGxE = 0)
polygenicmain = GenoMatSD %*% rnorm(causalnum, sd = sqrt(1/causalnum)) %>% as.numeric()
polygenicG2E = GenoMatSD %*% rnorm(causalnum, sd = sqrt(1/causalnum)) %>% as.numeric()

nSub = 25000; nFam = 2500
randMat = data.table::fread("your_working_directory/randMat10Members.txt")
bvectoGxE = c(as.numeric(t(randMat[sample(1e6, nFam),])), rnorm(nSub))

longpheno = pheno_simu(nSub = nSub, nFam = nFam, FamMode = "10-members", 
                       medianObs = 4, Env = "binary", MedFreq = 0.7, MedPrpo = 0.3, betaE = 1, betaG = 1, betaGxE = 0,  
                       Geno = NULL, bVectoE = polygenicG2E, bVectomain = polygenicmain, bVectoGxE = bvectoGxE, rho = 0.8)

unrelated_subject = c(paste0("F10-", rep(1:2500, each = 4), "-", 1:4), paste0("U-", 1:25000))

nullmodel_full = lmer(pheno ~ xone + xtwo + xthree + Env + (Env|SubjID), data = longpheno,
                      control = lmerControl(optimizer ='optimx', optCtrl = list(method = 'nlminb')))

nullmodel_unrelat = lmer(pheno ~ xone + xtwo + xthree + Env + (Env|SubjID),
                         data = longpheno %>% filter(SubjID %in% unrelated_subject),
                         control = lmerControl(optimizer ='optimx', optCtrl = list(method = 'nlminb')))

rm(GenoMat, GenoMatSD)
gc()


SPAGRMGE_NULL_Model = SPAGRMGE.NullModel(NullModel = nullmodel_full,
                                         PhenoFile = longpheno,
                                         SubjIDColname = "SubjID",
                                         PhenoColname = "pheno",
                                         CovaColname = c("xone", "xtwo", "xthree"),
                                         Envcolname = "Env",
                                         PlinkFile = "your_working_directory/nSub_25000_nFam_2500",
                                         SparseGRMFile = "your_working_directory/SparseGRM_0.05.txt",
                                         PairwiseIBDFile = "your_working_directory/PairwiseIBD.txt",
                                         control = list(ControlOutlier = FALSE))


GRAB.Marker(objNull = SPAGRMGE_NULL_Model,
            GenoFile = paste0("your_working_directory/random", causalnum, "SNPs.bed"),
            OutputFile = paste0("your_working_directory/EwithG_TIE_", Envtype, causalnum, "/SPAGRMGE.", nrep, ".txt"),
            control = list(SPA_Cutoff = 2,
                           zeta = 0,
                           tol = 1e-4,
                           impute_method = "mean",
                           missing_cutoff = 0.05,
                           min_maf_marker = 2e-4,
                           min_mac_marker = 20,
                           nMarkersEachChunk = 10000))

gc()

SAGELD_NULL_Model = SAGELD.NullModel(NullModel = nullmodel_full,
                                     PlinkFile = "your_working_directory/nSub_25000_nFam_2500",
                                     SparseGRMFile = "your_working_directory/SparseGRM_0.05.txt",
                                     PairwiseIBDFile = "your_working_directory/PairwiseIBD.txt",
                                     control = list(ControlOutlier = FALSE))

GRAB.Marker(objNull = SAGELD_NULL_Model,
            GenoFile = paste0("your_working_directory/random", causalnum, "SNPs.bed"),
            OutputFile = paste0("your_working_directory/EwithG_TIE_", Envtype, causalnum, "/SAGELD.", nrep, ".txt"),
            control = list(SPA_Cutoff = 2,
                           zeta = 0,
                           tol = 1e-4,
                           impute_method = "mean",
                           missing_cutoff = 0.05,
                           min_maf_marker = 2e-4,
                           min_mac_marker = 20,
                           nMarkersEachChunk = 10000))

