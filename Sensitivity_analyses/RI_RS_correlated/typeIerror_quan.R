
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
library(MASS)
source("your_working_directory/SimulatePhenotype_RI_RS_correlated.R")

scr = 1; Envtype = "quan"

tables = data.frame(nrep = rep(1:100, times = 3),
                    rho = rep(c(0, 25, 50), each = 100))

nrep = tables$nrep[n.cpu]; rho = tables$rho[n.cpu]

nSub = 25000; nFam = 2500
randMat = data.table::fread("your_working_directory/randMat10Members.txt")
bvectomain = c(as.numeric(t(randMat[sample(1e6, nFam),])), rnorm(nSub))
bvectoGxE = c(as.numeric(t(randMat[sample(1e6, nFam),])), rnorm(nSub))

N = nSub + 10 * nFam
randomeffect = mvrnorm(N, c(0, 0), Sigma = matrix(c(1, rho/100, rho/100, 1), 2, 2))
randomemain = randomeffect[,1]
randomGxE = randomeffect[,2]

GenoMat = GRAB::GRAB.ReadGeno(GenoFile = "your_working_directory/causalSNPs.bed",
                              control = list(AllMarkers = TRUE,
                                             imputeMethod = "mean"))

maf = colMeans(GenoMat$GenoMat)/2
GenoMatSD=t((t(GenoMat$GenoMat) - 2*maf)/sqrt(2*maf*(1-maf)))

Geno = GenoMatSD %*% rep(0.1, 10) %>% as.numeric()

longpheno = pheno_simu(nSub = nSub, nFam = nFam, FamMode = "10-members", 
                       medianObs = 4, Env = "quan", betaG = 1, betaGxE = 0,  
                       Geno = Geno, randomemain = randomemain, randomGxE = randomGxE,
                       bVectomain = bvectomain, bVectoGxE = bvectoGxE, rho = 0.8)

nullmodel = lmer(pheno ~ xone + xtwo + xthree + Env + (Env|SubjID), data = longpheno, 
                 control = lmerControl(optimizer ='optimx', optCtrl = list(method = 'nlminb')))


SAGELD_NULL_Model = SAGELD.NullModel(NullModel = nullmodel,
                                     PlinkFile = "your_working_directory/nSub_25000_nFam_2500",
                                     SparseGRMFile = "your_working_directory/SparseGRM_0.05.txt",
                                     PairwiseIBDFile = "your_working_directory/PairwiseIBD.txt",
                                     control = list(ControlOutlier = FALSE))

GRAB.Marker(objNull = SAGELD_NULL_Model,
            GenoFile = "your_working_directory/causalSNPs.bed",
            OutputFile = paste0("your_working_directory/", Envtype, rho, "/SAGELD.", nrep, ".txt"),
            control = list(SPA_Cutoff = 2,
                           zeta = 0,
                           tol = 1e-4,
                           impute_method = "mean",
                           missing_cutoff = 0.05,
                           min_maf_marker = 2e-4,
                           min_mac_marker = 20,
                           nMarkersEachChunk = 10000))
