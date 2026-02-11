
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
source("your_working_directory/GALLOP.R")
source("your_working_directory/SimulatePhenotype_high_relatedness.R")

scr = 1; Envtype = "quan"

tables = data.frame(nrep = rep(1:100, times = 3),
                    heritability = rep(c("low", "moderate", "high"), each = 100))

nrep = tables$nrep[n.cpu]; heritability = tables$heritability[n.cpu]

nFam = 2500
nsibling = 2500
n_monozygotic_twins = 2500

# simulate polygenic effects
randMat = data.table::fread("your_working_directory/randMat4Members.txt")
bvectomain = c(as.numeric(t(randMat[sample(25e5, nFam),])),
               as.numeric(t(randMat[sample(25e5, nsibling), c("V3", "V4")])),
               rep(rnorm(n_monozygotic_twins), each = 2))

bvectoGxE = c(as.numeric(t(randMat[sample(25e5, nFam),])),
              as.numeric(t(randMat[sample(25e5, nsibling), c("V3", "V4")])),
              rep(rnorm(n_monozygotic_twins), each = 2))

longpheno = pheno_simu(nFam = nFam, nsibling = nsibling, n_monozygotic_twins = n_monozygotic_twins,
                       medianObs = 4, Env = "quan", betaG = 0, betaGxE = 0, Geno = NULL, 
                       bVectomain = bvectomain, bVectoGxE = bvectoGxE, heritability = heritability, rho = 0.8)

nullmodel = lmer(pheno ~ xone + xtwo + xthree + Env + (Env|SubjID), data = longpheno, 
                 control = lmerControl(optimizer ='optimx', optCtrl = list(method = 'nlminb')))


gallop(nullmodel = nullmodel,
       Phenofile = longpheno,
       IDcolname = "SubjID",
       Phenocolname = "pheno",
       Covacolname = c("xone", "xtwo", "xthree"),
       Timecolname = "Env",
       Genofile = "your_working_directory/20000nuclearfamilies",
       SNPIDfile = paste0("your_working_directory/SNPs_group", 1:100, ".txt"),
       MinMAF = 0,
       Outputfile = paste0("your_working_directory/", Envtype, heritability, "/GALLOP.", nrep, ".txt"))


SAGELD_NULL_Model = SAGELD.NullModel(NullModel = nullmodel,
                                     PlinkFile = "your_working_directory/20000nuclearfamilies",
                                     SparseGRMFile = "your_working_directory/20000nuclear_SparseGRM_0.05.txt",
                                     PairwiseIBDFile = "your_working_directory/20000nuclear_PairwiseIBD.txt",
                                     control = list(ControlOutlier = FALSE))

GRAB.Marker(objNull = SAGELD_NULL_Model,
            GenoFile = "your_working_directory/20000nuclearfamilies.bed",
            OutputFile = paste0("your_working_directory/", Envtype, heritability, "/SAGELD.", nrep, ".txt"),
            control = list(SPA_Cutoff = 2,
                           zeta = 0,
                           tol = 1e-4,
                           impute_method = "mean",
                           missing_cutoff = 0.05,
                           min_maf_marker = 2e-4,
                           min_mac_marker = 20,
                           nMarkersEachChunk = 10000))
