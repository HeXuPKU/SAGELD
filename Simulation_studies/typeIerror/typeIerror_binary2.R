
# Please run this code on the slurm system
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
source("your_working_directory/SCEBE.R")
source("your_working_directory/SimulatePhenotype.R")

scr = 2; Envtype = "binary"

tables = data.frame(nrep = rep(1:100, times = 3),
                    medianObs = rep(c(2, 4, 8), each = 100))

nrep = tables$nrep[n.cpu]; medianObs = tables$medianObs[n.cpu]

nSub = 25000; nFam = 2500
randMat = data.table::fread("your_working_director/randMat10Members.txt") # pre-simulated random effects.
bvectomain = c(as.numeric(t(randMat[sample(1e6, nFam),])), rnorm(nSub))
bvectoGxE = c(as.numeric(t(randMat[sample(1e6, nFam),])), rnorm(nSub))

longpheno = pheno_simu(nSub = nSub, nFam = nFam, FamMode = "10-members", 
                       medianObs = medianObs, Env = "binary", MedFreq = 0.7, MedPrpo = 0.3, betaG = 1, betaGxE = 0,  
                       Geno = NULL, bVectomain = bvectomain, bVectoGxE = bvectoGxE, rho = 0.8)

nullmodel = lmer(pheno ~ xone + xtwo + xthree + Env + (Env|SubjID), data = longpheno, 
                 control = lmerControl(optimizer ='optimx', optCtrl = list(method = 'nlminb')))


gallop(nullmodel = nullmodel,
       Phenofile = longpheno,
       IDcolname = "SubjID",
       Phenocolname = "pheno",
       Covacolname = c("xone", "xtwo", "xthree"),
       Timecolname = "Env",
       Genofile = "your_working_director/nSub_25000_nFam_2500",
       SNPIDfile = paste0("your_working_director/SNPs_group", 1:100, ".txt"),
       MinMAF = 0,
       Outputfile = paste0("your_working_director", scr, "/", Envtype, medianObs, "/GALLOP.", nrep, ".txt"))


scebe(nullmodel = nullmodel,
      Phenofile = longpheno,
      IDcolname = "SubjID",
      Phenocolname = "pheno",
      Covacolname = c("xone", "xtwo", "xthree"),
      Timecolname = "Env",
      Genofile = "your_working_director/nSub_25000_nFam_2500",
      SNPIDfile = paste0("your_working_director/SNPs_group", 1:100, ".txt"),
      MinMAF = 0,
      Outputfile = paste0("your_working_director", scr, "/", Envtype, medianObs, "/SCEBE.", nrep, ".txt"))


SAGELD_NULL_Model = SAGELD.NullModel(NullModel = nullmodel,
                                     PlinkFile = "your_working_director/nSub_25000_nFam_2500",
                                     SparseGRMFile = "your_working_director/SparseGRM_0.05.txt",
                                     PairwiseIBDFile = "your_working_director/PairwiseIBD.txt")

GRAB.Marker(objNull = SAGELD_NULL_Model,
            GenoFile = "your_working_director/nSub_25000_nFam_2500.bed",
            OutputFile = paste0("your_working_director", scr, "/", Envtype, medianObs, "/SAGELD.", nrep, ".txt"),
            control = list(SPA_Cutoff = 2,
                           zeta = 0,
                           tol = 1e-4,
                           impute_method = "mean",
                           missing_cutoff = 0.05,
                           min_maf_marker = 2e-4,
                           min_mac_marker = 20,
                           nMarkersEachChunk = 10000))


GALLOP = data.table::fread(paste0("your_working_director", scr, "/", Envtype, medianObs, "/GALLOP.", nrep, ".txt"))
GALLOP = GALLOP %>% filter(LTpval < 5e-5)

if(nrow(GALLOP) != 0)
{
  temp_SNPIDs = paste0("your_working_director", scr, "/", Envtype, medianObs, "/Wald.", nrep, "SNPIDs.txt")
  
  data.table::fwrite(data.frame(GALLOP$SNPID), file = temp_SNPIDs,
                     row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  
  longpheno = longpheno %>% arrange(SubjID)
  frq = longpheno %>% group_by(SubjID) %>% summarize(n = n())
  
  GenoMat = GRAB::GRAB.ReadGeno(GenoFile = "your_working_director/nSub_25000_nFam_2500.bed",
                                SampleIDs = frq$SubjID,
                                control = list(IDsToIncludeFile = temp_SNPIDs,
                                               ImputeMethod = "mean"))
  
  Wald_results = c()
  
  for(i in 1:ncol(GenoMat$GenoMat))
  {
    cat(i, "\n")
    
    temp_longpheno = longpheno %>% mutate(Geno = rep(GenoMat$GenoMat[,i], times = frq$n))
    
    fullmodel = lmer(pheno ~ xone + xtwo + xthree + Env + Geno + Geno * Env + (Env|SubjID), data = temp_longpheno,
                     control = lmerControl(optimizer ='optimx', optCtrl = list(method = 'nlminb')))
    
    Wald_results = rbind(Wald_results,
                         data.frame(SNPID = colnames(GenoMat$GenoMat)[i],
                                    mainPvalue = summary(fullmodel)$coefficients[6,5],
                                    interPvalue = summary(fullmodel)$coefficients[7,5]))
  }
  
  data.table::fwrite(Wald_results,
                     file = paste0("your_working_director", scr, "/", Envtype, medianObs, "/Wald.", nrep, ".txt"),
                     row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  
  file.remove(temp_SNPIDs);
  rm(GenoMat); gc()
}


# cross-sectional analysis
reduce_pheno1 = longpheno %>% group_by(SubjID) %>% arrange(SubjID) %>%
  summarize(slope = mean(pheno[Env == 1]) - mean(pheno[Env == 0]),
            xone = first(xone), xtwo = first(xtwo), xthree = mean(xthree)) %>%
  mutate(FID = SubjID, IID = SubjID) %>%
  select(FID, IID, slope, xone, xtwo, xthree) %>% drop_na

reduce_pheno2 = longpheno %>% group_by(SubjID) %>% arrange(SubjID) %>%
  summarize(pheno = first(pheno), xone = first(xone), xtwo = first(xtwo), xthree = first(xthree), Env = first(Env)) %>%
  mutate(FID = SubjID, IID = SubjID) %>%
  select(FID, IID, pheno, xone, xtwo, xthree, Env)

randomslope = ranef(nullmodel)$SubjID
randomslope$SubjID = rownames(randomslope)

data.table::fwrite(reduce_pheno1 %>% select(FID, IID, slope),
                   file = paste0("your_working_director", scr, "/", Envtype, medianObs, "/Slope.", nrep, ".txt"),
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

data.table::fwrite(reduce_pheno1 %>% select(FID, IID, xone, xtwo, xthree),
                   file = paste0("your_working_director", scr, "/", Envtype, medianObs, "/Slope_cov.", nrep, ".txt"),
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

data.table::fwrite(reduce_pheno2 %>% select(FID, IID, pheno),
                   file = paste0("your_working_director", scr, "/", Envtype, medianObs, "/Crosssectional.", nrep, ".txt"),
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

data.table::fwrite(reduce_pheno2 %>% select(FID, IID, xone, xtwo, xthree),
                   file = paste0("your_working_director", scr, "/", Envtype, medianObs, "/Crosssectional_cov.", nrep, ".txt"),
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

data.table::fwrite(reduce_pheno2 %>% select(FID, IID, Env),
                   file = paste0("your_working_director", scr, "/", Envtype, medianObs, "/Crosssectional_Env.", nrep, ".txt"),
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

data.table::fwrite(reduce_pheno2 %>% select(FID, IID, pheno, xone, xtwo, xthree, Env),
                   file = paste0("your_working_director", scr, "/", Envtype, medianObs, "/Crosssectionalun.", nrep, ".txt"),
                   row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)

data.table::fwrite(randomslope %>% mutate(FID = SubjID, IID = SubjID) %>% select(FID, IID, Env),
                   file = paste0("your_working_director", scr, "/", Envtype, medianObs, "/RandomSlope.", nrep, ".txt"),
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

fastGWAcommoand1 = paste0("your_working_director/gcta-1.94.1 ",
                          "--bfile your_working_director/nSub_25000_nFam_2500 ",
                          "--grm-sparse your_working_director/fastGWA ",
                          "--fastGWA-mlm --thread-num 1 --nofilter ",
                          "--pheno your_working_director", scr, "/", Envtype, medianObs, "/Slope.", nrep, ".txt ",
                          "--qcovar your_working_director", scr, "/", Envtype, medianObs, "/Slope_cov.", nrep, ".txt ",
                          "--out your_working_director", scr, "/", Envtype, medianObs, "/Slope.", nrep)

system(fastGWAcommoand1)

fastGWAcommoand2 = paste0("your_working_director/gcta-1.94.1 ",
                          "--bfile your_working_director/nSub_25000_nFam_2500 ",
                          "--grm-sparse your_working_director/fastGWA ",
                          "--fastGWA-mlm --thread-num 1 --nofilter ",
                          "--pheno your_working_director", scr, "/", Envtype, medianObs, "/Crosssectional.", nrep, ".txt ",
                          "--qcovar your_working_director", scr, "/", Envtype, medianObs, "/Crosssectional_cov.", nrep, ".txt ",
                          "--envir your_working_director", scr, "/", Envtype, medianObs, "/Crosssectional_Env.", nrep, ".txt ",
                          "--out your_working_director", scr, "/", Envtype, medianObs, "/Crosssectional.", nrep)

system(fastGWAcommoand2)

fastGWAcommoand3 = paste0("your_working_director/gcta-1.94.1 ",
                          "--bfile your_working_director/nSub_25000_nFam_2500 ",
                          "--grm-sparse your_working_director/fastGWA ",
                          "--fastGWA-mlm --thread-num 1 --nofilter ",
                          "--pheno your_working_director", scr, "/", Envtype, medianObs, "/RandomSlope.", nrep, ".txt ",
                          "--out your_working_director", scr, "/", Envtype, medianObs, "/RandomSlope.", nrep)

system(fastGWAcommoand3)


file.remove(paste0("your_working_director", scr, "/", Envtype, medianObs, "/Slope.", nrep, ".txt"))
file.remove(paste0("your_working_director", scr, "/", Envtype, medianObs, "/Slope_cov.", nrep, ".txt"))
file.remove(paste0("your_working_director", scr, "/", Envtype, medianObs, "/Crosssectional.", nrep, ".txt"))
file.remove(paste0("your_working_director", scr, "/", Envtype, medianObs, "/Crosssectional_cov.", nrep, ".txt"))
file.remove(paste0("your_working_director", scr, "/", Envtype, medianObs, "/Crosssectional_Env.", nrep, ".txt"))
file.remove(paste0("your_working_director", scr, "/", Envtype, medianObs, "/Crosssectionalun.", nrep, ".txt"))
file.remove(paste0("your_working_director", scr, "/", Envtype, medianObs, "/RandomSlope.", nrep, ".txt"))
