
###### part1 preconditioning ######
# library(dplyr)
# library(tidyr)
# library(lme4)
# library(optimx)
# library(GRAB)
# 
# longpheno = data.table::fread("your_working_directory/BP.txt")
# 
# SubjIDinPLINK = data.table::fread("your_working_directory/markers_for_SPAGRMGE_lambda.fam")
# 
# genofile = data.table::fread("your_working_directory/ukb22828_c1_b0_v3_s487203.sample")
# IDs_have_genotypes = genofile$ID_1
# 
# longpheno = longpheno %>%
#   drop_na() %>% filter(bmi < 40) %>%
#   filter(feid %in% SubjIDinPLINK$V1) %>%
#   filter(feid %in% IDs_have_genotypes) %>%
#   select(feid, sbp, dbp, pp, age, sex_genetic, medication, PC1:PC10, bmi) %>%
#   mutate(sbp = scale(sbp), age = scale(age), age2 = scale(age^2), bmi = scale(bmi)) %>%
#   mutate(age_sex = age * sex_genetic, age_med = age * medication)
# 
# nullmodel = lmer(sbp ~ sex_genetic + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + medication + age + age2 + age_sex + age_med + bmi + (bmi | feid),
#                  data = longpheno, control = lmerControl(optimizer ='optimx', optCtrl = list(method = 'nlminb')))
# 
# SAGELD_NULL_Model = SAGELD.NullModel(NullModel = nullmodel,
#                                      PlinkFile = "your_working_directory/markers_for_SPAGRMGE_lambda",
#                                      SparseGRMFile = "your_working_directory/UKB410K_sparseGRM_0.05.txt",
#                                      PairwiseIBDFile = "your_working_directory/UKB410K_pairwiseIBD.txt")
# 
# save(SAGELD_NULL_Model,
#      file = "your_working_directory/SBP_bmi.RData")


###### part2 analyzing ######
# Please run this code on the slurm system
args=commandArgs(TRUE)
print(args)
print(sessionInfo())
n.cpu=as.numeric(args)
print(n.cpu)

cat("This is BMI.")

library(GRAB)

GenoFile = paste0("your_working_directory/ukb22828_c", n.cpu, "_b0_v3.bgen")
bgenbgiFile = paste0("your_working_directory/ukb_imp_chr", n.cpu, "_v3.bgen.bgi")
bgensampleFile = "your_working_directory/ukb22828_c1_b0_v3_s487203.sample"

OutputDir = paste0("your_working_directory/SAGELD", n.cpu)

SNPprefilter = data.table::fread(paste0("your_working_directory/ukb_mfi_chr", n.cpu, "_v3.txt"))
SNPprefilter = SNPprefilter %>% filter(V6 > 5e-3 & V8 >= 0.6)

IDsToIncludeFile = paste0(OutputDir, ".SNPprefilter.txt")

data.table::fwrite(data.table::data.table(SNPprefilter$V2), file = IDsToIncludeFile,
                   col.names = F, row.names = F, quote = F, sep = "\t")

precond.SAGELD = "your_working_directory/SBP_bmi.RData"
load(precond.SAGELD)

GRAB.Marker(objNull = SAGELD_NULL_Model,
            GenoFile = GenoFile,
            GenoFileIndex = c(bgenbgiFile, bgensampleFile),
            OutputFile = paste0(OutputDir, ".txt"),
            control = list(IDsToIncludeFile = IDsToIncludeFile,
                           SPA_Cutoff = 2,
                           zeta = 0,
                           tol = 1e-4,
                           impute_method = "mean",
                           missing_cutoff = 0.05,
                           min_maf_marker = 0.01,
                           min_mac_marker = 20,
                           nMarkersEachChunk = 10000))

file.remove(IDsToIncludeFile)