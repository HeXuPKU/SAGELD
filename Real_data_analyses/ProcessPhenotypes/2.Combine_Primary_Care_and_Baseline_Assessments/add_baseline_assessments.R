
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
source("your_working_directory/extractPheno.R")

### extract the age.
raw_age = extractPhenoFromUKBB(fieldID = 21003)

cleaned_age = raw_age$fieldData %>% 
  pivot_longer(cols = - feid, names_to = "observation", values_to = "age") %>% drop_na %>%
  mutate(observation = gsub("f21003", "", observation)) 


### extract SBP, DBP, and PP (with QC)
raw_SBP = extractPhenoFromUKBB(fieldID = 4080)

raw_DBP = extractPhenoFromUKBB(fieldID = 4079)

# read in medication information. For details, see the code repository of https://www.nature.com/articles/s41467-019-09572-5
antihypertensive = data.table::fread("your_working_directory/antihypertensive.txt")

antihypertensive = antihypertensive %>% 
  rename("_0_0" = "antiht1", "_1_0" = "antiht2", "_2_0" = "antiht3", "_3_0" = "antiht4") %>%
  pivot_longer(cols = - feid, names_to = "observation", values_to = "medication") %>% drop_na

# QC for SBP, DBP and PP
cleaned_SBP = raw_SBP$fieldData %>% 
  select(feid, f4080_0_0, f4080_1_0, f4080_2_0, f4080_3_0) %>% 
  pivot_longer(cols = - feid, names_to = "observation", values_to = "sbp") %>% drop_na %>%
  mutate(observation = gsub("f4080", "", observation)) %>%
  inner_join(cleaned_age) %>%
  filter(sbp > (140 - 4 * 20) & sbp < (140 + 4 * 20))

cleaned_DBP = raw_DBP$fieldData %>% 
  select(feid, f4079_0_0, f4079_1_0, f4079_2_0, f4079_3_0) %>% 
  pivot_longer(cols = - feid, names_to = "observation", values_to = "dbp") %>% drop_na %>%
  mutate(observation = gsub("f4079", "", observation)) %>%
  inner_join(cleaned_age) %>%
  filter(dbp > (80 - 4 * 10.5) & dbp < (80 + 4 * 10.5))

cleaned_BP = inner_join(cleaned_SBP, cleaned_DBP) %>% mutate(pp = sbp - dbp) %>%
  filter(pp > (57.5 - 4 * 15) & pp < (57.5 + 4 * 15)) %>%
  left_join(antihypertensive) %>% mutate(medication = as.numeric(medication)) %>%
  select(feid, sbp, dbp, pp, age, sex_genetic, medication, 
         PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, UNRELATED_WB)

# combine with gp_clinical data
BP_clinical = data.table::fread("your_working_directory/BP_medication.txt")

BP_clinical = BP_clinical %>% select(eid, SBP, DBP, PP, age, sex_genetic, medication, 
                                     PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, UNRELATED_WB) %>%
  rename(feid = eid, sbp = SBP, dbp = DBP, pp = PP)

final_BP = rbind(cleaned_BP %>% mutate(gp = 0), BP_clinical %>% mutate(gp = 1)) %>% arrange(feid, age)

data.table::fwrite(final_BP, "your_working_directory/BP.txt",
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


### extract LDL, HDL, TC, and TG
raw_LDL = extractPhenoFromUKBB(fieldID = 30780)

raw_HDL = extractPhenoFromUKBB(fieldID = 30760)

raw_TC = extractPhenoFromUKBB(fieldID = 30690)

raw_TG = extractPhenoFromUKBB(fieldID = 30870)

# read in medication information
antihyperlipidemia = data.table::fread("your_working_directory/antihyperlipidemia.txt")

antihyperlipidemia = antihyperlipidemia %>% 
  rename("_0_0" = "antihl1", "_1_0" = "antihl2", "_2_0" = "antihl3", "_3_0" = "antihl4") %>%
  pivot_longer(cols = - feid, names_to = "observation", values_to = "medication") %>% drop_na

# QC for LDL, HDL, TC, TG
cleaned_LDL = raw_LDL$fieldData %>% 
  pivot_longer(cols = - feid, names_to = "observation", values_to = "ldl") %>% drop_na %>%
  mutate(observation = gsub("f30780", "", observation)) %>%
  inner_join(cleaned_age) %>%
  filter(ldl > (3.6 - 4 * 0.87) & ldl < (3.6 + 4 * 0.87)) %>%
  left_join(antihyperlipidemia) %>% mutate(medication = as.numeric(medication)) %>%
  inner_join(covariate) %>% select(feid, ldl, age, sex_genetic, medication, 
                                   PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, UNRELATED_WB)

cleaned_HDL = raw_HDL$fieldData %>% 
  pivot_longer(cols = - feid, names_to = "observation", values_to = "hdl") %>% drop_na %>%
  mutate(observation = gsub("f30760", "", observation)) %>%
  inner_join(cleaned_age) %>%
  filter(hdl > (1.45 - 4 * 0.38) & hdl < (1.45 + 4 * 0.38)) %>%
  left_join(antihyperlipidemia) %>% mutate(medication = as.numeric(medication)) %>%
  inner_join(covariate) %>% select(feid, hdl, age, sex_genetic, medication, 
                                   PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, UNRELATED_WB)

cleaned_TC = raw_TC$fieldData %>% 
  pivot_longer(cols = - feid, names_to = "observation", values_to = "tc") %>% drop_na %>%
  mutate(observation = gsub("f30690", "", observation)) %>%
  inner_join(cleaned_age) %>%
  filter(tc > (5.7 - 4 * 1.15) & tc < (5.7 + 4 * 1.15)) %>%
  left_join(antihyperlipidemia) %>% mutate(medication = as.numeric(medication)) %>%
  inner_join(covariate) %>% select(feid, tc, age, sex_genetic, medication, 
                                   PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, UNRELATED_WB)

cleaned_TG = raw_TG$fieldData %>% 
  pivot_longer(cols = - feid, names_to = "observation", values_to = "tri") %>% drop_na %>%
  mutate(observation = gsub("f30870", "", observation)) %>%
  inner_join(cleaned_age) %>%
  filter(tri > (1.75 - 4 * 1) & tri < (1.75 + 4 * 1)) %>%
  left_join(antihyperlipidemia) %>% mutate(medication = as.numeric(medication)) %>%
  inner_join(covariate) %>% select(feid, tri, age, sex_genetic, medication, 
                                   PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, UNRELATED_WB)

# combine with gp_clinical data
LDL_clinical = data.table::fread("your_working_directory/LDL.txt")

LDL_clinical = LDL_clinical %>% select(eid, LDL, age, sex_genetic, medication, PC1, PC2, PC3, 
                                       PC4, PC5, PC6, PC7, PC8, PC9, PC10, UNRELATED_WB) %>% 
  rename(feid = eid, ldl = LDL)

final_LDL = rbind(cleaned_LDL %>% mutate(gp = 0), LDL_clinical %>% mutate(gp = 1)) %>% arrange(feid, age)

data.table::fwrite(final_LDL, "your_working_directory/LDL.txt",
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


HDL_clinical = data.table::fread("your_working_directory/HDL.txt")

HDL_clinical = HDL_clinical %>% select(eid, HDL, age, sex_genetic, medication, PC1, PC2, PC3, 
                                       PC4, PC5, PC6, PC7, PC8, PC9, PC10, UNRELATED_WB) %>% 
  rename(feid = eid, hdl = HDL)

final_HDL = rbind(cleaned_HDL %>% mutate(gp = 0), HDL_clinical %>% mutate(gp = 1)) %>% arrange(feid, age)

data.table::fwrite(final_HDL, "your_working_directory/HDL.txt",
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


TC_clinical = data.table::fread("your_working_directory/TC.txt")

TC_clinical = TC_clinical %>% select(eid, TC, age, sex_genetic, medication, PC1, PC2, PC3, 
                                     PC4, PC5, PC6, PC7, PC8, PC9, PC10, UNRELATED_WB) %>% 
  rename(feid = eid, tc = TC)

final_TC = rbind(cleaned_TC %>% mutate(gp = 0), TC_clinical %>% mutate(gp = 1)) %>% arrange(feid, age)

data.table::fwrite(final_TC, "your_working_directory/TC.txt",
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


TG_clinical = data.table::fread("your_working_directory/TG.txt")

TG_clinical = TG_clinical %>% select(eid, tri, age, sex_genetic, medication, PC1, PC2, PC3, 
                                     PC4, PC5, PC6, PC7, PC8, PC9, PC10, UNRELATED_WB) %>% 
  rename(feid = eid)

final_TG = rbind(cleaned_TG %>% mutate(gp = 0), TG_clinical %>% mutate(gp = 1)) %>% arrange(feid, age)

data.table::fwrite(final_TG, "your_working_directory/TG.txt",
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


### extract RBG, HbA1c
raw_RBG = extractPhenoFromUKBB(fieldID = 30740)

raw_HbA1c = extractPhenoFromUKBB(fieldID = 30750)

# read in medication information
antihyperglycemia = data.table::fread("your_working_directory/antihyperglycemia.txt")

antihyperglycemia = antihyperglycemia %>% 
  rename("_0_0" = "antihg1", "_1_0" = "antihg2", "_2_0" = "antihg3", "_3_0" = "antihg4") %>%
  pivot_longer(cols = - feid, names_to = "observation", values_to = "medication") %>% drop_na

# QC for RBG, HbA1c
cleaned_RBG = raw_RBG$fieldData %>% 
  pivot_longer(cols = - feid, names_to = "observation", values_to = "rbg") %>% drop_na %>%
  mutate(observation = gsub("f30740", "", observation)) %>%
  inner_join(cleaned_age) %>%
  filter(rbg > (5.12 - 4 * 1.25) & rbg < (5.12 + 4 * 1.25)) %>%
  left_join(antihyperglycemia) %>% mutate(medication = as.numeric(medication)) %>%
  inner_join(covariate) %>% select(feid, rbg, age, sex_genetic, medication, 
                                   PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, UNRELATED_WB)

cleaned_HbA1c = raw_HbA1c$fieldData %>% 
  pivot_longer(cols = - feid, names_to = "observation", values_to = "a1c") %>% drop_na %>%
  mutate(observation = gsub("f30750", "", observation)) %>%
  inner_join(cleaned_age) %>%
  mutate(a1c = round(a1c/10.929 + 2.15, 1)) %>%
  filter(a1c > (5.5 - 4 * 0.6) & a1c < (5.5 + 4 * 0.6)) %>%
  left_join(antihyperglycemia) %>% mutate(medication = as.numeric(medication)) %>%
  inner_join(covariate) %>% select(feid, a1c, age, sex_genetic, medication, 
                                   PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, UNRELATED_WB)


# combine with gp_clinical data
RBG_clinical = data.table::fread("your_working_directory/randomglu_medication.txt")

RBG_clinical = RBG_clinical %>% select(eid, rglu, age, sex_genetic, medication, PC1, PC2, PC3, 
                                       PC4, PC5, PC6, PC7, PC8, PC9, PC10, UNRELATED_WB) %>% 
  rename(feid = eid, rbg = rglu)

final_RBG = rbind(cleaned_RBG %>% mutate(gp = 0), RBG_clinical %>% mutate(gp = 1)) %>% arrange(feid, age)

data.table::fwrite(final_RBG, "your_working_directory/RBG.txt",
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


HbA1c_clinical = data.table::fread("your_working_directory/hba1c_medication.txt")

HbA1c_clinical = HbA1c_clinical %>% select(eid, hba1c, age, sex_genetic, medication, PC1, PC2, PC3, 
                                           PC4, PC5, PC6, PC7, PC8, PC9, PC10, UNRELATED_WB) %>% 
  rename(feid = eid, a1c = hba1c)

final_HbA1c = rbind(cleaned_HbA1c %>% mutate(gp = 0), HbA1c_clinical %>% mutate(gp = 1)) %>% arrange(feid, age)

data.table::fwrite(final_HbA1c, "your_working_directory/HbA1c.txt",
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


### extract eGFR (with QC)
raw_eGFR = extractPhenoFromUKBB(fieldID = 30700)

# QC for eGFR
cleaned_eGFR = raw_eGFR$fieldData %>% 
  select(feid, f30700_0_0, f30700_1_0) %>% 
  pivot_longer(cols = - feid, names_to = "observation", values_to = "crea") %>% drop_na %>%
  mutate(observation = gsub("f30700", "", observation)) %>%
  inner_join(cleaned_age) %>%
  inner_join(covariate) %>% 
  rowwise() %>%
  # mutate(egfr = ifelse(sex_genetic == 1, 141*min(crea/88.4/0.9, 1)^(-0.411)*max(crea/88.4/0.9, 1)^(-1.209)*0.993^age,
  #                      141*min(crea/88.4/0.7, 1)^(-0.329)*max(crea/88.4/0.7, 1)^(-1.209)*0.993^age*1.018)) %>%
  mutate(egfr = ifelse(sex_genetic == 1, 142*min(crea/88.4/0.9, 1)^(-0.302)*max(crea/88.4/0.9, 1)^(-1.2)*0.9938^age,
                       142*min(crea/88.4/0.7, 1)^(-0.241)*max(crea/88.4/0.7, 1)^(-1.2)*0.9938^age*1.012)) %>%
  mutate(egfr = round(egfr, 1)) %>%
  filter(egfr > (90 - 4 * 15) & egfr < (90 + 4 * 15)) %>%
  select(feid, egfr, age, sex_genetic, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, UNRELATED_WB)

# combine with gp_clinical data
eGFR_clinical = data.table::fread("/gdata01/user/xuhe/SPA-GRM/real_data-2022-12-29/final_pheno/eGFR.txt")

eGFR_clinical = eGFR_clinical %>% 
  select(eid, crea, age, sex_genetic, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, UNRELATED_WB) %>%
  rename(feid = eid)

# final_eGFR = rbind(cleaned_eGFR, eGFR_clinical) %>% as.data.frame() %>% arrange(feid, age)

final_eGFR = rbind(raw_eGFR$fieldData %>% 
                     select(feid, f30700_0_0, f30700_1_0) %>% 
                     pivot_longer(cols = - feid, names_to = "observation", values_to = "crea") %>% drop_na %>%
                     mutate(observation = gsub("f30700", "", observation)) %>%
                     inner_join(cleaned_age) %>%
                     inner_join(covariate) %>% 
                     select(feid, crea, age, sex_genetic, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, UNRELATED_WB) %>% mutate(gp = 0),
                   eGFR_clinical %>% mutate(gp = 1)) %>%
  arrange(feid, age) %>%
  rowwise() %>% # CKD-EPI Creatinine Equation (2021)
  mutate(egfr = ifelse(sex_genetic == 1, 142*min(crea/88.4/0.9, 1)^(-0.302)*max(crea/88.4/0.9, 1)^(-1.2)*0.9938^age,
                       142*min(crea/88.4/0.7, 1)^(-0.241)*max(crea/88.4/0.7, 1)^(-1.2)*0.9938^age*1.012)) %>%
  ungroup %>%
  filter(egfr > 84.5 - 4 * 16.5 & egfr < 84.5 + 4 * 16.5)

data.table::fwrite(final_eGFR, "your_working_directory/eGFR.txt",
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


### extract BMI (with QC)
raw_BMI = extractPhenoFromUKBB(fieldID = 21001)

# QC for BMI
cleaned_BMI = raw_BMI$fieldData %>% 
  select(feid, f21001_0_0, f21001_1_0, f21001_2_0, f21001_3_0) %>% 
  pivot_longer(cols = - feid, names_to = "observation", values_to = "bmi") %>% drop_na %>%
  mutate(observation = gsub("f21001", "", observation)) %>%
  inner_join(cleaned_age) %>%
  inner_join(covariate) %>% 
  filter(bmi > (28 - 4 * 5) & bmi < (28 + 4 * 5)) %>%
  select(feid, bmi, age, sex_genetic, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, UNRELATED_WB)

# combine with gp_clinical data
BMI_clinical = data.table::fread("your_working_directory/bmi.txt")

BMI_clinical = BMI_clinical %>% 
  select(eid, bmi, age, sex_genetic, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, UNRELATED_WB) %>%
  rename(feid = eid)

final_BMI = rbind(cleaned_BMI %>% mutate(gp = 0), BMI_clinical %>% mutate(gp = 1)) %>% as.data.frame() %>% arrange(feid, age)

data.table::fwrite(final_BMI, "your_working_directory/BMI.txt",
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
