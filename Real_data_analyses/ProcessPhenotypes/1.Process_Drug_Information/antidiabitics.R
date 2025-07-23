
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)

add_treatment = function(gp_clinical, # cleaned measurements.
                         gp_script)   # cleaned medication.
{
  # take LDL and statins as example.
  LDL_records = gp_clinical %>% group_by(eid) %>% summarize(n = n())
  
  LDL_treat = c()
  
  for(i in 1:nrow(LDL_records))
  {
    if(i %% 1000 == 0) cat(i, "\n");
    
    if(LDL_records$eid[i] %in% gp_script$eid)
    {
      temp_LDL = gp_clinical %>% filter(eid == LDL_records$eid[i])
      
      temp_statins = gp_script %>% filter(eid == LDL_records$eid[i])
      
      for(j in 1:LDL_records$n[i])
      {
        temp_med = c();
        
        for(k in 1:nrow(temp_statins))
        {
          if(temp_LDL$event_dt[j] >= temp_statins$med_start[k] & temp_LDL$event_dt[j] <= temp_statins$med_end[k])
          {
            temp_med = c(temp_med, TRUE);
            break;
            
          }else
          {
            temp_med = c(temp_med, FALSE)
          }
        }
        
        temp_med = as.numeric(any(temp_med))
        
        LDL_treat = c(LDL_treat, temp_med)
      }
      
    }else
    {
      LDL_treat = c(LDL_treat, rep(0, times = LDL_records$n[i]))
    }
  }
  
  return(gp_clinical %>% mutate(medication = LDL_treat))
}


### read in longitudinal phenotype and medication data.
# The gp clinical data have been extracted. For details, see the code repository of https://www.nature.com/articles/s41467-025-56669-1
randomglu = data.table::fread("your_working_directory/RBG.txt")

hba1c = data.table::fread("your_working_directory/HbA1c.txt")

# The gp scripts have been extracted. For details, see the code repository of https://www.nature.com/articles/s41467-025-58152-3
Diabetes_lowering_drugs = data.table::fread("your_working_directory/antidiabetics.txt")

### merge longitudinal phenotype and medicaiton data, including
# select metformin/insulin/sulfonylureas related drugs.
# Remove duplicated records in the same day;
# Treatment starts 7 days after the first dose;
# Treatment ends 6 months after the last dose.
# 
# metformin = Diabetes_lowering_drugs %>%
#   filter(category == "biguanides" | category == "") %>%
#   filter(grepl(paste0("Metformin|Glucophage|Orabet|Glyformin|Ledermetin|Glucamet|", 
#                       "Milform|Metsol|Bolamyn|Metabet|Glucient|Diagemet|Sukkarto|Meijumet|Yaltormin|Metuxtan"), drug_name, ignore.case = TRUE)) %>%
#   distinct(eid, issue_date, .keep_all = TRUE)
# 
# records = metformin %>% group_by(eid) %>% summarize(n = n())
# 
# metformin_treatment = c()
# 
# for(i in 1:nrow(records))
# {
#   if(i %% 1000 == 0) cat(i, "\n");
# 
#   if(records$n[i] == 1) next;
# 
#   temp_metformin = metformin %>% filter(eid == records$eid[i])
# 
#   med_start = temp_metformin$issue_date[1]; med_end = temp_metformin$issue_date[1]
# 
#   for(j in 2:records$n[i])
#   {
#     if(temp_metformin$issue_date[j - 1] + 180 >= temp_metformin$issue_date[j])
#     {
#       med_end = temp_metformin$issue_date[j]
# 
#       if(j == records$n[i])
#       {
#         metformin_treatment = rbind(metformin_treatment,
#                                     data.frame(eid = records$eid[i],
#                                                med_start = med_start + 7,
#                                                med_end = med_end + 180))
#       }
# 
#     }else
#     {
#       med_end = temp_metformin$issue_date[j - 1];
# 
#       if(med_end > med_start)
#       {
#         metformin_treatment = rbind(metformin_treatment,
#                                     data.frame(eid = records$eid[i],
#                                                med_start = med_start + 7,
#                                                med_end = med_end + 180))
#       }
# 
#       if(j < records$n[i])
#       {
#         med_start = temp_metformin$issue_date[j]; med_end = temp_metformin$issue_date[j]
#       }
#     }
#   }
# }
# 
# data.table::fwrite(metformin_treatment,
#                    file = "your_working_directory/metformin_treatment.txt",
#                    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
# 
# 
# insulin = Diabetes_lowering_drugs %>%
#   filter(category %in% c("insulin", "insulin_pen", "insulin_short", "")) %>%
#   filter(grepl(paste0("insulin|lantus|novomix|novofine|glargine|levemir|humulin|", 
#                       "insulatard|mixtard|humalog|tresiba|xultophy|insuman|degludec|", 
#                       "velosulin|lispro|novorapid|actrapid|ultratard|apidra|monotard|",
#                       "hypurin|protaphane|toujeo|penmix"), drug_name, ignore.case = TRUE)) %>%
#   distinct(eid, issue_date, .keep_all = TRUE)
# 
# records = insulin %>% group_by(eid) %>% summarize(n = n())
# 
# insulin_treatment = c()
# 
# for(i in 1:nrow(records))
# {
#   if(i %% 1000 == 0) cat(i, "\n");
#   
#   if(records$n[i] == 1) next;
#   
#   temp_insulin = insulin %>% filter(eid == records$eid[i])
#   
#   med_start = temp_insulin$issue_date[1]; med_end = temp_insulin$issue_date[1]
#   
#   for(j in 2:records$n[i])
#   {
#     if(temp_insulin$issue_date[j - 1] + 180 >= temp_insulin$issue_date[j])
#     {
#       med_end = temp_insulin$issue_date[j]
#       
#       if(j == records$n[i])
#       {
#         insulin_treatment = rbind(insulin_treatment,
#                                   data.frame(eid = records$eid[i],
#                                              med_start = med_start + 7,
#                                              med_end = med_end + 180))
#       }
#       
#     }else
#     {
#       med_end = temp_insulin$issue_date[j - 1];
#       
#       if(med_end > med_start)
#       {
#         insulin_treatment = rbind(insulin_treatment,
#                                   data.frame(eid = records$eid[i],
#                                              med_start = med_start + 7,
#                                              med_end = med_end + 180))
#       }
#       
#       if(j < records$n[i])
#       {
#         med_start = temp_insulin$issue_date[j]; med_end = temp_insulin$issue_date[j]
#       }
#     }
#   }
# }
# 
# data.table::fwrite(insulin_treatment,
#                    file = "your_working_directory/insulin_treatment.txt",
#                    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
# 
# 
# sulfonylureas = Diabetes_lowering_drugs %>%
#   filter(category %in% c("sulfonylureas", "")) %>%
#   filter(grepl(paste0("gliclazide|glimepiride|glibenclamide|tolbutamide|glipizide|", 
#                       "diamicron|amaryl|daonil|chlorpropamide|diabinese|minodiab|",
#                       "tolazamide|tolanase|diabetamide|euglucon"), drug_name, ignore.case = TRUE)) %>%
#   distinct(eid, issue_date, .keep_all = TRUE)
# 
# records = sulfonylureas %>% group_by(eid) %>% summarize(n = n())
# 
# sulfonylureas_treatment = c()
# 
# for(i in 1:nrow(records))
# {
#   if(i %% 1000 == 0) cat(i, "\n");
#   
#   if(records$n[i] == 1) next;
#   
#   temp_sulfonylureas = sulfonylureas %>% filter(eid == records$eid[i])
#   
#   med_start = temp_sulfonylureas$issue_date[1]; med_end = temp_sulfonylureas$issue_date[1]
#   
#   for(j in 2:records$n[i])
#   {
#     if(temp_sulfonylureas$issue_date[j - 1] + 180 >= temp_sulfonylureas$issue_date[j])
#     {
#       med_end = temp_sulfonylureas$issue_date[j]
#       
#       if(j == records$n[i])
#       {
#         sulfonylureas_treatment = rbind(sulfonylureas_treatment,
#                                         data.frame(eid = records$eid[i],
#                                                    med_start = med_start + 7,
#                                                    med_end = med_end + 180))
#       }
#       
#     }else
#     {
#       med_end = temp_sulfonylureas$issue_date[j - 1];
#       
#       if(med_end > med_start)
#       {
#         sulfonylureas_treatment = rbind(sulfonylureas_treatment,
#                                         data.frame(eid = records$eid[i],
#                                                    med_start = med_start + 7,
#                                                    med_end = med_end + 180))
#       }
#       
#       if(j < records$n[i])
#       {
#         med_start = temp_sulfonylureas$issue_date[j]; med_end = temp_sulfonylureas$issue_date[j]
#       }
#     }
#   }
# }
# 
# data.table::fwrite(sulfonylureas_treatment,
#                    file = "your_working_directory/sulfonylureas_treatment.txt",
#                    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
# 
# 
# DPP4_inhibitors = Diabetes_lowering_drugs %>%
#   filter(category == "other" | category == "") %>%
#   filter(grepl("saxagliptin|sitagliptin|januvia|linagliptin|alogliptin|vildagliptin|onglyza", drug_name, ignore.case = TRUE)) %>%
#   distinct(eid, issue_date, .keep_all = TRUE)
# 
# GLP1_agonists = Diabetes_lowering_drugs %>%
#   filter(category == "other" | category == "") %>%
#   filter(grepl("exenatide|byetta|bydureon|liraglutide|victoza|dulaglutide|trulicity|lixisenatide|lyxumia|albiglutide", drug_name, ignore.case = TRUE)) %>%
#   distinct(eid, issue_date, .keep_all = TRUE)
# 
# TZDs = Diabetes_lowering_drugs %>%
#   filter(category == "other" | category == "") %>%
#   filter(grepl("pioglitazone|rosiglitazone", drug_name, ignore.case = TRUE)) %>%
#   distinct(eid, issue_date, .keep_all = TRUE)
# 
# SGLT2_inhibitors = Diabetes_lowering_drugs %>%
#   filter(category == "other" | category == "") %>%
#   filter(grepl("dapagliflozin|forxiga|canagliflozin|invokana|empagliflozin|jardiance", drug_name, ignore.case = TRUE)) %>%
#   distinct(eid, issue_date, .keep_all = TRUE)
# 
# DM_combined = rbind(metformin, sulfonylureas, insulin, DPP4_inhibitors, GLP1_agonists, TZDs, SGLT2_inhibitors) %>%
#   arrange(eid, issue_date) %>% distinct(eid, issue_date, .keep_all = TRUE)
# 
# records = DM_combined %>% group_by(eid) %>% summarize(n = n())
# 
# DM_combined_treatment = c()
# 
# for(i in 1:nrow(records))
# {
#   if(i %% 1000 == 0) cat(i, "\n");
# 
#   if(records$n[i] == 1) next;
# 
#   temp_DM_combined = DM_combined %>% filter(eid == records$eid[i])
# 
#   med_start = temp_DM_combined$issue_date[1]; med_end = temp_DM_combined$issue_date[1]
# 
#   for(j in 2:records$n[i])
#   {
#     if(temp_DM_combined$issue_date[j - 1] + 180 >= temp_DM_combined$issue_date[j])
#     {
#       med_end = temp_DM_combined$issue_date[j]
# 
#       if(j == records$n[i])
#       {
#         DM_combined_treatment = rbind(DM_combined_treatment,
#                                       data.frame(eid = records$eid[i],
#                                                  med_start = med_start + 7,
#                                                  med_end = med_end + 180))
#       }
# 
#     }else
#     {
#       med_end = temp_DM_combined$issue_date[j - 1];
# 
#       if(med_end > med_start)
#       {
#         DM_combined_treatment = rbind(DM_combined_treatment,
#                                       data.frame(eid = records$eid[i],
#                                                  med_start = med_start + 7,
#                                                  med_end = med_end + 180))
#       }
# 
#       if(j < records$n[i])
#       {
#         med_start = temp_DM_combined$issue_date[j]; med_end = temp_DM_combined$issue_date[j]
#       }
#     }
#   }
# }
# 
# data.table::fwrite(DM_combined_treatment,
#                    file = "your_working_directory/DM_combined_treatment.txt",
#                    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


# Diabetes-lowering drugs (combined)
DM_combined_treatment = data.table::fread("your_working_directory/DM_combined_treatment.txt")

randomglu_treatment = add_treatment(gp_clinical = randomglu,
                                    gp_script = DM_combined_treatment)

data.table::fwrite(randomglu_treatment,
                   file = "your_working_directory/randomglu_medication.txt",
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

hba1c_treatment = add_treatment(gp_clinical = hba1c,
                                gp_script = DM_combined_treatment)

data.table::fwrite(hba1c_treatment,
                   file = "your_working_directory/hba1c_medication.txt",
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
