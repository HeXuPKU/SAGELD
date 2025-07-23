
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
BP = data.table::fread("your_working_directory/bp.txt")

# The gp scripts have been extracted. For details, see the code repository of https://www.nature.com/articles/s41467-025-58152-3
antihypertensive = data.table::fread("your_working_directory/antihypertensive.txt")

### merge longitudinal phenotype and medicaiton data, including
# Remove duplicated records in the same day;
# Treatment starts 7 days after the first dose;
# Treatment ends 6 months after the last dose.
# 
# betablocker = antihypertensive %>%
#   filter(type == "betablocker" | type == "") %>%
#   filter(grepl(paste0("acebutolol|atenolol|bisoprolol|carvedilol|celiprol|ollabetalol|",
#                       "metoprolol|nadolol|nebivolol|oxprenolol|pindolol|propranolol|",
#                       "sotalol|timolol|tenidone|inderal|tenoret|tenormin|tenif|monocor|",
#                       "labetalol|cardicor|kalten|adalat|celectol|nebilet|kalten|atenix|",
#                       "betaloc|emcor|sotacor|antipressan|trasidrex|lopresor|cardone|trasicor|sectral"), drug_name)) %>%
#   distinct(eid, issue_date, .keep_all = TRUE)
# 
# records = betablocker %>% group_by(eid) %>% summarize(n = n())
# 
# betablocker_treatment = c()
# 
# for(i in 1:nrow(records))
# {
#   if(i %% 1000 == 0) cat(i, "\n");
#   
#   if(records$n[i] == 1) next;
#   
#   temp_betablocker = betablocker %>% filter(eid == records$eid[i])
#   
#   med_start = temp_betablocker$issue_date[1]; med_end = temp_betablocker$issue_date[1]
#   
#   for(j in 2:records$n[i])
#   {
#     if(temp_betablocker$issue_date[j - 1] + 180 >= temp_betablocker$issue_date[j])
#     {
#       med_end = temp_betablocker$issue_date[j]
#       
#       if(j == records$n[i])
#       {
#         betablocker_treatment = rbind(betablocker_treatment,
#                                       data.frame(eid = records$eid[i],
#                                                  med_start = med_start + 7,
#                                                  med_end = med_end + 180))
#       }
#       
#     }else
#     {
#       med_end = temp_betablocker$issue_date[j - 1];
#       
#       if(med_end > med_start)
#       {
#         betablocker_treatment = rbind(betablocker_treatment,
#                                       data.frame(eid = records$eid[i],
#                                                  med_start = med_start + 7,
#                                                  med_end = med_end + 180))
#       }
#       
#       if(j < records$n[i])
#       {
#         med_start = temp_betablocker$issue_date[j]; med_end = temp_betablocker$issue_date[j]
#       }
#     }
#   }
# }
# 
# data.table::fwrite(betablocker_treatment,
#                    file = "your_working_directory/betablocker_treatment.txt",
#                    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
# 
# 
# ACEi = antihypertensive %>% 
#   filter(type == "ACEi" | type == "") %>%
#   filter(grepl(paste0("ramipril|lisinopril|perindopril|enalapril|captopril|trandolapril|",
#                       "fosinopril|quinapril|imidapril|moexipril|cilazapril|zestril|",
#                       "zestoretic|capozide|innovace|zidocapt|capoten|coversyl|staril|",
#                       "tritace|carace|accupro|gopten|innozide"), drug_name)) %>%
#   distinct(eid, issue_date, .keep_all = TRUE)
# 
# records = ACEi %>% group_by(eid) %>% summarize(n = n())
# 
# ACEi_treatment = c()
# 
# for(i in 1:nrow(records))
# {
#   if(i %% 1000 == 0) cat(i, "\n");
#   
#   if(records$n[i] == 1) next;
#   
#   temp_ACEi = ACEi %>% filter(eid == records$eid[i])
#   
#   med_start = temp_ACEi$issue_date[1]; med_end = temp_ACEi$issue_date[1]
#   
#   for(j in 2:records$n[i])
#   {
#     if(temp_ACEi$issue_date[j - 1] + 180 >= temp_ACEi$issue_date[j])
#     {
#       med_end = temp_ACEi$issue_date[j]
#       
#       if(j == records$n[i])
#       {
#         ACEi_treatment = rbind(ACEi_treatment,
#                                data.frame(eid = records$eid[i],
#                                           med_start = med_start + 7,
#                                           med_end = med_end + 180))
#       }
#       
#     }else
#     {
#       med_end = temp_ACEi$issue_date[j - 1];
#       
#       if(med_end > med_start)
#       {
#         ACEi_treatment = rbind(ACEi_treatment,
#                                data.frame(eid = records$eid[i],
#                                           med_start = med_start + 7,
#                                           med_end = med_end + 180))
#       }
#       
#       if(j < records$n[i])
#       {
#         med_start = temp_ACEi$issue_date[j]; med_end = temp_ACEi$issue_date[j]
#       }
#     }
#   }
# }
# 
# data.table::fwrite(ACEi_treatment,
#                    file = "your_working_directory/ACEi_treatment.txt",
#                    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
# 
# 
# ARB = antihypertensive %>%
#   filter(type == "ARB" | type == "") %>%
#   filter(grepl(paste0("losartan|candesartan|irbesartan|valsartan|telmisartan|olmesartan|",
#                       "eprosartan|cozaar|aprovel|diovan|amias|olmetec|micardis"), drug_name)) %>%
#   distinct(eid, issue_date, .keep_all = TRUE)
# 
# records = ARB %>% group_by(eid) %>% summarize(n = n())
# 
# ARB_treatment = c()
# 
# for(i in 1:nrow(records))
# {
#   if(i %% 1000 == 0) cat(i, "\n");
#   
#   if(records$n[i] == 1) next;
#   
#   temp_ARB = ARB %>% filter(eid == records$eid[i])
#   
#   med_start = temp_ARB$issue_date[1]; med_end = temp_ARB$issue_date[1]
#   
#   for(j in 2:records$n[i])
#   {
#     if(temp_ARB$issue_date[j - 1] + 180 >= temp_ARB$issue_date[j])
#     {
#       med_end = temp_ARB$issue_date[j]
#       
#       if(j == records$n[i])
#       {
#         ARB_treatment = rbind(ARB_treatment,
#                               data.frame(eid = records$eid[i],
#                                          med_start = med_start + 7,
#                                          med_end = med_end + 180))
#       }
#       
#     }else
#     {
#       med_end = temp_ARB$issue_date[j - 1];
#       
#       if(med_end > med_start)
#       {
#         ARB_treatment = rbind(ARB_treatment,
#                               data.frame(eid = records$eid[i],
#                                          med_start = med_start + 7,
#                                          med_end = med_end + 180))
#       }
#       
#       if(j < records$n[i])
#       {
#         med_start = temp_ARB$issue_date[j]; med_end = temp_ARB$issue_date[j]
#       }
#     }
#   }
# }
# 
# data.table::fwrite(ARB_treatment,
#                    file = "your_working_directory/ARB_treatment.txt",
#                    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
# 
# 
# CCB = antihypertensive %>%
#   filter(type == "CCB" | type == "") %>%
#   filter(grepl(paste0("amlodipine|nifedipine|diltiazem|felodipine|lercanidipine|",
#                       "lacidipine|nicardipine|isradipine|mibefradil|adalat|tildiem|",
#                       "adizem|coracten|adipine|cardioplen|coracten|dilzem|istin|",
#                       "slozem|angitil|viazem|motens|zemtard|zanidip|plendil|amlostin|",
#                       "tensipine|cardene|plendil|nifedipress|fortipine|dilcardia|angiopine|",
#                       "vascalpha|calchan|nifensar|nivaten"), drug_name)) %>%
#   distinct(eid, issue_date, .keep_all = TRUE)
# 
# records = CCB %>% group_by(eid) %>% summarize(n = n())
# 
# CCB_treatment = c()
# 
# for(i in 1:nrow(records))
# {
#   if(i %% 1000 == 0) cat(i, "\n");
# 
#   if(records$n[i] == 1) next;
# 
#   temp_CCB = CCB %>% filter(eid == records$eid[i])
# 
#   med_start = temp_CCB$issue_date[1]; med_end = temp_CCB$issue_date[1]
# 
#   for(j in 2:records$n[i])
#   {
#     if(temp_CCB$issue_date[j - 1] + 180 >= temp_CCB$issue_date[j])
#     {
#       med_end = temp_CCB$issue_date[j]
# 
#       if(j == records$n[i])
#       {
#         CCB_treatment = rbind(CCB_treatment,
#                               data.frame(eid = records$eid[i],
#                                          med_start = med_start + 7,
#                                          med_end = med_end + 180))
#       }
# 
#     }else
#     {
#       med_end = temp_CCB$issue_date[j - 1];
# 
#       if(med_end > med_start)
#       {
#         CCB_treatment = rbind(CCB_treatment,
#                               data.frame(eid = records$eid[i],
#                                          med_start = med_start + 7,
#                                          med_end = med_end + 180))
#       }
# 
#       if(j < records$n[i])
#       {
#         med_start = temp_CCB$issue_date[j]; med_end = temp_CCB$issue_date[j]
#       }
#     }
#   }
# }
# 
# data.table::fwrite(CCB_treatment,
#                    file = "your_working_directory/CCB_treatment.txt",
#                    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
# 
# 
# thiazide_diuretics = antihypertensive %>%
#   filter(type == "thiazide_diuretics" | type == "") %>%
#   filter(grepl(paste0("bendroflumethiazide|indapamide|xipamide|chlortalidone|hydrochlorothiazide|",
#                       "chlorothiazide|metolazone|mefruside|cyclopenthiazide|bendrofluazide|",
#                       "natrilix|aprinox|diurexan|nindaxa|saluric|berkozide|hygroton"), drug_name)) %>%
#   distinct(eid, issue_date, .keep_all = TRUE)
# 
# records = thiazide_diuretics %>% group_by(eid) %>% summarize(n = n())
# 
# thiazide_diuretics_treatment = c()
# 
# for(i in 1:nrow(records))
# {
#   if(i %% 1000 == 0) cat(i, "\n");
# 
#   if(records$n[i] == 1) next;
# 
#   temp_thiazide_diuretics = thiazide_diuretics %>% filter(eid == records$eid[i])
# 
#   med_start = temp_thiazide_diuretics$issue_date[1]; med_end = temp_thiazide_diuretics$issue_date[1]
# 
#   for(j in 2:records$n[i])
#   {
#     if(temp_thiazide_diuretics$issue_date[j - 1] + 180 >= temp_thiazide_diuretics$issue_date[j])
#     {
#       med_end = temp_thiazide_diuretics$issue_date[j]
# 
#       if(j == records$n[i])
#       {
#         thiazide_diuretics_treatment = rbind(thiazide_diuretics_treatment,
#                                              data.frame(eid = records$eid[i],
#                                                         med_start = med_start + 7,
#                                                         med_end = med_end + 180))
#       }
# 
#     }else
#     {
#       med_end = temp_thiazide_diuretics$issue_date[j - 1];
# 
#       if(med_end > med_start)
#       {
#         thiazide_diuretics_treatment = rbind(thiazide_diuretics_treatment,
#                                              data.frame(eid = records$eid[i],
#                                                         med_start = med_start + 7,
#                                                         med_end = med_end + 180))
#       }
# 
#       if(j < records$n[i])
#       {
#         med_start = temp_thiazide_diuretics$issue_date[j]; med_end = temp_thiazide_diuretics$issue_date[j]
#       }
#     }
#   }
# }
# 
# data.table::fwrite(thiazide_diuretics_treatment,
#                    file = "your_working_directory/thiazide_diuretics_treatment.txt",
#                    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
# 
# 
# loop_diuretics = antihypertensive %>%
#   filter(type == "loop_diuretics" | type == "") %>%
#   filter(grepl(paste0("furosemide|bumetanide|torasemide|piretanide|amilofruse|burinex|", 
#                       "frumil|frusemide|lasix|froop|lasoride"), drug_name)) %>%
#   distinct(eid, issue_date, .keep_all = TRUE)
# 
# records = loop_diuretics %>% group_by(eid) %>% summarize(n = n())
# 
# loop_diuretics_treatment = c()
# 
# for(i in 1:nrow(records))
# {
#   if(i %% 1000 == 0) cat(i, "\n");
#   
#   if(records$n[i] == 1) next;
#   
#   temp_loop_diuretics = loop_diuretics %>% filter(eid == records$eid[i])
#   
#   med_start = temp_loop_diuretics$issue_date[1]; med_end = temp_loop_diuretics$issue_date[1]
#   
#   for(j in 2:records$n[i])
#   {
#     if(temp_loop_diuretics$issue_date[j - 1] + 180 >= temp_loop_diuretics$issue_date[j])
#     {
#       med_end = temp_loop_diuretics$issue_date[j]
#       
#       if(j == records$n[i])
#       {
#         loop_diuretics_treatment = rbind(loop_diuretics_treatment,
#                                              data.frame(eid = records$eid[i],
#                                                         med_start = med_start + 7,
#                                                         med_end = med_end + 180))
#       }
#       
#     }else
#     {
#       med_end = temp_loop_diuretics$issue_date[j - 1];
#       
#       if(med_end > med_start)
#       {
#         loop_diuretics_treatment = rbind(loop_diuretics_treatment,
#                                              data.frame(eid = records$eid[i],
#                                                         med_start = med_start + 7,
#                                                         med_end = med_end + 180))
#       }
#       
#       if(j < records$n[i])
#       {
#         med_start = temp_loop_diuretics$issue_date[j]; med_end = temp_loop_diuretics$issue_date[j]
#       }
#     }
#   }
# }
# 
# data.table::fwrite(loop_diuretics_treatment,
#                    file = "your_working_directory/loop_diuretics_treatment.txt",
#                    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
# 
# 
# potassium_sparing_diuretics = antihypertensive %>%
#   filter(type == "pot_sparing_diuretics" | type == "") %>%
#   filter(grepl("spironolactone|eplerenone|amiloride|aldactone|spiroctan|triamterene", drug_name)) %>%
#   distinct(eid, issue_date, .keep_all = TRUE)
# 
# diuretics = rbind(thiazide_diuretics, loop_diuretics, potassium_sparing_diuretics) %>%
#   arrange(eid, issue_date) %>% distinct(eid, issue_date, .keep_all = TRUE)
# 
# records = diuretics %>% group_by(eid) %>% summarize(n = n())
# 
# diuretics_treatment = c()
# 
# for(i in 1:nrow(records))
# {
#   if(i %% 1000 == 0) cat(i, "\n");
# 
#   if(records$n[i] == 1) next;
# 
#   temp_diuretics = diuretics %>% filter(eid == records$eid[i])
# 
#   med_start = temp_diuretics$issue_date[1]; med_end = temp_diuretics$issue_date[1]
# 
#   for(j in 2:records$n[i])
#   {
#     if(temp_diuretics$issue_date[j - 1] + 180 >= temp_diuretics$issue_date[j])
#     {
#       med_end = temp_diuretics$issue_date[j]
# 
#       if(j == records$n[i])
#       {
#         diuretics_treatment = rbind(diuretics_treatment,
#                                              data.frame(eid = records$eid[i],
#                                                         med_start = med_start + 7,
#                                                         med_end = med_end + 180))
#       }
# 
#     }else
#     {
#       med_end = temp_diuretics$issue_date[j - 1];
# 
#       if(med_end > med_start)
#       {
#         diuretics_treatment = rbind(diuretics_treatment,
#                                              data.frame(eid = records$eid[i],
#                                                         med_start = med_start + 7,
#                                                         med_end = med_end + 180))
#       }
# 
#       if(j < records$n[i])
#       {
#         med_start = temp_diuretics$issue_date[j]; med_end = temp_diuretics$issue_date[j]
#       }
#     }
#   }
# }
# 
# data.table::fwrite(diuretics_treatment,
#                    file = "your_working_directory/diuretics_treatment.txt",
#                    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
# 
# 
# bp_combined = rbind(betablocker, ACEi, ARB, CCB, diuretics) %>%
#   arrange(eid, issue_date) %>% distinct(eid, issue_date, .keep_all = TRUE)
# 
# records = bp_combined %>% group_by(eid) %>% summarize(n = n())
# 
# bp_combined_treatment = c()
# 
# for(i in 1:nrow(records))
# {
#   if(i %% 1000 == 0) cat(i, "\n");
# 
#   if(records$n[i] == 1) next;
# 
#   temp_bp_combined = bp_combined %>% filter(eid == records$eid[i])
# 
#   med_start = temp_bp_combined$issue_date[1]; med_end = temp_bp_combined$issue_date[1]
# 
#   for(j in 2:records$n[i])
#   {
#     if(temp_bp_combined$issue_date[j - 1] + 180 >= temp_bp_combined$issue_date[j])
#     {
#       med_end = temp_bp_combined$issue_date[j]
# 
#       if(j == records$n[i])
#       {
#         bp_combined_treatment = rbind(bp_combined_treatment,
#                                       data.frame(eid = records$eid[i],
#                                                  med_start = med_start + 7,
#                                                  med_end = med_end + 180))
#       }
# 
#     }else
#     {
#       med_end = temp_bp_combined$issue_date[j - 1];
# 
#       if(med_end > med_start)
#       {
#         bp_combined_treatment = rbind(bp_combined_treatment,
#                                       data.frame(eid = records$eid[i],
#                                                  med_start = med_start + 7,
#                                                  med_end = med_end + 180))
#       }
# 
#       if(j < records$n[i])
#       {
#         med_start = temp_bp_combined$issue_date[j]; med_end = temp_bp_combined$issue_date[j]
#       }
#     }
#   }
# }
# 
# data.table::fwrite(bp_combined_treatment,
#                    file = "your_working_directory/bp_combined_treatment.txt",
#                    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


# antihypertensive drugs (combined)
BP_combined_treatment = data.table::fread("your_working_directory/bp_combined_treatment.txt")

BP_treatment = add_treatment(gp_clinical = BP,
                             gp_script = BP_combined_treatment)

data.table::fwrite(BP_treatment,
                   file = "your_working_directory/BP_medication.txt",
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
