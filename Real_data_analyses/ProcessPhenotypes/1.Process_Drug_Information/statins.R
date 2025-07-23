
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
LDL = data.table::fread("your_working_directory/LDL.txt")

HDL = data.table::fread("your_working_directory/HDL.txt")

TC = data.table::fread("your_working_directory/totchol.txt")

TG = data.table::fread("your_working_directory/triglyceride.txt")

# The gp scripts have been extracted. For details, see the code repository of https://www.nature.com/articles/s41467-025-58152-3
statins = data.table::fread("your_working_directory/statins.txt")

### merge longitudinal phenotype and medicaiton data, including
# Remove duplicated records in the same day;
# Treatment starts 7 days after the first dose;
# Treatment ends 6 months after the last dose.
# 
# statins = statins %>% 
#   filter(drug != "") %>% 
#   distinct(eid, issue_date, .keep_all = TRUE)
# 
# records = statins %>% group_by(eid) %>% summarize(n = n())
# 
# statins_treatment = c()
# 
# for(i in 1:nrow(records))
# {
#   if(i %% 1000 == 0) cat(i, "\n");
# 
#   if(records$n[i] == 1) next;
# 
#   temp_statins = statins %>% filter(eid == records$eid[i])
# 
#   med_start = temp_statins$issue_date[1]; med_end = temp_statins$issue_date[1]
# 
#   for(j in 2:records$n[i])
#   {
#     if(temp_statins$issue_date[j - 1] + 180 >= temp_statins$issue_date[j])
#     {
#       med_end = temp_statins$issue_date[j]
# 
#       if(j == records$n[i])
#       {
#         statins_treatment = rbind(statins_treatment,
#                                   data.frame(eid = records$eid[i],
#                                              med_start = med_start + 7,
#                                              med_end = med_end + 180))
#       }
# 
#     }else
#     {
#       med_end = temp_statins$issue_date[j - 1];
# 
#       if(med_end > med_start)
#       {
#         statins_treatment = rbind(statins_treatment,
#                                   data.frame(eid = records$eid[i],
#                                              med_start = med_start + 7,
#                                              med_end = med_end + 180))
#       }
# 
#       if(j < records$n[i])
#       {
#         med_start = temp_statins$issue_date[j]; med_end = temp_statins$issue_date[j]
#       }
#     }
#   }
# }
# 
# data.table::fwrite(statins_treatment,
#                    file = "your_working_directory/statins_treatment.txt",
#                    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
# 
# 
# atorvastatin = statins %>%
#   filter(drug == "atorvastatin") %>%
#   distinct(eid, issue_date, .keep_all = TRUE)
# 
# records = atorvastatin %>% group_by(eid) %>% summarize(n = n())
# 
# atorvastatin_treatment = c()
# 
# for(i in 1:nrow(records))
# {
#   if(i %% 1000 == 0) cat(i, "\n");
# 
#   if(records$n[i] == 1) next;
# 
#   temp_atorvastatin = atorvastatin %>% filter(eid == records$eid[i])
# 
#   med_start = temp_atorvastatin$issue_date[1]; med_end = temp_atorvastatin$issue_date[1]
# 
#   for(j in 2:records$n[i])
#   {
#     if(temp_atorvastatin$issue_date[j - 1] + 180 >= temp_atorvastatin$issue_date[j])
#     {
#       med_end = temp_atorvastatin$issue_date[j]
# 
#       if(j == records$n[i])
#       {
#         atorvastatin_treatment = rbind(atorvastatin_treatment,
#                                        data.frame(eid = records$eid[i],
#                                                   med_start = med_start + 7,
#                                                   med_end = med_end + 180))
#       }
# 
#     }else
#     {
#       med_end = temp_atorvastatin$issue_date[j - 1];
# 
#       if(med_end > med_start)
#       {
#         atorvastatin_treatment = rbind(atorvastatin_treatment,
#                                        data.frame(eid = records$eid[i],
#                                                   med_start = med_start + 7,
#                                                   med_end = med_end + 180))
#       }
# 
#       if(j < records$n[i])
#       {
#         med_start = temp_atorvastatin$issue_date[j]; med_end = temp_atorvastatin$issue_date[j]
#       }
#     }
#   }
# }
# 
# data.table::fwrite(atorvastatin_treatment,
#                    file = "your_working_directory/atorvastatin_treatment.txt",
#                    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
# 
# 
# simvastatin = statins %>%
#   filter(drug == "simvastatin") %>%
#   distinct(eid, issue_date, .keep_all = TRUE)
# 
# records = simvastatin %>% group_by(eid) %>% summarize(n = n())
# 
# simvastatin_treatment = c()
# 
# for(i in 1:nrow(records))
# {
#   if(i %% 1000 == 0) cat(i, "\n");
#   
#   if(records$n[i] == 1) next;
#   
#   temp_simvastatin = simvastatin %>% filter(eid == records$eid[i])
#   
#   med_start = temp_simvastatin$issue_date[1]; med_end = temp_simvastatin$issue_date[1]
#   
#   for(j in 2:records$n[i])
#   {
#     if(temp_simvastatin$issue_date[j - 1] + 180 >= temp_simvastatin$issue_date[j])
#     {
#       med_end = temp_simvastatin$issue_date[j]
#       
#       if(j == records$n[i])
#       {
#         simvastatin_treatment = rbind(simvastatin_treatment,
#                                       data.frame(eid = records$eid[i],
#                                                  med_start = med_start + 7,
#                                                  med_end = med_end + 180))
#       }
#       
#     }else
#     {
#       med_end = temp_simvastatin$issue_date[j - 1];
#       
#       if(med_end > med_start)
#       {
#         simvastatin_treatment = rbind(simvastatin_treatment,
#                                       data.frame(eid = records$eid[i],
#                                                  med_start = med_start + 7,
#                                                  med_end = med_end + 180))
#       }
#       
#       if(j < records$n[i])
#       {
#         med_start = temp_simvastatin$issue_date[j]; med_end = temp_simvastatin$issue_date[j]
#       }
#     }
#   }
# }
# 
# data.table::fwrite(simvastatin_treatment,
#                    file = "your_working_directory/simvastatin_treatment.txt",
#                    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
# 
# 
# pravastatin = statins %>%
#   filter(drug == "pravastatin") %>%
#   distinct(eid, issue_date, .keep_all = TRUE)
# 
# records = pravastatin %>% group_by(eid) %>% summarize(n = n())
# 
# pravastatin_treatment = c()
# 
# for(i in 1:nrow(records))
# {
#   if(i %% 1000 == 0) cat(i, "\n");
#   
#   if(records$n[i] == 1) next;
#   
#   temp_pravastatin = pravastatin %>% filter(eid == records$eid[i])
#   
#   med_start = temp_pravastatin$issue_date[1]; med_end = temp_pravastatin$issue_date[1]
#   
#   for(j in 2:records$n[i])
#   {
#     if(temp_pravastatin$issue_date[j - 1] + 180 >= temp_pravastatin$issue_date[j])
#     {
#       med_end = temp_pravastatin$issue_date[j]
#       
#       if(j == records$n[i])
#       {
#         pravastatin_treatment = rbind(pravastatin_treatment,
#                                       data.frame(eid = records$eid[i],
#                                                  med_start = med_start + 7,
#                                                  med_end = med_end + 180))
#       }
#       
#     }else
#     {
#       med_end = temp_pravastatin$issue_date[j - 1];
#       
#       if(med_end > med_start)
#       {
#         pravastatin_treatment = rbind(pravastatin_treatment,
#                                       data.frame(eid = records$eid[i],
#                                                  med_start = med_start + 7,
#                                                  med_end = med_end + 180))
#       }
#       
#       if(j < records$n[i])
#       {
#         med_start = temp_pravastatin$issue_date[j]; med_end = temp_pravastatin$issue_date[j]
#       }
#     }
#   }
# }
# 
# data.table::fwrite(pravastatin_treatment,
#                    file = "your_working_directory/pravastatin_treatment.txt",
#                    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
# 
# 
# rosuvastatin = statins %>%
#   filter(drug == "rosuvastatin") %>%
#   distinct(eid, issue_date, .keep_all = TRUE)
# 
# records = rosuvastatin %>% group_by(eid) %>% summarize(n = n())
# 
# rosuvastatin_treatment = c()
# 
# for(i in 1:nrow(records))
# {
#   if(i %% 1000 == 0) cat(i, "\n");
#   
#   if(records$n[i] == 1) next;
#   
#   temp_rosuvastatin = rosuvastatin %>% filter(eid == records$eid[i])
#   
#   med_start = temp_rosuvastatin$issue_date[1]; med_end = temp_rosuvastatin$issue_date[1]
#   
#   for(j in 2:records$n[i])
#   {
#     if(temp_rosuvastatin$issue_date[j - 1] + 180 >= temp_rosuvastatin$issue_date[j])
#     {
#       med_end = temp_rosuvastatin$issue_date[j]
#       
#       if(j == records$n[i])
#       {
#         rosuvastatin_treatment = rbind(rosuvastatin_treatment,
#                                        data.frame(eid = records$eid[i],
#                                                   med_start = med_start + 7,
#                                                   med_end = med_end + 180))
#       }
#       
#     }else
#     {
#       med_end = temp_rosuvastatin$issue_date[j - 1];
#       
#       if(med_end > med_start)
#       {
#         rosuvastatin_treatment = rbind(rosuvastatin_treatment,
#                                        data.frame(eid = records$eid[i],
#                                                   med_start = med_start + 7,
#                                                   med_end = med_end + 180))
#       }
#       
#       if(j < records$n[i])
#       {
#         med_start = temp_rosuvastatin$issue_date[j]; med_end = temp_rosuvastatin$issue_date[j]
#       }
#     }
#   }
# }
# 
# data.table::fwrite(rosuvastatin_treatment,
#                    file = "your_working_directory/rosuvastatin_treatment.txt",
#                    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

# combined
statins_treatment = data.table::fread("your_working_directory/statins_treatment.txt")

LDL_treatment = add_treatment(gp_clinical = LDL,
                              gp_script = statins_treatment)

data.table::fwrite(LDL_treatment,
                   file = "your_working_directory/LDL.txt",
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

HDL_treatment = add_treatment(gp_clinical = HDL,
                              gp_script = statins_treatment)

data.table::fwrite(HDL_treatment,
                   file = "your_working_directory/HDL.txt",
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

TC_treatment = add_treatment(gp_clinical = TC,
                              gp_script = statins_treatment)

data.table::fwrite(TC_treatment,
                   file = "your_working_directory/TC.txt",
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

TG_treatment = add_treatment(gp_clinical = TG,
                              gp_script = statins_treatment)

data.table::fwrite(TG_treatment,
                   file = "your_working_directory/TG.txt",
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
