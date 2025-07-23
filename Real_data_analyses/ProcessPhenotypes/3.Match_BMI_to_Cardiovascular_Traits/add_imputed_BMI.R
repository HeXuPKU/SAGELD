
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(data.table)

impute_bmi = function(BMI,
                      trait)
{
  trait[, row_id := .I]
  
  temp_trait = trait %>% select(feid, age, row_id) 
  
  temp_trait[, `:=`(min_age = age - 1, max_age = age)]
  
  temp_result = BMI[temp_trait, on = .(feid, age >= min_age, age <= max_age), allow.cartesian = TRUE,
                    .(feid, trait_age = i.age, bmi_age = x.age, bmi = x.bmi, row_id = i.row_id)]
  
  temp_result[, age_diff := abs(trait_age - bmi_age)]
  
  temp_result = temp_result[order(age_diff, -bmi_age), .SD[1], by = .(row_id)]
  
  output = inner_join(trait, temp_result %>% select(feid, bmi, age_diff, row_id))
  
  return(output)
}

### BMI
BMI = data.table::fread("your_working_directory/BMI.txt")

BMI = BMI %>% select(feid, bmi, age)
setkey(BMI, feid, age)

### BP
BP = data.table::fread("your_working_directory/BP.txt")

final_BP = impute_bmi(BMI, BP)

data.table::fwrite(final_BP, 
                   file = "your_working_directory/BP.txt",
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


### LDL
LDL = data.table::fread("your_working_directory/LDL.txt")

final_LDL = impute_bmi(BMI, LDL)

data.table::fwrite(final_LDL, 
                   file = "your_working_directory/LDL.txt",
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


### TC
TC = data.table::fread("your_working_directory/TC.txt")

final_TC = impute_bmi(BMI, TC)

data.table::fwrite(final_TC, 
                   file = "your_working_directory/TC.txt",
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


### HbA1c
HbA1c = data.table::fread("your_working_directory/HbA1c.txt")

final_HbA1c = impute_bmi(BMI, HbA1c)

data.table::fwrite(final_HbA1c, 
                   file = "your_working_directory/HbA1c.txt",
                   row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
