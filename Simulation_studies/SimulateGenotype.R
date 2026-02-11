
library(dplyr)
library(tidyr)

# I used to simulate 100,000 common (MAF > 0.05) and rare variants (MAF < 0.05 & MAC > 20) for 50,000 related subjects 
# see also https://github.com/HeXuPKU/SPAGRM/tree/main/simulation/1.simulate_genotype
# Here I randomly extract 100,000 common variants with MAF > 0.01 from these two data sets.

### step1 merge this two plink file.
setwd("your_working_director")

mergelist = data.frame(c("your_working_director/nSub_25000_nFam_2500.bed",
                         "your_working_director/nSub_25000_nFam_2500.bed"),
                       c("your_working_director/nSub_25000_nFam_2500.bim",
                         "your_working_director/nSub_25000_nFam_2500.bim"),
                       c("your_working_director/nSub_25000_nFam_2500.fam",
                         "your_working_director/nSub_25000_nFam_2500.fam"))

data.table::fwrite(mergelist, "mergelist.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

plinkcommand = "plink --merge-list mergelist.txt --make-bed --maf 0.01 --out merged_data"
system(plinkcommand)

### step2 extract 100,000 variants from merged data.
SNPID = data.table::fread("your_working_director/merged_data.bim")

select_SNP = sample(SNPID$V2, 1e5)

data.table::fwrite(data.frame(select_SNP), "SNPIDs.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

plinkcommand = "plink --bfile merged_data --extract SNPIDs.txt --make-bed --out nSub_25000_nFam_2500"
system(plinkcommand)

plinkcommand = "plink --bfile nSub_25000_nFam_2500 --freq --out nSub_25000_nFam_2500"
system(plinkcommand)

file.remove("merged_data.bed")
file.remove("merged_data.bim")
file.remove("merged_data.fam")
file.remove("merged_data.log")

### step3 divide SNPS into 100 groups.
SNPIDs = data.table::fread("SNPIDs.txt", header = FALSE); SNPIDs = SNPIDs$V1

SNPIDsGroups = split(SNPIDs, cut(seq_along(SNPIDs), 100, labels = FALSE))

for(i in 1:100)
{
  cat(i, "\n")
  
  data.table::fwrite(data.frame(SNPIDsGroups[[i]]), paste0("SNPs_group", i, ".txt"), 
                     row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
}

