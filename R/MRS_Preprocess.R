# We will first pre-process the MRS data
# Then we will use that data in the Machine learning models

library(xlsx)
library(data.table)
library(feather)

# main xlsx sheet
path <- "G:/PGC ML/MRS/MRS_EWAS_sample_secondary_analysis_final_v3.xlsx"
pheno <- read.xlsx(path, sheetIndex = 4)
head(pheno[, 1:5])

sample_sheet <- read.xlsx(path, sheetIndex = 2)
head(sample_sheet)

table(pheno$EWAS_id %in% sample_sheet$EWAS_id)

mrs_pheno <- merge(sample_sheet, pheno, by="EWAS_id")
mrs_pheno$Study <- "MRS"
mrs_pheno$BaseName <- paste0(mrs_pheno$sentrix_id,"_",mrs_pheno$sentrix_position)


# check unique
length(unique(mrs_pheno$BaseName))
any(duplicated(mrs_pheno$BaseName))


# file with cell estimations
mrs_cell_est <- read.csv("G:/PGC ML/MRS/MRS_Pheno_PCs_CellTypes.csv")
mrs_pheno <- merge(mrs_pheno, mrs_cell_est,
                   by.x = "BaseName",
                   by.y = "SampleID")

# As some columns are common in both and while merging will add .x and .y to the columns.
# lets remove the columns with .y and then remove .x from the remaining columns
mrs_pheno <- mrs_pheno[, which(!grepl("\\.y$", colnames(mrs_pheno),
                                      ignore.case = T))]
dim(mrs_pheno)

# Now remove the .x part from the column names
colnames(mrs_pheno) <- gsub("\\.x$", "", colnames(mrs_pheno))

write.csv(mrs_pheno, "G:/PGC ML/MRS/MRS_Pheno.csv", row.names = F)


# We will also convert CSV methylation file to feather file
# So that we can read it fast in python
mrs_meth <- fread("G:/PGC ML/MRS/MRS_noob_qcd_crossReactiveProbesRemoved_combat_CP_wcovar_age_ptsd_allPreAsControls.csv",
                  data.table = F)

write_feather(mrs_meth, "G:/PGC ML/MRS/MRS_noob_qcd_crossReactiveProbesRemoved_combat_CP_wcovar_age_ptsd_allPreAsControls.feather")
