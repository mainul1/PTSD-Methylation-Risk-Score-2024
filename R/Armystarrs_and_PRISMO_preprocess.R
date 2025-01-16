
# We will first pre-process the Armystarrs data
# Then we will use that data in the Machine learning models

library(xlsx)
library(data.table)
library(feather)
library(tidyr)
library(dplyr)

# main xlsx sheet
armys_path <- "G:/PGC ML/ArmySTARRS/ArmySTARRS_PPDS_EWAS_secondary_analysis_final_v6.xlsx"
prismo_path <- "G:/PGC ML/PRISMO/PRISMO_EWAS_secondary_analysis_final_v4.xlsx"
paths <- c(armys_path, prismo_path)

pheno <- lapply(paths, function(x) read.xlsx(x, sheetIndex = 4))
df_names <- c("Armystarrs", "Prismo")
names(pheno) <- df_names

View(pheno$Armystarrs)
View(pheno$Prismo)

# create new id so that we can match tp 1 methylation data with tp 2 phenotype data
pheno$Armystarrs <- pheno$Armystarrs %>% mutate(EWAS_id = sub("\\_.*", "", EWAS_id))

# Some measurements are assessed only at tp 0, so we need to use those measures
# for other time points
army_cols <- c("MaltreatmentGlobal", "nondeploy_trauma_exposed_critA",
               "deploy_trauma_exposed_critA", "trauma_exposed_critA")

# We created new_id to group them and fill na using other timepoint

pheno$Armystarrs <- pheno$Armystarrs %>%
  arrange(visitkey) %>%
  group_by(EWAS_id) %>%
  fill(!!! army_cols) %>%
  mutate(n = n())

prismo_cols <- c("ETItot", "ETIalg", "ETIlich", "ETIgeest", "ETIsex", "Pes_number")
pheno$Prismo <- pheno$Prismo %>%
  arrange(visitkey) %>%
  group_by(EWAS_id) %>%
  fill(!!! prismo_cols) %>%
  mutate(n = n())

View(pheno$Armystarrs)
View(pheno$Prismo)


lapply(pheno, function(x) head(x[, 1:5]))
lapply(pheno, dim)

# In Armystarrs, specimens were collected at 0&1 timepoint (tp) only,
# but 1 has no phenotype data
# so we will use tp 1 methylation data and tp 2 phenotype data
pheno$Armystarrs <- pheno$Armystarrs %>% subset(visit==2) # phenotype data
pheno$Prismo <- pheno$Prismo %>% subset(visit=='2_epic') # phenotype data

# function to read the sheets and combine sentrix_id and position
sample_sheets <- lapply(seq_along(paths), function(i) {
  sheet <- read.xlsx(paths[i], sheetIndex = 2)
  if(df_names[i] == "Armystarrs"){
    sheet <- sheet %>% subset(visit== 1)
  }
  else{
    sheet <- sheet %>% subset(visit == '2_epic')
  }

  sheet <- sheet[!is.na(sheet$sentrix_id), ]
  sheet$BaseName <- paste0(sheet$sentrix_id,"_",sheet$sentrix_position)
  sheet$Study <- df_names[i]
  sheet
  })

names(sample_sheets) <- df_names
lapply(sample_sheets, function(x) head(x[, 1:5]))
lapply(sample_sheets, dim)


# now create new ids for armystarrs sample sheet too,
# so that we can match it with the phenotype
sample_sheets$Armystarrs <- sample_sheets$Armystarrs %>%
  mutate(EWAS_id = gsub("\\_.*", "", EWAS_id))
View(sample_sheets$Armystarrs)
View(sample_sheets$Prismo)

lapply(sample_sheets, function(x) length(unique(x$EWAS_id)))
lapply(pheno, function(x) length(unique(x$EWAS_id)))

table(sample_sheets$Armystarrs$EWAS_id %in% pheno$Armystarrs$EWAS_id)
table(sample_sheets$Prismo$EWAS_id %in% pheno$Prismo$EWAS_id)

# combine phenotypes and sample sheet
pheno_comb <- Map(function(x,y) merge(x,y, by="EWAS_id"), pheno, sample_sheets)
lapply(pheno_comb, dim)

View(pheno_comb$Armystarrs)
View(pheno_comb$Prismo)


View(pheno_comb$Armystarrs %>% select(contains(c("EWAS_id", "pcl"))))


# Some sentrix ids and positions are duplicated
# lets see how many
lapply(pheno_comb, function(x) length(which(duplicated(x$BaseName))))
duplic <- lapply(pheno_comb, function(x) x[duplicated(x$BaseName), ])
lapply(duplic, dim)
View(duplic$Armystarrs)
View(duplic$Prismo)




# Get the samples that have duplicates and no duplicates
all_duplicates <-  lapply(pheno_comb, function(x) x %>% group_by(BaseName) %>% filter(n()>1))
lapply(all_duplicates, dim)
lapply(all_duplicates, function(x) length(unique(x$BaseName))) # count unique from the set of duplicated


no_duplicates <- lapply(pheno_comb, function(x) x %>% group_by(BaseName) %>% filter(n() == 1))
lapply(no_duplicates, dim)
View(no_duplicates$Armystarrs)


# How many can be unique. The id can still be duplicated
uniq <- lapply(pheno_comb, function(x) x[!duplicated(x$BaseName), ])
lapply(uniq, dim)


# Check ewas ids are unique
lapply(pheno_comb, function(x) length(unique(x$EWAS_id)))


# Get only the samples on Epic array
epic_samps <- lapply(pheno_comb, function(x) {
  s <- x[grepl("epic", x$array, ignore.case = T), ]
  s <- s[!is.na(s$sentrix_id), ]
})
lapply(epic_samps, dim)
View(epic_samps$Armystarrs)



# get only 1 & 2 visit, removing visit 3
# epic_samps$Armystarrs <- epic_samps$Armystarrs[epic_samps$Armystarrs$visit !=3, ]
View(epic_samps$Armystarrs %>% select(contains(c("EWAS_id", "pcl"))))

# Map(function(x,y) table(x$EWAS_id %in% y$EWAS_id), pheno, pheno_comb)

# pheno_comb <- Map(function(x,y) merge(x, y, by="EWAS_id"), pheno_comb, pheno)

# epic_samps <- armystrs_pheno %>% filter(grepl("Epic", array, ignore.case = T))



# file with cell estimations
armystrs_cell_est <- "G:/PGC ML/ArmySTARRS/STARRS_Pheno_PCs_CellTypes_With_smoking_scores.csv"
prismo_cell_est <- "G:/PGC ML/PRISMO/Prismo_Pheno_051220_With_smoking_scores.csv"
cell_est_paths <- c(armystrs_cell_est, prismo_cell_est)
cell_ests <- lapply(cell_est_paths, read.csv)
names(cell_ests) <- df_names

table(epic_samps$Armystarrs$BaseName %in% cell_ests$Armystarrs$SID)
lapply(epic_samps, dim)


armystrs_pheno <- merge(epic_samps$Armystarrs, cell_ests$Armystarrs,
                   by.x = "BaseName",
                   by.y = "SID")

prismo_pheno <- merge(epic_samps$Prismo, cell_ests$Prismo,
                      by.x = "BaseName",
                      by.y = "ID")



# As some columns are common in both and while merging will add .x and .y to the columns.
# lets remove the columns with .y and then remove .x from the remaining columns
phenos <- list(armystrs_pheno, prismo_pheno)

lapply(phenos, function(x) table(duplicated(x$BaseName)))

# if there were same columns in two dfs, .x/.y is added on merge
# lets remove if there are any
phenos <- lapply(phenos, function(x) {
  x[, which(!grepl("\\.y$", colnames(x),
                                ignore.case = T))]
})
names(phenos) <- df_names

lapply(phenos, dim)

# Now remove the .x part from the column names
phenos <- lapply(phenos, function(x) {
  colnames(x) <- gsub(".x$", "", colnames(x))
  x
})

write.csv(phenos$Armystarrs, "G:/PGC ML/ArmySTARRS/armystarrs_Pheno_ML_updated.csv", row.names = F)

write.csv(phenos$Prismo, "G:/PGC ML/PRISMO//prismo_Pheno_ML_updated.csv", row.names = F)

# We will also convert CSV methylation file to feather file
# So that we can read it fast in python
armystrs_meth <- fread("G:/PGC ML/ArmySTARRS/Starrs_noob_qcd_crossReactiveProbesRemoved_combat_CP_wcovar_age2TP_ptsd_allPreAsControls.csv",
                  data.table = F)

write_feather(armystrs_meth, "G:/PGC ML/ArmySTARRS/Starrs_noob_qcd_crossReactiveProbesRemoved_combat_CP_wcovar_age2TP_ptsd_allPreAsControls.feather")


prismo_meth <- fread("G:/PGC ML/PRISMO/Prismo_noob_qcd_crossReactiveProbesRemoved_combat_CP_wcovar_age_ptsd_allPreAsControls.csv",
                     data.table = F)
write_feather(prismo_meth, "G:/PGC ML/PRISMO/Prismo_noob_qcd_crossReactiveProbesRemoved_combat_CP_wcovar_age_ptsd_allPreAsControls.feather" )


