# Lets select other phenotypes we need from our data

# 1. Symptom clusters

dnhs_pheno <- read.csv("E:/DNHS_EWAS_DATA/DNHSEWAS492_RProject/Data/Long_Stacked Combined DF without Imputation_Including_Smoking.csv")

cluster_cols <- c("w1c1_life_worst_intrusion", "w2c2_life_worst_intrusion",
                  "w1c1_life_worst_avoidance", "w2c2_life_worst_avoidance",
                  "w1c1_life_worst_hyperarousal", "w2c2_life_worst_hyperarousal")


# we have some variables in the other file
master_pheno <- read.csv("E:/DNHS_EWAS/DN170703_Uddin_copy.csv")

# MDD
phq9_cols <- c("w1c1_phq9sum", "w2c2_phq9sum")
gad_cols <- c("w1c1_gad7sum", "w1c1_gad7cat")

master_sub <- master_pheno[, c("RESP", cluster_cols, phq9_cols, gad_cols)]


# we have w2c2 gad7sum in a different file
# so lets get that and combine it with main file
w2c2 <- read.csv("E:/DNHS SAS Datasets/w2_survey_cohort2.csv")
w2c2_cols <- c("RESP", "gad7sum", "gad7cat")

w2c2_sub <- w2c2[, w2c2_cols]

comb <- merge(master_pheno, w2c2_sub, by= "RESP", all = T)

names(comb)[names(comb) == "gad7sum"] <- "w2c2_gad7sum"
names(comb)[names(comb) == "gad7cat"] <- "w2c2_gad7cat"

write.csv(comb, "E:/DNHS_EWAS/DN170703_Uddin_update_gad.csv")

comb_gad <- comb[, c("RESP", "w2c2_gad7sum", "w2c2_gad7cat")]


comb_all <- merge(master_sub, comb_gad, by = "RESP")


# Now combine the values from columns from w1c1 and w2c2 into one
pattern <- c("life_worst_intrusion", "life_worst_avoidance",
             "life_worst_hyperarousal",  "phq9sum",
             "gad7sum", "gad7cat")

for (i in seq_along(pattern)) {
  print(i)
  print(pattern[i])
  cols <- colnames(comb_all)[grepl(pattern[i], colnames(comb_all))]
  print(cols)
  df <- comb_all[, cols]
  print(head(df))
  comb_all[[pattern[i]]] <- ifelse(apply(is.na(df),1,all),NA,rowSums(df,na.rm=TRUE))
}



# Now we need to combine this information with the main file to
# store if for future use
table(master_pheno$RESP %in% comb_all$RESP)
table(dnhs_pheno$RESP %in% comb_all$RESP)

get_cols <- c("RESP", "life_worst_intrusion", "life_worst_avoidance",
              "life_worst_hyperarousal", "phq9sum", "gad7sum",
              "gad7cat", "w2c2_gad7cat")

final <- merge(dnhs_pheno, comb_all[, get_cols], by = "RESP")

write.csv(final, "E:/DNHS_EWAS_DATA/DNHSEWAS492_RProject/Data/Long_Stacked Combined DF without Imputation_Including_Smoking.csv")


# Now add the columns to Batch1 and 2 epic as well

pheno_b1_b2 <- read.csv("G:/DNHS 2nd Batch/DNHS2ndBatachAnalysis/data/DNHS_Pheno_Batch1&2_with_childhoodT.csv")
cols <- get_cols[-c(7,8)]
b1_b2_final <- merge(pheno_b1_b2, comb_all[, cols], by="RESP")

write.csv(b1_b2_final[, -2], "G:/DNHS 2nd Batch/DNHS2ndBatachAnalysis/data/DNHS_Pheno_Batch1&2_with_childhoodT_gad_phq.csv")
