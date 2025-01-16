

# This code is to create sample data for Kurt and Seyma
# to run the risk score analysis on an independent cohort.
# this dataset is after updating the code with exposure variables
# in the model and then adjusted for exposure variables


library(openxlsx)
library(dplyr)

# data when using random forest to select the features
# fn <- ("G:/PGC ML/Combined Data/2022-03-30_15-07-11/Current_ptsd_important_features.xlsx")

# Features using elastic net. When we remove exposure variables in the model
# the CpGs are same in the model but the weights may change
# so load the features we have from the model and use updated weights
fn <- ("G:/PGC ML/Combined Data/2022-03-30_15-07-11/ElasticNet_Current_ptsd_important_features.xlsx")

# with adjusted exposure variables
fn_adj_exps <- "G:/PGC ML/Combined Data/2023-03-08_21-41-11/ElasticNet_Current_ptsd_important_features_Adj_wd_Exp_Vars.xlsx"

# function to read all sheets
read_allsheets <- function(path){
  sheet_nms <- getSheetNames(path)
  sheets <- lapply(sheet_nms, read.xlsx, xlsxFile = path)
  names(sheets) <- sheet_nms
  return(sheets)
}

# read sheets
samp_data_old <- read_allsheets(path = fn)

samp_data_wd_adj_exps <- read_allsheets(path = fn_adj_exps)


# we'll get only nonCpGs and XY probes removed data
samp_data <- list("Without Exposures in model" = samp_data_old$`Without NonCpGXY Probes`,
                  "With Adjusted Exposures" = samp_data_wd_adj_exps$`Without NonCpGXY Probes`)

# dim
lapply(samp_data, dim)
lapply(samp_data, function(x) head(x[, 1:5]))

# get some random 20 samples
set.seed(42)
samp_data_rnd <- lapply(samp_data, function(x) x[sample(nrow(x), 20), ])
View(samp_data_rnd$`Without Exposures in model`)

# how many random numbers for cpgs we need
# columns of df - 4 cols (id, childhood and others) * 20(samples)
wo_exp <- (ncol(samp_data$`Without Exposures in model`)-4)*20
adj_exp <- (ncol(samp_data$`With Adjusted Exposures`)-2)*20 # we have only 2 cols extra

rand_mat_wo_exp <- matrix(runif(wo_exp), nrow = 20) #
dim(rand_mat_wo_exp)

rand_mat_adj_exp <- matrix(runif(adj_exp), nrow = 20)
dim(rand_mat_adj_exp)

# insert random data for CpGs
samp_data_rnd$`Without Exposures in model`[, 2:3729] <- rand_mat_wo_exp
samp_data_rnd$`With Adjusted Exposures`[, 2:4151] <- rand_mat_adj_exp

View(samp_data_rnd$`Without Exposures in model`)

View(samp_data_rnd$`With Adjusted Exposures`)

lapply(samp_data_rnd, dim)

# Now load phenotype data to add phenotype info
pheno <- read.csv("E:/DNHS_EWAS_DATA/DNHSEWAS492_RProject/Data For Traditional analysis/DNHS_Pheno_wd_PCs_outlierInfo.csv")
head(pheno)

pheno2 <- read.csv("G:/DNHS 2nd Batch/DNHS2ndBatachAnalysis/data/pheno_PCs_QC_with_smoking_scores.csv")
View(pheno2)

# get smoking info
smok_df <- pheno2[, c("X", "SmoS")]

# merge pheno with smoking scores
pheno <-  merge(pheno, smok_df,
                by.x = "SampleID",
                by.y = "X")


# remove remitted
pheno_clean <- pheno[!(pheno$PTSDLife == 1 & pheno$PTSDpm == 0), ]

pheno_clean <- pheno_clean %>% filter(TraumaNum > 0)

unq_ppts <- pheno_clean[!duplicated(pheno_clean$RESP), ]
dim(unq_ppts)

samples <- unq_ppts %>% group_by(PTSDpm) %>% slice(1:10) #20 = 10 from cases, 10 from controls
dim(samples)


# now add cumulative trauma, Gender, lifetime/current ptsd, childhood maltreatment
cols <- c("SampleID", "Gender", "Age", "PTSDLife", "PTSDpm",
          "TraumaNum", "CD8T", "CD4T","NK", "Bcell",
          "Mono", "Comp.2", "Comp.3", "SmoS",
          "childhood_cum_trauma")

# to make names more readable, assign new names
new_names <- c("SampleID", "Gender", "Age", "LifetimePTSD", "CurrentPTSD",
               "TraumaNumber", "CD8T", "CD4T","NK", "Bcell",
               "Mono", "Comp.2", "Comp.3", "SmoS",
               "ChildhoodMaltreatment")

pheno_sub <- samples[, cols]

colnames(pheno_sub) <- new_names

# add phenotype information to each df (adj_exp and noncpgxy)
rand_data <- lapply(samp_data_rnd, function(x) cbind(x, pheno_sub))
lapply(rand_data, dim)

# lets remove some duplicate columns
dup_cols <- c("Childhood_Mt", "Traumanum", "current_ptsd", "SampleID")
rand_data_final <- lapply(rand_data, function(x){
  x %>%
    select(-c(colnames(.)[colnames(.) %in% dup_cols]))
})

lapply(rand_data_final, dim)

# check if we have all important cpgs
cpgs <- lapply(samp_data_rnd, function(x){
  x %>%
    select(starts_with('cg'))
})
lapply(cpgs, dim)

# all should be true
Map(function(x, y) table(colnames(x) %in% colnames(y)), cpgs, rand_data_final)


p <- "data/Sample_data_updated_for_adj_exp.xlsx"

write.xlsx(rand_data_final,
           file = p)


# lets also combine the weights of without exposures in the model and adjusted exposures
weights_wo_exp_inmodel <- read.csv("G:/PGC ML/Combined Data/2022-03-30_15-07-11/Important_features_wo_non_CpGsXY_EN_WO_Exposure_Vars_selected_wd_EN_l1_r_0.1.csv")
dim(weights_wo_exp_inmodel)

weights_wd_adj_exps <- read.csv("G:/PGC ML/Combined Data/2023-03-08_21-41-11/Important_features_wo_non_CpGsXY_EN_selected_wd_EN_l1_r_0.1.csv")
dim(weights_wd_adj_exps)


weights_ls <- list("Without Exposures in model" = weights_wo_exp_inmodel,
     "With Adjusted Exposures" = weights_wd_adj_exps)

lapply(weights_ls, dim)

# save
write.xlsx(weights_ls, file = "Data/Wo_exp_and_adj_exps_selected_wd_EN_l1_r_0.1.xlsx")
