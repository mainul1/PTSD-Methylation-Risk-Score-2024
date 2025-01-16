

# This code is to create sample data for Kurt and Seyma
# to run the risk score analysis on an independent cohort.

library(openxlsx)

# data when using random forest to select the features
# fn <- ("G:/PGC ML/Combined Data/2022-03-30_15-07-11/Current_ptsd_important_features.xlsx")

# Features using elastic net
fn <- ("G:/PGC ML/Combined Data/2022-03-30_15-07-11/ElasticNet_Current_ptsd_important_features.xlsx")

# function to read all sheets
read_allsheets <- function(path){
  sheet_nms <- getSheetNames(path)
  sheets <- lapply(sheet_nms, read.xlsx, xlsxFile = path)
  names(sheets) <- sheet_nms
  return(sheets)
}

# read sheets
samp_data <- read_allsheets(path = fn)

# dim
lapply(samp_data, dim)
lapply(samp_data, function(x) head(x[, 1:5]))

# get some random 20 samples
set.seed(42)
samp_data_rnd <- lapply(samp_data, function(x) x[sample(nrow(x), 20), ])
View(samp_data_rnd$`Without NonCpGXY Probes`)

# how many random numbers for cpgs we need
# columns of df - 4 cols (id, childhood and others) * 20(samples)
noncpg <- (ncol(samp_data$`Without NonCpG Probes`)-4)*20
noncpgxy <- (ncol(samp_data$`Without NonCpGXY Probes`)-4)*20

rand_mat_noncpg <- matrix(runif(noncpg), nrow = 20) #
dim(rand_mat_noncpg)

rand_mat_noncpgxy <- matrix(runif(noncpgxy), nrow = 20)
dim(rand_mat_noncpgxy)

samp_data_rnd$`Without NonCpG Probes`[, 2:2779] <- rand_mat_noncpg
samp_data_rnd$`Without NonCpGXY Probes`[, 2:3729] <- rand_mat_noncpgxy

View(samp_data_rnd$`Without NonCpG Probes`)

View(samp_data_rnd$`Without NonCpGXY Probes`)

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

# add phenotype information to each df (noncpg and noncpgxy)
rand_data <- lapply(samp_data_rnd, function(x) cbind(x, pheno_sub))
lapply(rand_data, dim)

# lets remove some duplicate columns
dup_cols <- c("Childhood_Mt", "Traumanum", "current_ptsd", "SampleID")
rand_data_final <- lapply(rand_data, function(x){
  x %>%
    select(-c(dup_cols))
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


p <- "data/Sample_data.xlsx"

write.xlsx(rand_data_final,
           file = p)


