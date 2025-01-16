
# ------------------ Run this file first ----------------

# The code is to show how to perform covariate adjustment
# and calculate risk scores using data
# containing important features and the weights

# Based on the feedback from the group we tested models on two types of datasets
# 1. Removing all non-CpG probes
# 2. Removing all non-CpG and X or Y chromosome probes together


# Step 1: ----------------------------------------
# Load the data containing important features.
# This should be the data from your cohort processed using PGC pipeline.
# For demonstration purpose a sample data is given
# containing important cpgs and covariates and phenotypes
# For your data set, you need to pull out the
# features (columns) given in the sample dataset ('Sample_data.xlsx')

# If some important features are missing from your data, drop those features
#

# load libraries
library(openxlsx)
library(dplyr)


# Change path as per the location on your computer
# As an example sample data is provided
# But this data needs to be from your cohort containing
# important features and PTSD variables
# given in the "EN_selected_wd_EN_l1_r_0.1.xlsx"

fname <- "data/Sample_data.xlsx" # path to the file


# function to read all sheets
# Input is file path as given above
read_allsheets <- function(path){
  sheet_nms <- getSheetNames(path)
  sheets <- lapply(sheet_nms, read.xlsx, xlsxFile = path)
  names(sheets) <- sheet_nms
  return(sheets)
}

# function call to read data from sample sheets
sample_df <- read_allsheets(path = fname)


# display dimension and head of each df
lapply(sample_df, dim)
lapply(sample_df, function(x) head(x[, 1:5]))


# Now name female:2, male:1 in the gender column
sample_df <- lapply(sample_df, function(x){
  x$Gender <- ifelse(x$Gender == 'F', 2, 1) # replace the name according to your data
  x
})


# Step 2: ---------------------------------
# Covariate adjustment for methylation data
# use the listed covariates to adjust data.
# add "Gender" as well

covars <- c("Bcell", "CD4T", "CD8T",
            "Mono", "NK", "SmoS", "Comp.2",
            "Comp.3", "Age", "Gender")

# pull required columns
covars_df <- lapply(sample_df, function(x) x[, covars])

# Get the beta values only
beta <- lapply(sample_df, function(x) x %>% select(starts_with(c("cg", "ch."))))
lapply(beta, dim)

# convert to m-values
m_vals <- lapply(beta, function(x) log(x/(1-x)))
lapply(m_vals, function(x) head(x[, 1:5]))


# Apply linear model to each cpg (cpg as dependent variable) and
# replace the beta values with the residuals
adj_df <- lapply(seq_along(m_vals), function(i) { # loop over dfs
  x <- lapply(seq_along(1:ncol(m_vals[[i]])), function(j){ # loop over each column (cpg)
    cpg <- colnames(m_vals[[i]])[j]
    f <- as.formula(paste(cpg,
                          paste(covars, collapse = "+"),
                          sep = "~"))

    res <- lm(formula = f, data = cbind(m_vals[[i]], covars_df[[i]]))
    residual <- round(resid(res), 5) # upto 5 significant digits

  })
})

names(adj_df) <- names(m_vals)


# add column names (cpg names) back to data
adj_df <- lapply(seq_along(m_vals), function(i){
  names(adj_df[[i]]) <- colnames(m_vals[[i]])
  adj_df[[i]]
})
names(adj_df) <- names(m_vals) # name to dfs


# make dataframe
residual_df <- lapply(adj_df, function(x) data.frame(do.call(cbind, x)))
View(residual_df$`Without NonCpG Probes`)
View(residual_df$`Without NonCpGXY Probes`)

# Let's convert to beta again
beta_adj <- lapply(residual_df, function(x) 1/(1+(1/exp(x))))
View(beta_adj$`Without NonCpG Probes`)



# Step 3: ---------------------------------
# After covariate adjustment, let's scale trauma and
# childhood maltreatment columns using min-max approach

# function to scale the data between 0-1
maxmin <- function(x, na.rm=TRUE){
  if(is.vector(x)==TRUE){
    maxs <- max(x, na.rm = na.rm)
    mins <- min(x, na.rm = na.rm)
    scale(x,center=mins,scale=maxs-mins)
  } else {
    maxs <- apply(x, 2,max)
    mins <- apply(x, 2,min)
    scale(x, center = mins, scale = maxs - mins)
  }
}


# We will scale trauma and childhoodmaltreatment columns
# and combine beta adjusted, gender and scaled columns

scale_cols <- c("TraumaNumber", "ChildhoodMaltreatment")
beta_adj_pheno <- lapply(seq_along(sample_df), function(i){
  scale_df <- sample_df[[i]][, which(names(sample_df[[i]]) %in% (scale_cols))]
  scale_df <- data.frame(maxmin(x = scale_df, na.rm = T))
  cbind(beta_adj[[i]], scale_df)
})
names(beta_adj_pheno) <- names(beta_adj)



# beta_adj_pheno$Gender <- sample_df$Gender



# Step 4: ------------------------------------
# Now to calculate risk scores load the weights
# given with the dataset

weights_fname <- "data/EN_selected_wd_EN_l1_r_0.1.xlsx"
weights <- read_allsheets(path = weights_fname)
names(weights) <- names(beta_adj_pheno)
View(weights$`Without NonCpG Probes`)

# change the name of two features to make it uniform
weights <- lapply(weights, function(x){
  x$Feature <- ifelse(x$Feature == "Traumanum", 'TraumaNumber',
         ifelse(x$Feature == "Childhood_Mt",
                "ChildhoodMaltreatment", x$Feature))
  x
})


View(weights$`Without NonCpGXY Probes`)


# Order the columns in adjusted data according to features
beta_adj_pheno <- lapply(seq_along(beta_adj_pheno), function(i){
  beta_adj_pheno[[i]][, order(match(colnames(beta_adj_pheno[[i]]),
                               weights[[i]]$Feature))]
})
names(beta_adj_pheno) <- names(beta_adj)

# check order, all should be true
table(colnames(beta_adj_pheno$`Without NonCpG Probes`) == weights$`Without NonCpG Probes`$Feature)
table(colnames(beta_adj_pheno$`Without NonCpGXY Probes`) == weights$`Without NonCpGXY Probes`$Feature)

# stopifnot(all(colnames(beta_adj_pheno) == weights$Feature))


# Risk scores are the weighted sum for each sample
# using feature importance as weights
beta_adj_pheno <- lapply(seq_along(beta_adj_pheno), function(i){

  # Methylation and exposure risk scores
  print(paste("Ordered: ",
              table(colnames(beta_adj_pheno[[i]]) ==
                      weights[[i]]$Feature),collapse = '\n'))
  beta_adj_pheno[[i]]$risk_score <- apply(beta_adj_pheno[[i]], 1,
                                          function(x) sum(x * weights[[i]]$Importance))

  # Methylation only risk scores
  meth_fea_weights <- weights[[i]] %>%
    filter(grepl("^cg|^ch.", Feature)) # get cpg weights only

  df_x <- beta_adj_pheno[[i]]
  meth <- df_x[, which(colnames(df_x) %in% meth_fea_weights$Feature)] # get cpgs

  ord_meth <- meth[, order(match(colnames(meth), meth_fea_weights$Feature))] # order

  print(paste("Ordered: ",
              table(colnames(ord_meth) ==
                      meth_fea_weights$Feature),collapse = '\n'))

  # calculate risk scores
  beta_adj_pheno[[i]]$meth_only_risk_score <- apply(ord_meth, 1,
                                                    function(x) sum(x * meth_fea_weights$Importance))

  beta_adj_pheno[[i]] # return
})



# Now calculate the mean methylation of the CpGs only
beta_adj_pheno <- lapply(beta_adj_pheno, function(x){
  x %>%
    mutate(mean_meth = rowMeans(select(., starts_with(c('cg','ch.')))))
})



# Finally add PTSD information
# If you dont have lifetime ptsd, remove it from the list
ptsd_cols <- c("Gender", "LifetimePTSD","CurrentPTSD")
beta_adj_pheno <- lapply(seq_along(beta_adj_pheno), function(i){
  cbind(beta_adj_pheno[[i]], sample_df[[i]][, ptsd_cols])
})
names(beta_adj_pheno) <- names(beta_adj)


# save data and use in the next script
write.xlsx(beta_adj_pheno,
           file = "data/Covariate adjusted and scaled data.xlsx",
          rowNames = F)



