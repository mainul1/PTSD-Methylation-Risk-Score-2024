# The code is to show how to perform covariate adjustment
# and calculate risk scores using data
# containing important features and the weights


# Step 1: ----------------------------------------
# Load the data. This should be the data from your cohort.
# For demonstration purpose a sample data is given
# containing important cpgs and covariates and phenotypes
# For your larger data set, you need to pull out the
# features given in the sample dataset

sample_df <- read.csv("G:/PGC ML/Combined Data/Sample data/Sample_data.csv") # path to the file
head(sample_df[, 1:5])

# Now name female:2, male:1 in the gender column
sample_df$Gender <- ifelse(sample_df$Gender == 'F', 2, 1)


# Step 2: ---------------------------------
# Covariate adjustment for methylation data
# use the listed covariates.

covars <- c("Bcell", "CD4T", "CD8T",
            "Mono", "NK", "SmoS", "Comp.2",
            "Comp.3", "Age")

covars_df <- sample_df[, covars]

# Get the beta values
beta <- sample_df %>% select(starts_with(c("cg", "ch.")))
dim(beta)

# convert to m-values
m_vals <- log(beta/(1-beta))
head(m_vals[, 1:5])


# Apply linear model to each cpg and
# replace the beta values with the residuals
adj_df <- lapply(seq_along(1:ncol(m_vals)), function(i){
  cpg <- colnames(m_vals)[i]
  f <- as.formula(paste(cpg,
                        paste(covars, collapse = "+"),
                        sep = "~"))

  res <- lm(formula = f, data = cbind(m_vals, covars_df))
  residual <- round(resid(res), 5) # upto 5 significant digits

})


# add colnames
names(adj_df) <- colnames(m_vals)

# make dataframe
residual_df <- data.frame(do.call(cbind, adj_df))
View(residual_df)

# convert to beta again
beta_adj <- 1/(1+(1/exp(residual_df)))
View(beta_adj)



# Step 3: ---------------------------------
# After covariate adjustment, let scale trauma and
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

scale_cols <- c("TraumaNumber", "ChildhoodMaltreatment")
scale_df <- sample_df[, which(names(sample_df) %in% (scale_cols))]
scale_df <- data.frame(maxmin(x = scale_df, na.rm = T))

# other_df <- test_sc[, -which(names(test_sc) %in% (scale_cols))]

# Now combine beta adjusted, gender and scaled columns
beta_adj_pheno <- cbind(beta_adj, scale_df)
beta_adj_pheno$Gender <- sample_df$Gender

# test_sc <- adj_df[, -which(names(adj_df) %in% (c("risk_score", "CurrentPTSD")))]
# dim(test_sc)



# Step 4: ------------------------------------
# Now to calculate risk scores load the weights
# given with the dataset

weights <- read.csv("G:/PGC ML/Combined Data/2021-11-27_19-16-35/Important_features.csv")
View(weights)

# change the name of two features to make it uniform
weights$Feature <- ifelse(weights$Feature == "Traumanum", 'TraumaNumber',
                          ifelse(weights$Feature == "Childhood_Mt",
                                 "ChildhoodMaltreatment", weights$Feature))

# save beta values and adjusted beta values
# pheno_cols <- c("Gender", "TraumaNumber",
#                 "ChildhoodMaltreatment")

ptsd_cols <- c("LifetimePTSD","CurrentPTSD")


# Now add more phenotypes
# beta_pheno <- cbind(beta, sample_df[, pheno_cols])
# beta_adj_pheno <- cbind(beta_adj, sample_df[, pheno_cols])

# Sort according to feature weights
# beta_pheno <- beta_pheno[, order(match(colnames(beta_pheno), weights$Feature))]

beta_adj_pheno <- beta_adj_pheno[, order(match(colnames(beta_adj_pheno),
                                               weights$Feature))]

# check order, all should be true
# table(colnames(beta_pheno) == weights$Feature)
table(colnames(beta_adj_pheno) == weights$Feature)


# Calculate weighted sum
# beta_pheno$risk_score <- apply(beta_pheno, 1, function(x) sum(x * weights$importance))
beta_adj_pheno$risk_score <- apply(beta_adj_pheno, 1, function(x) sum(x * weights$Importance))

# calculate mean methylation of unadjusted and adjusted data
# beta_pheno <- beta_pheno %>%
#   mutate(mean_meth =  rowMeans(select(., starts_with(c('cg','ch.')))))

beta_adj_pheno <- beta_adj_pheno %>%
  mutate(mean_meth = rowMeans(select(., starts_with(c('cg','ch.')))))


# finally add PTSD

# beta_pheno <- cbind(beta_pheno, sample_df[, ptsd_cols])
ptsd_cols <- c("LifetimePTSD","CurrentPTSD")
beta_adj_pheno <- cbind(beta_adj_pheno, sample_df[, ptsd_cols])


# save data
write.csv(beta_pheno, "G:/PGC ML/Combined Data/Sample data/Unadjusted data.csv",
          row.names = F)


write.csv(beta_adj_pheno, "G:/PGC ML/Combined Data/Sample data/Covariate adjusted data.csv",
          row.names = F)



