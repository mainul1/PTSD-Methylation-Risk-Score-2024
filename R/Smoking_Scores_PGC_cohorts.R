#########################################################################
### PGC Pipeline - Script 3.2 - Calculate smoking score
#########################################################################
#' remove all objects from your workspace
rm(list=ls())
gc()

#' Load all packages, if needed
library(data.table)
library(tibble)
library(tools)

#' This function will calculate the smoking score and combine it with the phenotype information
#' Higher scores are associated with smoking
#' input: Combat adjusted beta values, phenotype file, smoking coefficients,
#' Output: Phenotye information with smoking scores
#'
smoking_score <- function(beta, pheno, Scoefs, xcol){
  # Changing beta values of 0 to 0.0001
  beta[beta < 0.0001] <- 0.0001
  beta[beta > 0.9999] <- 0.9999

  # Convert to M-vals
  # --------------------------------
  # Do we need to convert beta to m values here
  # because they are not used after conversion
  # SK: They are used for matrix multiplication (only the selected probes)
  # SK: Technically we can convert to M-values after selecting the 39 smoking probes (Sprobes2), but I kept Mark's order here.
  beta <- log2(beta/(1-beta))

  # How many probes are we using?
  Sprobes <- beta[which(row.names(beta) %in% as.character(Scoefs$Marker)),]
  message("Number of probes used (out of 39): ", nrow(Sprobes))
  message(paste0("Missing probes: "),
          paste0(Scoefs$Marker[!Scoefs$Marker%in%row.names(Sprobes)], collapse = ","))
  Scoefs2 <- Scoefs[as.character(Scoefs$Marker) %in% row.names(Sprobes),]
  Sprobes2 <- Sprobes[as.character(Scoefs2$Marker), ]
  stopifnot(all(rownames(Sprobes2)==Scoefs2$Marker))

  ## Since EWAS was computed for M values, we need to transform.
  ## SK: This line is from Mark, I just copied and pasted, but I think he refers to previous M-value transformation
  ## SK: We should delete this line, it's confusing :)
  # ----------------------------------------- Not sure what you mean by the comment here
  Smo <- t(Scoefs2$Coefficient) %*% as.matrix(Sprobes2)
  Smo2 <- t(Smo)
  SmoScore <- data.frame(row.names(Smo2),Smo2)
  names(SmoScore) <- c("ID", "SmoS")

  # also removed all.x = True because by default merge will combine the ids that are common
  ## SK: The reason I keep all.x = T is to not lose phenotype information, and get NA for samples don't have methylation data
  ## SK: Because, I don't remove QC'ed out samples from my main phenotype document, just flag them; but remove them from methylation file.
  print(head(SmoScore))

  dat <- merge(pheno, SmoScore, by.x = xcol, by.y = "ID")
  return(dat)
}

## Load ComBat adjusted beta values
beta_file_paths <- c("G:/GTP Data/QCd Data/GTP_Noob_QCd_Combat_adj.feather",
                     "G:/PGC ML/MRS/MRS_noob_qcd_crossReactiveProbesRemoved_combat_CP_wcovar_age_ptsd_allPreAsControls.feather",
                     "G:/PGC ML/ArmySTARRS/Starrs_noob_qcd_crossReactiveProbesRemoved_combat_CP_wcovar_age2TP_ptsd_allPreAsControls.csv",
                     "G:/PGC ML/PRISMO/Prismo_noob_qcd_crossReactiveProbesRemoved_combat_CP_wcovar_age_ptsd_allPreAsControls.csv")

pheno_file_paths <- c("G:/GTP Data/QCd Data/Pheno_662_samps.csv",
                      "G:/PGC ML/MRS/MRS_Pheno.csv",
                      "G:/PGC ML/ArmySTARRS/STARRS_Pheno_PCs_CellTypes.csv",
                      "G:/PGC ML/PRISMO/Prismo_Pheno_051220.csv")

df_names <- c("GTP", "MRS", "ArmyStarrs", "Prismo")


# read beta files
meth_data <- lapply(beta_file_paths, function(x){
  message("Reading file: ", basename(x))
  if(file_ext(x) == 'feather'){
    d <- read_feather(x)
  }else{
    d <- fread(x, data.table = F)
  }
  d <- column_to_rownames(d ,var = colnames(d)[1])
})

names(meth_data) <- df_names

lapply(meth_data, function(x) head(x[, 1:5]))
lapply(meth_data, dim)


## Load the phenotype file
pheno_data <- lapply(pheno_file_paths, read.csv)
names(pheno_data) <- df_names
lapply(pheno_data, function(x) head(x[, 1:5]))
lapply(pheno_data, dim)


#Now get the samples that are in methylation data
pheno_data <- Map(function(x, y) x[which(x[, 1] %in% colnames(y)), ], pheno_data, meth_data)

# Lets check if the methylation and phenotype data is in order
Map(function(x,y) table(colnames(x) == y[, 1]), meth_data, pheno_data)

pheno_ordered <- Map(function(x,y) x[order(match(x[, 1], colnames(y))), ], pheno_data, meth_data)
Map(function(x,y) table(colnames(x) == y[, 1]), meth_data, pheno_ordered)


## Load smoking score coefficients shared in the pipeline package
Scoefs <- read.csv("G:/DNHS_2nd_Batch(using PGC pipeline)/Smoking_probes_betas.csv",stringsAsFactors=FALSE)
names(Scoefs)<-c("Marker","Coefficient")


## Run the function to compute smoking score
scores <- Map(function(x,y) smoking_score(beta = x, pheno = y,
                                          Scoefs = Scoefs,
                                          xcol = colnames(y)[1]),
              meth_data, pheno_ordered)

lapply(scores, dim)

lapply(scores, function(x) head(x[, 1:5]))

Map(function(x,y) ncol(x) == nrow(y), meth_data, scores) # check order

# File names to save the phenotype data with smoking scores
out_file_names <- paste0(gsub(".csv", "", pheno_file_paths), "_With_smoking_scores.csv")

# write data
Map(function(x,y) write.csv(x, y, row.names = F), scores, out_file_names)

