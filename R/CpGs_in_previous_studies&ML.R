# code to investigate previously identified CpGs with 
# the set found using PGCML
library(xlsx)
library(stringr)

cpgs <- read.xlsx("../data/CpGs previosly found PGC.xlsx", 
                  sheetIndex = 1)
head(cpgs)
dim(cpgs)


ml_cpgs <- read.csv("G:/PGC ML/Combined Data/2021-11-27_19-16-35/Important_features.csv")
head(ml_cpgs)


# How many are matching 
cpgs[cpgs$CPG %in% ml_cpgs$Feature, ]
table(cpgs$CPG %in% ml_cpgs$Feature)

gene_ls <- str_trim(cpgs$Gene)
gene_ls <- paste0("^",gene_ls)
gene_ls

# check gene-wise
annotated_cpgs <- read.csv("G:/PGC ML/Combined Data/2021-11-27_19-16-35/Annotated_important_cpgs.csv")
dim(annotated_cpgs)

genes_found <- annotated_cpgs[which(grepl(paste(c(gene_ls), collapse = '|'), annotated_cpgs$UCSC_RefGene_Name)), ]
genes_found


