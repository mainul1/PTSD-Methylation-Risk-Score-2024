# We will check if there is an overlap between MRS CpGs
# and 12 meta-analysis CpGs that were submitted with EWAS grant report

library(openxlsx)
library(dplyr)
library(data.table)

# Use the cpgs from models 1, 2 and 3. Models 1&2 have same CpGs but different weights

# 12 meta-analysis cpgs, reported in EWAS report
cpgs_12 <- read.xlsx("D:/HT12_Hits_Metaanalysis/Data/12CpGs_in_grant.xlsx")

# All in meta-analysis
all_cpgs <- read.xlsx("D:/HT12_Hits_Metaanalysis/Data/Copy of CpGlist.xlsx")

# 11 cpgs in Alicia's SOBP slides
cpgs_11 <- read.xlsx("D:/HT12_Hits_Metaanalysis/Data/11CpGs_in_SOBP_slides.xlsx")

meta_analysis <- list("12CpGs" = cpgs_12,
                      "all_CpGs" = all_cpgs,
                      "11CpGs" = cpgs_11)

lapply(meta_analysis, dim)


# ----------------------
# function to read all sheets
read_allsheets <- function(path){
  sheet_nms <- getSheetNames(path)
  sheets <- lapply(sheet_nms, read.xlsx, xlsxFile = path)
  names(sheets) <- sheet_nms
  return(sheets)
}

# MRS Cpgs
fname <- "G:/PGC ML/Combined Data/2022-03-30_15-07-11/Elasticnet risk scores ptsdpm test data.xlsx"
sheet_names <- getSheetNames(fname)
mrs_cpgs <- read_allsheets(path = fname)
lapply(mrs_cpgs, dim)
View(mrs_cpgs$`Without NonCpGXY Probes`)

# check how many are overlapping
lapply(meta_analysis, function(x){
  table(x$IlmnID %in% colnames(mrs_cpgs$`Without NonCpGXY Probes`))
})

# get those that are overlapping
overlap_wd_unadj <- lapply(meta_analysis, function(x){
  x %>%
    filter(IlmnID %in% colnames(mrs_cpgs$`Without NonCpGXY Probes`))
})

overlap_wd_unadj


# ------------------------------------------------
# MRS cpgs with adjusted exposure variables
fname_adj_exp <- "G:/PGC ML/Combined Data/2023-03-08_21-41-11/Elasticnet risk scores ptsdpm test data.xlsx"
mrs_cpgs_adj_exp <- read_allsheets(path = fname_adj_exp)
lapply(mrs_cpgs_adj_exp, dim)
names(mrs_cpgs_adj_exp)

# check how many are overlapping
lapply(meta_analysis, function(x){
  table(x$IlmnID %in% colnames(mrs_cpgs_adj_exp$`Without NonCpGXY Probes`))
})

# get those that are overlapping
overlap_wd_adj <- lapply(meta_analysis, function(x){
  x %>%
    filter(IlmnID %in% colnames(mrs_cpgs_adj_exp$`Without NonCpGXY Probes`))
})

overlap_wd_adj


# Check how many are matching between two, exposure unadjusted and adjusted
Map(function(x, y) table(x$IlmnID %in% y$IlmnID), overlap_wd_unadj, overlap_wd_adj)

final <- list("12CpGs_MRSCpGs" = overlap_wd_unadj$`12CpGs`,
              "allCpGs_MRSCpGs" = overlap_wd_unadj$all_CpGs,
              "12CpGs_Adj_Exp_MRSCpGs" = overlap_wd_adj$`12CpGs`,
              "allCpGs_Adj_Exp_MRSCpGs" = overlap_wd_adj$all_CpGs,
              "11CpGs_MRSCpGs" = overlap_wd_unadj$`11CpGs`,
              "11CpGs_Adj_Exp_MRSCpGs" = overlap_wd_adj$`11CpGs`)


# both  datasets, without adjusted exposures and with adjusted exposures
# have same overlap (3 Cpgs)
# so lets save
write.xlsx(final, "D:/HT12_Hits_Metaanalysis/Data/Overlap_between_MRS_CpGs_and_metaanalysis_CpGs.xlsx")



# Lets look at the overlap between Cpgs(gene wise) with genes identified in Freeze3 PGC PTSD data

# We need to annotate cpgs so lets load the annotation file
annot_file <- fread("D:/EWAS meta-analysis(Janelle)/data/infinium-methylationepic-v-1-0-b5-manifest-file.csv",
                    data.table = F, skip = 7, fill = T)

dim(annot_file)

annot_file_sub <- annot_file %>% dplyr::select(IlmnID, Infinium_Design_Type,
                                               Genome_Build, CHR, MAPINFO,
                                               UCSC_RefGene_Name, UCSC_RefGene_Group,
                                               Relation_to_UCSC_CpG_Island)
head(annot_file_sub)


# get the genes of Cpgs from model 1 and 2

genes_m1m2 <- annot_file_sub %>%
  filter(IlmnID %in% colnames(mrs_cpgs$`Without NonCpGXY Probes`))
dim(genes_m1m2)
table(genes_m1m2$IlmnID %in% colnames(mrs_cpgs$`Without NonCpGXY Probes`))


# model 3
genes_m3 <- annot_file_sub %>%
  filter(IlmnID %in% colnames(mrs_cpgs_adj_exp$`Without NonCpGXY Probes`))
dim(genes_m3)

# load Freeze 3 data
f3 <- read.xlsx("../../../Manuscript/Previous published/Freeze_3_Manuscript_Apri25_2023/Freeze_3_Manuscript_Apri25_2023/SupplementaryTables.xlsx",
                sheet = "9. Gene-mapping", startRow = 3, colNames = TRUE)
head(f3)

# overlap betweem genes from models 1 and 2 and Freeze 3
get_overlap <- function(model_ls, comp_gene_ls){
  gene_ls <- strsplit(model_ls$UCSC_RefGene_Name, ';')
  unq_genes <- unique(unlist(gene_ls))
  overlap <- comp_gene_ls[comp_gene_ls$symbol %in% unq_genes, ]
  message("overlaping genes :")
  print(dim(overlap))
  print(table(comp_gene_ls$symbol %in% unq_genes))
  overlap
}

models_1and2_overl <- get_overlap(genes_m1m2, f3)

# overlap between genes from model 3 and Freeze 3
model3_overl <- get_overlap(genes_m3, f3)


overlap_wd_f3 <- list("m1m2" = models_1and2_overl,
                      "m3" = model3_overl)

write.xlsx(overlap_wd_f3, file = "data/Overlaping_genes_from_m1_2&3_with_freeze3.xlsx")


# Overlap wih Dean et al and Schulterbrauks et al papers
dean_cpgs <- read.xlsx("../../Manuscript/Previous published/BiomarkerPanel_343 (1).xlsx")
head(dean_cpgs)

table(grepl('^cg', dean_cpgs$ProbeID))

dean_cpgs <- dean_cpgs[grepl('^cg', dean_cpgs$ProbeID), ]
dean_cpgs$Description <- gsub("Gene:", '', dean_cpgs$Description)
genes <- dean_cpgs$Description
genes <- genes[!is.na(genes)]

length(unique(genes))

table(colnames(mrs_cpgs$`Without NonCpGXY Probes`) %in% dean_cpgs$ProbeID)
table(colnames(mrs_cpgs_adj_exp$`Without NonCpGXY Probes`) %in% dean_cpgs$ProbeID)

dean_m1m2 <- dean_cpgs[dean_cpgs$ProbeID %in% colnames(mrs_cpgs$`Without NonCpGXY Probes`), ]
dean_m1m2

dean_m3 <- dean_cpgs[dean_cpgs$ProbeID %in% colnames(mrs_cpgs_adj_exp$`Without NonCpGXY Probes`), ]
dean_m3


schult_cpgs <- read.xlsx("../../../Manuscript/Previous published/Hits from Katharina Schultebraucks 2020 .xlsx")
head(schult_cpgs)
table(colnames(mrs_cpgs$`Without NonCpGXY Probes`) %in% schult_cpgs$CpGs)
table(colnames(mrs_cpgs_adj_exp$`Without NonCpGXY Probes`) %in% schult_cpgs$CpGs)


dean <- list("dean_wd_m1m1Cpgs" = dean_m1m2,
             "dean_wd_m3Cpgs" = dean_m3)

write.xlsx(dean, "../data/Overlap_of Cpgs_m1m2m3_wd_dean_paper.xlsx")
