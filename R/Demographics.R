
# This code is to get demographics of the samples used in ML

library(dplyr)
library(tidyr)
library(kableExtra)
library(xlsx)

# load all data
all_cohorts <- read.csv("G:/PGC ML/Pre_Processed Data/2021-11-15_21-41-53/DNHS_GTP_MRS_ArmyS_Prismo_Pheno.csv")

# remove some columns that we alread have
all_cohorts <- all_cohorts %>% select(-c("Childhood_Mt", "Gender", "Traumanum"))
View(all_cohorts)


# Data used in training and testing
path <- "G:/PGC ML/Combined Data/2021-11-27_19-16-35/"
current <- read.csv(paste0(path, "Current_ptsd_important_features.csv"))

lifetime <- read.csv(paste0(path, "Lifetime_ptsd_important_features.csv"))

table(colnames(current) == colnames(lifetime))


# Now combine current and lifetime to get phenotype information
current_comb <- merge(current, all_cohorts,
                      by = "Basename")
dim(current_comb)


lifetime_comb <- merge(lifetime, all_cohorts,
                       by = "Basename")
dim(lifetime_comb)


# Get the number of cases and controls in each
table(current_comb$current_ptsd)
table(current_comb$Gender)

table(lifetime_comb$lifetime_ptsd)
table(lifetime_comb$Gender)

# function to get cases and control count for individual cohorts
get_count <- function(df, group, col, wide = FALSE,
                      all= TRUE){
  ptsd <- df %>% group_by(across(all_of(group))) %>%
    count(across(all_of(col))) %>%
    mutate(sum = sum(n)) %>%
    ungroup() %>%
    mutate(total = sum(n))
  print(ptsd)

  # information of all discovery cohorts together
  all_info <- df %>%
    count(across(all_of(col))) %>%
    mutate(sum = sum(n),
           total = sum(n)) %>%
    mutate(Study ='All cohorts ', .before = 1)
    print(all_info)

  ptsd <- rbind(ptsd, all_info)

  # make resulting df wide
  if(wide){
    ptsd <- ptsd %>% pivot_wider(names_from = col,
                                 values_from = n,
                                 names_sort = TRUE)
    colnames(ptsd) <- c("Study", paste(col, colnames(ptsd)[-1], sep = "_"))
  }
  ptsd

}


# call function for current and lifetime ptsd

cur_ptsd <- get_count(df = current_comb, group = "Study",
          col = "current_ptsd", wide = TRUE)
cur_ptsd <- cur_ptsd[, c(1,5,4,2,3)]
cur_ptsd

life_ptsd <- get_count(df = lifetime_comb, group = "Study",
                       col = "lifetime_ptsd", wide = TRUE)
life_ptsd <- life_ptsd[, c(1,5,4,2,3)]
life_ptsd

total <- merge(cur_ptsd, life_ptsd,
               by = 1, all.x = T)
total

# added "_updated" to file after add trauma variable
# f_path <- "data/Demographics_updated.xlsx"
f_path <- "../data/Demographics_updated_5_26_23.xlsx" # added all discovery cohorts combined info
write.xlsx(total, file = f_path, sheetName = "Number")


# Function to get mean and sd
get_mean_sd <- function(df, group, col, pval = FALSE,
                        wide = FALSE, outcome = ""){
  if(pval){
    mean_sd_func <- function(){
      df_grpd %>%
        summarize(mean_sd = paste0(round(across(all_of(col),
                                              mean, na.rm = T), 1),
                                 " (",
                                 round(across(all_of(col),
                                              sd, na.rm = T), 1), ")" ),
                t.test_p = if(length(unique(get(col))) != 1){
                  t.test(get(col) ~ get(outcome))$p.value} else NA,
                mann_whitney_p = if(length(unique(get(col))) != 1){
                  wilcox.test(get(col) ~ get(outcome))$p.value} else NA)
    }

  }else{

    mean_sd_func <- function(){
      df_grpd %>%
      summarize(mean_sd = paste0(round(across(all_of(col),
                                              mean, na.rm = T), 1),
                                     " (",
                                     round(across(all_of(col),
                                                  sd, na.rm = T), 1), ")" ))
    }
  }


  # all cohorts together
  if(!pval){
    df_grpd <- df %>% group_by(across(all_of(contains(outcome))))
  }else{df_grpd <- df}

  all_info <- mean_sd_func() %>%
    mutate(Study = "All cohorts", .before=1)
  # print(all_info)


  # individual cohorts
  df_grpd <- df %>%
    group_by(across(all_of(contains(group))))
  indv_info <- mean_sd_func()

  df <- rbind(indv_info, all_info)

  if(pval){
    df <- df %>% rename(!!paste0(outcome,"_mean_sd") := mean_sd,
                        !!paste0(outcome, "_t.test_p") := t.test_p,
                        !!paste0(outcome, "_mann_whitney_p") := mann_whitney_p)
  }


  if(wide){
    df <- df %>%
      pivot_wider(names_from = outcome,
                  names_prefix = outcome,
                  values_from = mean_sd)
  }
  df

}

# test_df <- current_life$Ptsdpm[current_life$Ptsdpm$Study == 'MRS', ]
# m1 <- glm(Traumanum ~ Ptsdpm, family="poisson", data=test_df)
# summary(m1)

# Age ---------------------------------------

current_life <- list(current_comb, lifetime_comb)
names(current_life) <- c("Ptsdpm", "Ptsdlife")
vars <- c("Age", "Pts_Severity", "Smos",
          "Childhood_Mt", "Traumanum")

# Get mean and sd for groups for cases and controls separately
sep <- lapply(seq_along(current_life), function(i) {
  m <- lapply(seq_along(vars), function(j) {
    outcome <- unlist(names(current_life)[i])
    group <- c("Study", outcome)
    get_mean_sd(df = current_life[[i]],
                group = group,
                col = vars[j],
                wide = TRUE,
                outcome = outcome)
  })
  names(m) <- vars
  m
    })

sep

# Mean and sd together, together for cases and controls with pvalue
tog <- lapply(seq_along(current_life), function(i) {
  m <- lapply(seq_along(vars), function(j) {
    outcome <- unlist(names(current_life)[i])
    get_mean_sd(df = current_life[[i]],
                group = "Study",
                col = vars[j],
                pval = TRUE,
                outcome = outcome)
  })
  names(m) <- vars
  m
})

tog

# age mean and std for all samples together
mean(current_comb$Age)
sd(current_comb$Age)

# Combine current ptsd
cur_mean_sd <- Map(function(x, y) merge(x, y, by =1, all.x = T ),
                   sep[[1]], tog[[1]])
cur_mean_sd


# combine lifetime ptsd
life_mean_sd <- Map(function(x, y) merge(x, y, by =1, all.x = T ),
                    sep[[2]], tog[[2]])
life_mean_sd

# combine current and lifetime ptsd
comb <- Map(function(x, y) merge(x , y , by = 1, all.x = T),
            cur_mean_sd, life_mean_sd)
comb


# sort columns
col_order <- c("Study", "Ptsdpm1", "Ptsdpm0",
               "Ptsdpm_t.test_p", "Ptsdpm_mann_whitney_p",
               "Ptsdpm_mean_sd", "Ptsdlife1", "Ptsdlife0",
               "Ptsdlife_t.test_p", "Ptsdlife_mann_whitney_p",
               "Ptsdlife_mean_sd")



sorted <- lapply(seq_along(comb), function(i){
  n <- names(comb)[i]
  x <- comb[[i]][, col_order]
  x$Variable <- n
  x
})
sorted

sorted <- do.call(rbind, sorted)
sorted

write.xlsx(sorted, file = f_path, append = TRUE, sheetName = "stats")



# gender ----------------------

get_category <- function(df, group1, count_col,
                       group2){

  # function to get information of independent cohorts
  get_cur_gen <- function(){
    df %>%
      group_by(across(all_of(group1))) %>%
      count(across(all_of(count_col))) %>%
      group_by(across(all_of(group2))) %>%
      mutate(percentage = round(prop.table(n), 3)*100)
  }
  cur_gen <- get_cur_gen()

  # group on PTSD and get information of all cohorts together
  group1 <- group1[2] # current ptsd
  all_cur_gen <- df %>%
    group_by(across(all_of(group1))) %>%
    count(across(all_of(count_col))) %>%
    mutate(percentage = round(prop.table(n), 3)*100) %>%
    mutate("Study" = "All cohorts")

  cur_gen <- rbind(cur_gen,all_cur_gen)

  print(cur_gen)
  print(all_cur_gen)

  count_all <- function(){
    df %>%
    count(across(all_of(count_col))) %>%
    mutate(percentage = round(prop.table(n), 3)*100) %>%
    mutate(total = paste0(n, " (", percentage, ")"))
  }

  # all cohorts
  all_tog <- count_all() %>%
    mutate("Study" = "All cohorts")

  # individual cohorts
  df <- df %>% group_by(across(all_of(group2)))
  all_indv <- count_all()

  total <- rbind(all_indv, all_tog)
  print(total)

  list(cur_gen, total)
}

# get total and percentage of males and females in whole dataset
table(current_comb$Gender)
table(current_comb$Gender)/nrow(current_comb)*100

group1 = c("Study", "current_ptsd")
cur_gen <- get_category(df = current_comb, group1 = group1,
           count_col = "Gender", group2 = "Study")
cur_gen

cur_gen[[1]] <- cur_gen[[1]] %>% filter(Gender == 1) %>%
  pivot_wider(names_from = c(current_ptsd, Gender),
              values_from = c(n, percentage))

colnames(cur_gen[[1]]) <- c("Study", paste("current",
                                           colnames(cur_gen[[1]])[-1], sep = "_"))

cur_gen[[1]] <- cur_gen[[1]] %>%
  mutate(M_case = paste0(current_n_1_1, " (", current_percentage_1_1, ")"),
         M_control = paste0(current_n_0_1, " (", current_percentage_0_1, ")"))
cur_gen[[2]] <- cur_gen[[2]] %>% filter(Gender == 1)

cur_gen <- merge(cur_gen[[1]], cur_gen[[2]][,c(1,5) ],
                 by = 1,
                 all.x = TRUE)
cur_gen

write.xlsx(cur_gen, file = f_path, append = TRUE, sheetName = "gender_ptsdpm")


# lifetime
group1 = c("Study", "lifetime_ptsd")
life_gen <- get_category(df = lifetime_comb, group1 = group1,
                      count_col = "Gender", group2 = "Study")
life_gen

life_gen[[1]] <- life_gen[[1]] %>% filter(Gender == 1) %>%
  pivot_wider(names_from = c(lifetime_ptsd, Gender),
              values_from = c(n, percentage))

colnames(life_gen[[1]]) <- c("Study", paste("current",
                                            colnames(life_gen[[1]])[-1], sep = "_"))

life_gen[[1]] <- life_gen[[1]] %>%
  mutate(M_case = paste0(current_n_1_1, " (", current_percentage_1_1, ")"),
         M_control = paste0(current_n_0_1, " (", current_percentage_0_1, ")"))
life_gen[[2]] <- life_gen[[2]] %>% filter(Gender == 1)

life_gen <- merge(life_gen[[1]], life_gen[[2]][,c(1,5) ],
                 by = 1,
                 all.x = TRUE)
life_gen

write.xlsx(life_gen, file = f_path, append = TRUE, sheetName = "gender_ptsdlife")


# Race ------------------------------
#
# Armystarrs -
# DNHS - 1	White, 2	Black or African American, 3	Other
# GTP - 2 African American, 1 others
# MRS - 1‐Black or African American, 2=American Indian or Alaskan Native,
#       3=Asian, 4=Native Hawaiian or Pacific Islander, 5=White, #N/A = missing

# PRISMO


get_race_cat <- function(df){
  df %>%
    mutate(
      Race_cat = case_when(
        Study == 'DNHS' & Race == 2 ~ 'African American',
        Study == 'DNHS' & Race != 2 ~ 'Other',
        Study == 'Armystarrs' & Race == 'White' ~ 'White',
        Study == 'Armystarrs' & Race == 'Black' ~ 'African American',
        Study == 'Armystarrs' & Race != 'White|Black' ~ 'Other',
        Study == 'GTP' & Race == '2.0' ~ 'African American',
        Study == 'GTP' & Race != '2.0' ~ 'Other',
        Study == 'MRS' & Race == 5 ~ 'White',
        Study == 'MRS' & Race == '1' ~ 'African American',
        Study == 'MRS' & Race != '1|5' ~ 'Other',
        Study == 'Prismo' & Race == 'eur' ~ 'White',
        Study == 'Prismo' & Race == 'aam' ~ 'African American',
        Study == 'Prismo' & Race != 'eur|aam'  ~ 'Other',

      ))
}

current <- get_race_cat(df = current_comb)
View(current)

# get race in all samples
table(current$Race_cat)
round(table(current$Race_cat)/nrow(current), 2)*100

lifetime <- get_race_cat(df = lifetime_comb)
View(lifetime)


# View(current %>% group_by(Study, current_ptsd) %>% count(Race_cat))
#
# lifetime_comb %>% group_by(Study) %>% count(Race)
#
# count_race <- function(df, group){
#   df %>%
#     group_by(across(all_of(group))) %>%
#     count(Race_cat) %>%
#     arrange(Study, desc(n)) %>%
#     mutate(percentage = round(prop.table(n), 3)*100) %>%
#     mutate(N_per = paste0(n, " (", percentage,")"))
# }
#
# current_race <- count_race(df = current,
#                            group = c("Study", "current_ptsd"))
# View(current_race)
#
# life_race <- count_race(df = lifetime,
#                         group = c("Study", "lifetime_ptsd") )
# View(life_race)
#
# comb_race <- merge(current_race, life_race,
#                    by = c("Study", "Race_cat"),
#                    all.x = T)
# View(comb_race)
#
#
# colnames(comb_race) <- sub("\\.x", "_current", colnames(comb_race))
# colnames(comb_race) <- sub("\\.y", "_lifetime", colnames(comb_race))
#
# comb_race <- comb_race %>% arrange(Study, desc(n_current))

group1 = c("Study", "current_ptsd")
cur_race <- get_category(df = current, group1 = group1,
                        count_col = "Race_cat", group2 = "Study")
cur_race


cur_race[[1]] <- cur_race[[1]] %>%
  pivot_wider(names_from = c(current_ptsd),
              values_from = c(n, percentage))

colnames(cur_race[[1]]) <- c("Study", paste("current",
                                           colnames(cur_race[[1]])[-1], sep = "_"))

cur_race[[1]] <- cur_race[[1]] %>%
  mutate(M_case = paste0(current_n_1, " (", current_percentage_1, ")"),
         M_control = paste0(current_n_0, " (", current_percentage_0, ")"))

cur_race <- merge(cur_race[[1]], cur_race[[2]][,c(1,2,5) ],
                 by = c(1,2),
                 all.x = TRUE)
cur_race <- cur_race %>% arrange(Study, desc(current_n_0))
cur_race


# lifetime race
grouplife = c("Study", "lifetime_ptsd")
lifetime_race <- get_category(df = lifetime, group1 = grouplife,
                         count_col = "Race_cat", group2 = "Study")
lifetime_race


lifetime_race[[1]] <- lifetime_race[[1]] %>%
  pivot_wider(names_from = c(lifetime_ptsd),
              values_from = c(n, percentage))

colnames(lifetime_race[[1]]) <- c("Study", paste("lifetime",
                                            colnames(lifetime_race[[1]])[-1], sep = "_"))

lifetime_race[[1]] <- lifetime_race[[1]] %>%
  mutate(M_case = paste0(lifetime_n_1, " (", lifetime_percentage_1, ")"),
         M_control = paste0(lifetime_n_0, " (", lifetime_percentage_0, ")"))

lifetime_race <- merge(lifetime_race[[1]], lifetime_race[[2]][,c(1,2,5) ],
                  by = c(1,2),
                  all.x = TRUE)
lifetime_race <- lifetime_race %>% arrange(Study, desc(lifetime_n_0))
lifetime_race

race <- merge(cur_race, lifetime_race,
              by = c(1,2),
              all.x = T) %>%
  arrange(Study, desc(current_n_0))
View(race)

# add empty rows so that we can easily copy the table in word
library(purrr)
race <- race %>%
  group_split(Study) %>%
  map_dfr(~ .x %>%
            add_row(Study = .x$Study[1], .before = 1))

# add pval colums as - for easy copying
race <- race %>%
  mutate(cur_pval = '-', .before = total.x) %>%
  mutate(life_pval = '-', .before = total.y)

write.xlsx(race, file = f_path, append = TRUE, sheetName = "race", showNA = FALSE)


# lets subset data with
current_comb_sub <- current_comb %>%
  select(Basename, Traumanum, Childhood_Mt, Study)
View(current_comb_sub)

grp_range <- current_comb_sub %>%
  group_by(Study) %>%
  summarise(mean_trauma = mean(Traumanum),
            std_trauma = sd(Traumanum),
            mean_chd_trauma = mean(Childhood_Mt),
            sd_chd_trauma = sd(Childhood_Mt)
            )
grp_range
