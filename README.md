# PTSD Methylation Risk Score 2024
Code in this repository is to pre-process data for machine learning and train models to predict PTSD and create a methylation risk score.

## Files:

### Helper functions
1. `install_needed_packages.R` install required packages.
2. `DNHS_more_pheno.R`get the required variables for DNHS.
3. `Smoking_Scores_PGC_cohorts.R` estimates smoking scores for each individual discovery cohort.
4. `MRS_Preprocess.R` pre-process Marine Resilience Study cohort to include in the training.
5. `Armystarrs_and_PRISMO_preprocess.R` to pre-process Army STARRS and PRISMO cohorts pre-post deployment samples to test risk scores.
6. `Check_after_updating_pheno.R` and `Check_after_updating_pheno.html` to check the updated phenotype file with the old file. 
7. `cpgassoc2.R` helper function to perform association analysis between each CpG and PTSD.
8. `Covariate_adjustment_1.R` example code to show covariate adjustment.
paper as we thought to make it an Epic data paper.
9. `Compare_Effect_Sizes.Rmd` and `Compare_Effect_Sizes.html` code to compare the effect sizes of discovery and Boston VA cohort for model1.
10. `Demographics.R` code to get demographic information for the manuscript.
11. `Cohort_Information.Rmd` and `Cohort_Information.html` code to get summary information from different cohorts, e.g., variables in each cohort to check data availability. 

### Key Instruction: The weights and features for each of the three risk scores are located in the Data folder. The files are named as follows:

1. eMRS_Model1.xlsx
2. MoRS_Model2.xlsx
3. MoRSAE_Model3.xlsx

### Python code to train and test models:
1. `makedirectory.py` Is to make a directory to store the outcome files from each run.
2. `Settings.ipynb` contains settings for packages and plots.
3. `Preprocess_data_updated_1.ipynb` preprocess all cohorts individually for machine learning. 
4. `pre_post_trauma_processing_v1.ipynb` Is to pre-process the cohorts with pre/post samples and choose post-trauma samples for machine learning.
5. `Imputation_Covariate_adjustment_2.1.ipyn` Code to perform imputation and covariate adjustment.
6. `Imputation_Covariate_adjustment_including_Expo_vaiables_2.1.ipynb` Code to perform imputation and covariate adjustment, including exposure variables.
7. `Feature_Selection_and_training_on_ptsdpm_3.3.ipynb` Feature selection using the covariate-adjusted data (output of step 3).
8. `Feature_Selection_and_training_on_ptsdpm_wd_exp_vars_adjustment_3.3.ipynb` Feature selection using the covariate-adjusted data for exposure variables
 (input is step 4 output).
9. `model_performance_5.5.ipynb` Running model and evaluating the performance (input is step 5 output). 
10. `model_performance_wd_exp_vaars_adjustment_5.5.ipynb` Running model and evaluating the performance with adjusted exposure variables (input is step 6 output).

### R code for downstream analysis:
1. `downstream_analysis_v5.qmd` To estimate risk scores for model 1 and 2 and test the risk scores using the test set in discovery cohorts. `downstream_analysis_v5.html` is the generated report. In steps 2 and 3, we test various data sets such as test set, civilians, military, and males and females to look at various scenarios.
2.  `downstream_analysis_adj_for_Exp_Vars_v5.qmd` is to estimate and test risk scores using model 3 on the test data set. `downstream_analysis_adj_for_Exp_Vars_v5.html` is the generated report.
4. `Test_RiskScores_with&without_exp_vars_wd_logit_6.Rmd` is a clean version of estimating and testing risk scores. It used the point-biserial correlation between binary and continuous variables. Also, we used the logit model to predict PTSD using risk scores. `Test_RiskScores_with&without_exp_vars_wd_logit_6.html` is the generated report. This file was used to generate density, distribution and correlation plots for discovery cohorts. 
5. `Pre_Post_Deployment_eMRS.qmd` and `Pre_Post_Deployment_eMRS.html` to test risk scores pre and post-deployment.
6. `Enrichment_analysis_1.qmd` to perform enrichment analysis of top CpGs from models 1, 2, and 3. Models 1 and 2 have the same set of CpGs.
7. `CpGs_in_previous_studies&ML.R` code to find overlap between identified significant CpGs and previous studies.
8. `Overlap_between_MRS_CpGs_metaanalysis_CpGs_Freeze3.R` to check overlap between identified significant CpGs and PGC EWAS meta-analysis and Freeze3 genes. 
9. `mQTL.qmd`and `mQTL.html` Comparing significant CpGs with [BIOS QTL browser CpGs](https://molgenis26.gcc.rug.nl/downloads/biosqtlbrowser/).

### Code for external cohorts is in R/Independent_Cohort:
1. `Create_sample_data.R` and `Create_sample_data without exp vars.R` code to create sample data with and without exposure variables as an example for external cohorts.
2. `Covariate_Adj_RiskScores_1.R` and `Covariate_Adj_RiskScores_without_exp_vars_1.R` code to estimate risk scores with and without exposure variables, respectively.
3. `Test_RiskScores_with&without_exp_vars_wd_logit_2.Rmd` code to test risk scores and generate plots. 



