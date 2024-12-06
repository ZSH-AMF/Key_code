#####################identification of prognostic bacteria biomarkers
###
library(vegan)
library(survival)
library(survminer)
library(randomForestSRC)
library(timeROC)
library(data.table)
library(tidyr)
library(dplyr)
library(tibble)
library(caret)
library(dplyr)
library(pec)
library(ggRandomForests)
library(rms)
library("timeROC")
library("compositions")
library(boot)
library(ggplot2)
##############
#######################################################
#######
data <- read.csv("E:/F盘/口腔菌群课题/整理文章结构/9.2/merged_data_species_RMNA.csv", header = TRUE,fileEncoding = "GBK")
data$sex <- as.factor(data$sex)
data$Neoadjuvant_chemoradiotherapy <- as.factor(data$Neoadjuvant_chemoradiotherapy)
data$Tumor_stage <- ifelse(data$Tumor_stage == "0", 0,
                           ifelse(data$Tumor_stage == "1", 0,
                                  ifelse(data$Tumor_stage == "2", 1,
                                         ifelse(data$Tumor_stage == "3", 2,
                                                ifelse(data$Tumor_stage == "4", 3, 9)))))
data$Lymphatic_invasion <- as.factor(data$Lymphatic_invasion)
data$Perineural_Invasion <- as.factor(data$Perineural_Invasion)
data$Tumor_stage <- as.factor(data$Tumor_stage)
microbe_data <- data[, 118:262]

# 
cox_results1 <- list()

# 
for (i in 1:ncol(microbe_data)) {
  microbe <- microbe_data[, i]
  microbe_name <- colnames(microbe_data)[i]
  
  # 
  surv_obj <- Surv(time = data$PFS_time, event = data$PFS)
  
  # 
  cox_model <- coxph(surv_obj ~ microbe, data = data)
  
  #
  cox_results1[[microbe_name]] <- summary(cox_model)
  
  
  # 
  cox_summary <- summary(cox_model)
  cox_results1[[microbe_name]] <- data.frame(
    Bacteria = microbe_name,
    HR = exp(cox_summary$coefficients[, "coef"]),
    HR_lower_95_CI = cox_summary$conf.int[, "lower .95"],
    HR_upper_95_CI = cox_summary$conf.int[, "upper .95"],
    p_value = cox_summary$coefficients[, "Pr(>|z|)"]
  )}
result_df <- do.call(rbind, cox_results1)
write.csv(result_df, "cox_results1.csv", row.names = FALSE)
##########
microbe_data <- data[, 118:262]

# 
cox_results1 <- list()

#
for (i in 1:ncol(microbe_data)) {
  microbe <- microbe_data[, i]
  microbe_name <- colnames(microbe_data)[i]
  
  # 
  surv_obj <- Surv(time = data$OS_time, event = data$OS)
  
  # 
  cox_model <- coxph(surv_obj ~ microbe, data = data)
  
  # 
  cox_results1[[microbe_name]] <- summary(cox_model)
  
  
  # 
  cox_summary <- summary(cox_model)
  cox_results1[[microbe_name]] <- data.frame(
    Bacteria = microbe_name,
    HR = exp(cox_summary$coefficients[, "coef"]),
    HR_lower_95_CI = cox_summary$conf.int[, "lower .95"],
    HR_upper_95_CI = cox_summary$conf.int[, "upper .95"],
    p_value = cox_summary$coefficients[, "Pr(>|z|)"]
  )}
result_df <- do.call(rbind, cox_results1)
write.csv(result_df, "cox_results_OS.csv", row.names = FALSE)
##############################
data <- read.csv("E:/F盘/口腔菌群课题/整理文章结构/9.2/merged_data_species_RMNA.csv", header = TRUE,fileEncoding = "GBK")
data$sex <- as.factor(data$sex)
data$Neoadjuvant_chemoradiotherapy <- as.factor(data$Neoadjuvant_chemoradiotherapy)
data$Tumor_stage <- ifelse(data$Tumor_stage == "0", 0,
                           ifelse(data$Tumor_stage == "1", 0,
                                  ifelse(data$Tumor_stage == "2", 1,
                                         ifelse(data$Tumor_stage == "3", 2,
                                                ifelse(data$Tumor_stage == "4", 3, 9)))))
data$Lymphatic_invasion <- as.factor(data$Lymphatic_invasion)
data$Perineural_Invasion <- as.factor(data$Perineural_Invasion)
data$Tumor_stage <- as.factor(data$Tumor_stage)
microbe_data <- data[, 118:215]
# 
cox_results1 <- list()

#
for (i in 1:ncol(microbe_data)) {
  microbe <- microbe_data[, i]
  microbe_name <- colnames(microbe_data)[i]
  
  # 
  surv_obj <- Surv(time = data$PFS_time, event = data$PFS)
  
  # 
  cox_model <- coxph(surv_obj ~ microbe+sex+age+Tumor_stage+Neoadjuvant_chemoradiotherapy, data = data)
  
  # 
  cox_results1[[microbe_name]] <- summary(cox_model)
  
  
  # 
  cox_summary <- summary(cox_model)
  cox_results1[[microbe_name]] <- data.frame(
    Bacteria = microbe_name,
    coef_names <- rownames(cox_summary$coefficients),
    HR = exp(cox_summary$coefficients[, "coef"]),
    HR_lower_95_CI = cox_summary$conf.int[, "lower .95"],
    HR_upper_95_CI = cox_summary$conf.int[, "upper .95"],
    p_value = cox_summary$coefficients[, "Pr(>|z|)"]
  )}
result_df <- do.call(rbind, cox_results1)
microbe_results <- result_df %>% filter(coef_names....rownames.cox_summary.coefficients. == "microbe")

# 
microbe_results <- microbe_results %>%
  mutate(FDR_p_value = p.adjust(p_value, method = "fdr"))


# 
final_result_df <- bind_rows(microbe_results)
write.csv(final_result_df, "species_multi_cox_results1_PFS_2_11.17.csv", row.names = FALSE)####
###############################
data <- read.csv("E:/F盘/口腔菌群课题/整理文章结构/9.2/merged_data_species_RMNA.csv", header = TRUE,fileEncoding = "GBK")
data$sex <- as.factor(data$sex)
data$Neoadjuvant_chemoradiotherapy <- as.factor(data$Neoadjuvant_chemoradiotherapy)
data$Tumor_stage <- ifelse(data$Tumor_stage == "0", 0,
                           ifelse(data$Tumor_stage == "1", 0,
                                  ifelse(data$Tumor_stage == "2", 1,
                                         ifelse(data$Tumor_stage == "3", 2,
                                                ifelse(data$Tumor_stage == "4", 3, 9)))))
data$Lymphatic_invasion <- as.factor(data$Lymphatic_invasion)
data$Perineural_Invasion <- as.factor(data$Perineural_Invasion)
data$Tumor_stage <- as.factor(data$Tumor_stage)
microbe_data <- data[, 118:215]
# 
cox_results1 <- list()

# 
for (i in 1:ncol(microbe_data)) {
  microbe <- microbe_data[, i]
  microbe_name <- colnames(microbe_data)[i]
  
  # 
  surv_obj <- Surv(time = data$OS_time, event = data$OS)
  
  # 
  cox_model <- coxph(surv_obj ~ microbe+sex+age+Tumor_stage+Neoadjuvant_chemoradiotherapy, data = data)
  
  # 
  cox_results1[[microbe_name]] <- summary(cox_model)
  

  # 
  cox_summary <- summary(cox_model)
  cox_results1[[microbe_name]] <- data.frame(
    Bacteria = microbe_name,
    coef_names <- rownames(cox_summary$coefficients),
    HR = exp(cox_summary$coefficients[, "coef"]),
    HR_lower_95_CI = cox_summary$conf.int[, "lower .95"],
    HR_upper_95_CI = cox_summary$conf.int[, "upper .95"],
    p_value = cox_summary$coefficients[, "Pr(>|z|)"]
  )}
result_df <- do.call(rbind, cox_results1)
microbe_results <- result_df %>% filter(coef_names....rownames.cox_summary.coefficients. == "microbe")

# 
microbe_results <- microbe_results %>%
  mutate(FDR_p_value = p.adjust(p_value, method = "fdr"))


# 
final_result_df <- bind_rows(microbe_results)
write.csv(final_result_df, "species_multi_cox_results1_OS_2_11.17.csv", row.names = FALSE)
###################
#####################################################################################
set.seed(12345)
data <- read.csv("E:/F盘/口腔菌群课题/整理文章结构/9.2/merged_data_species_RMNA.csv", header = TRUE,fileEncoding = "GBK")
data$sex <- as.factor(data$sex)
data$Neoadjuvant_chemoradiotherapy <- as.factor(data$Neoadjuvant_chemoradiotherapy)
data$Tumor_stage <- ifelse(data$Tumor_stage == 0, 0,
                           ifelse(data$Tumor_stage == "1", 0,
                                  ifelse(data$Tumor_stage == "2", 1,
                                         ifelse(data$Tumor_stage == "3", 2,
                                                ifelse(data$Tumor_stage == "4", 3, 9)))))
data$Tumor_stage <- as.factor(data$Tumor_stage)
data$Lymphatic_invasion <- as.factor(data$Lymphatic_invasion)
data$Perineural_Invasion <- as.factor(data$Perineural_Invasion)
num_iterations <- 1000  # Number of Monte Carlo iterations
results_list <- vector("list", num_iterations)
univariate_results_list <- vector("list", num_iterations)  # 存储单因素分析结果的列表
for (iter in 1:num_iterations) {
  cat("Iteration:", iter, "\n")
  
  # Stratified sampling
  data_train <- data %>% group_by(PFS) %>% sample_frac(size = .50)
  data_train <- as.data.frame(data_train)
  data_test <- anti_join(data, data_train)
  data_test <- as.data.frame(data_test)
  
  # Univariate Cox regression on training set
  microbe_data_train <- data_train[, 118:215]
  cox_results1 <- list()
  
  for (i in 1:ncol(microbe_data_train)) {
    microbe <- microbe_data_train[, i]
    microbe_name <- colnames(microbe_data_train)[i]
    
    surv_obj <- Surv(time = data_train$PFS_time, event = data_train$PFS)
    cox_model <- try(coxph(surv_obj ~ microbe, data = data_train), silent = TRUE)
    
    if (inherits(cox_model, "try-error")) next
    
    cox_summary <- summary(cox_model)
    cox_results1[[microbe_name]] <- data.frame(
      Bacteria = microbe_name,
      HR = exp(cox_summary$coefficients[, "coef"]),
      HR_lower_95_CI = cox_summary$conf.int[, "lower .95"],
      HR_upper_95_CI = cox_summary$conf.int[, "upper .95"],
      p_value = cox_summary$coefficients[, "Pr(>|z|)"]
    )
  }
  
  result_df <- do.call(rbind, cox_results1)
  
  # 
  result_df <- result_df %>%
    mutate(FDR_p_value = p.adjust(p_value, method = "fdr"))
  
  # 
  univariate_results_list[[iter]] <- result_df
  significant_bacteria <- filter(result_df, p_value < 0.05)
  
  if (nrow(significant_bacteria) == 0) {
    cat("No significant bacteria found in iteration:", iter, "\n")
    next
  }
  
  significant_bacteria_names <- significant_bacteria %>% pull(Bacteria)
  selected_data_train <- microbe_data_train %>% dplyr::select(all_of(significant_bacteria_names))
  
  # Multivariate Cox regression on training set
  cox_results <- list()
  
  for (i in 1:ncol(selected_data_train)) {
    microbe <- selected_data_train[, i]
    microbe_name <- colnames(selected_data_train)[i]
    
    surv_obj <- Surv(time = data_train$PFS_time, event = data_train$PFS)
    cox_model <- try(coxph(surv_obj ~ microbe + sex + age + Tumor_stage + Neoadjuvant_chemoradiotherapy, data = data_train), silent = TRUE)
    
    if (inherits(cox_model, "try-error")) next
    
    cox_summary <- summary(cox_model)
    cox_results[[microbe_name]] <- data.frame(
      Bacteria = microbe_name,
      coef_names = rownames(cox_summary$coefficients),
      HR = exp(cox_summary$coefficients[, "coef"]),
      HR_lower_95_CI = cox_summary$conf.int[, "lower .95"],
      HR_upper_95_CI = cox_summary$conf.int[, "upper .95"],
      p_value = cox_summary$coefficients[, "Pr(>|z|)"]
    )
  }
  
  multi_result_df <- do.call(rbind, cox_results)
  significant_bacteria <- filter(multi_result_df, p_value < 0.05)
  coef_names_values <- significant_bacteria$coef_names
  microbe_indices <- which(coef_names_values == "microbe")
  microbe_column_names <- significant_bacteria[microbe_indices, "Bacteria"]
  
  if (length(microbe_column_names) == 0) {
    cat("No significant microbes found for multivariate Cox regression in iteration:", iter, "\n")
    next
  }
  
  # Validation on test set
  microbe_data_test <- data_test[, 118:215]
  selected_data_test <- try(microbe_data_test %>% dplyr::select(all_of(microbe_column_names)), silent = TRUE)
  
  if (inherits(selected_data_test, "try-error")) {
    cat("Error selecting columns in iteration:", iter, "\n")
    next
  }
  
  cox_results_test <- list()
  
  for (i in 1:ncol(selected_data_test)) {
    microbe <- selected_data_test[, i]
    microbe_name <- colnames(selected_data_test)[i]
    
    surv_obj <- Surv(time = data_test$PFS_time, event = data_test$PFS)
    cox_model <- try(coxph(surv_obj ~ microbe + sex + age + Tumor_stage + Neoadjuvant_chemoradiotherapy, data = data_test), silent = TRUE)
    
    if (inherits(cox_model, "try-error")) next
    
    cox_summary <- summary(cox_model)
    cox_results_test[[microbe_name]] <- data.frame(
      Bacteria = microbe_name,
      coef_names = rownames(cox_summary$coefficients),
      HR = exp(cox_summary$coefficients[, "coef"]),
      HR_lower_95_CI = cox_summary$conf.int[, "lower .95"],
      HR_upper_95_CI = cox_summary$conf.int[, "upper .95"],
      p_value = cox_summary$coefficients[, "Pr(>|z|)"]
    )
  }
  
  result_df_test <- do.call(rbind, cox_results_test)
  
  #####
  #
  result_df_test <- result_df_test %>%
    mutate(FDR_p_value = p.adjust(p_value, method = "fdr"))
  
  # 
  final_result_df <- bind_rows(result_df_test)
  
  results_list[[iter]] <- final_result_df  # 保存测试集的多因素分析结
  
}

# Remove NULL elements
results_list <- results_list[!sapply(results_list, is.null)]
univariate_results_list <- univariate_results_list[!sapply(univariate_results_list, is.null)]  # 删除NULL元素

# Aggregate results from all iterations
final_results <- bind_rows(results_list)
final_univariate_results <- bind_rows(univariate_results_list)  # 汇总单因素分析结果

# Output results to CSV
write.csv(final_univariate_results, 
          file = "univariate_results_1000_.csv")