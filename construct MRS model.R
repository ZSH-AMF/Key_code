##############construct MRS model
data$sex <- as.factor(data$sex)
data$Neoadjuvant_chemoradiotherapy <- as.factor(data$Neoadjuvant_chemoradiotherapy)
data$Tumor_stage <- ifelse(data$Tumor_stage == 0, 0,
                           ifelse(data$Tumor_stage == "1",0,
                                  ifelse(data$Tumor_stage == "2",1,
                                         ifelse(data$Tumor_stage == "3",2,
                                                ifelse(data$Tumor_stage == "4",3,9)))))
data$Tumor_stage <- as.factor(data$Tumor_stage)
data_train <- data %>% group_by(PFS) %>% sample_frac(size = .50)
sample_ids_train <- unique(data_train$Class)
a<-data.frame(sample_ids_train)
data_train <- as.data.frame(data_train)
data_test <-  anti_join(data, data_train)
data_test  <- as.data.frame(data_test)
data$Binary_N_score <- ifelse(data$Binary_N > 0, 1, 0)
data$Binary_T_score <- ifelse(data$Binary_T > 0, 0, 1)
data$Binary_C_score <- ifelse(data$Binary_C > 0, 1, 0)

data$MHS <- data$Binary_N + data$Binary_T + 
  data$Binary_C
set.seed(100)

c_indices <- numeric(1000)

for (i in 1:1000) {

  data_train <- data %>% group_by(PFS) %>% sample_frac(size = 0.50)
  sample_ids_train <- unique(data_train$Class)
  

  a <- data.frame(sample_ids_train)
  data_train <- as.data.frame(data_train)
  data_test <- anti_join(data, data_train, by = "Class")
  data_test <- as.data.frame(data_test)

  cox_model <- coxph(Surv(PFS_time, PFS) ~Binary_N+Binary_T+Binary_C, 
                     data = data_train)
  

  predicted_risks_MHS <- predict(cox_model, newdata = data_train, type = "risk")
  

  surv_obj_test <- Surv(data_train$PFS_time, data_train$PFS)
  

  c_index <- survConcordance(surv_obj_test ~  predicted_risks_MHS)$concordance
  

  c_indices[i] <- c_index
  

  print(paste("C-index for sample", i, ":", c_index))
}


mean_c_index <- mean(c_indices)
sd_c_index <- sd(c_indices)
se_c_index <- sd_c_index / sqrt(length(c_indices))

min_c_index <- min(c_indices)
max_c_index <- max(c_indices)

ci_lower <- quantile(c_indices, 0.025)
ci_upper <- quantile(c_indices, 0.975)

print(paste("Average C-index over", length(c_indices), "samples:", mean_c_index))
print(paste("Standard deviation of C-index:", sd_c_index))
print(paste("95% CI for C-index: [", ci_lower, ",", ci_upper, "]"))
print(paste("Minimum C-index over", length(c_indices), "samples:", min_c_index))
print(paste("Maximum C-index over", length(c_indices), "samples:", max_c_index))
###################################
data$Binary_N <- ifelse(data$Binary_N > 0, 1, 0)
data$Binary_T <- ifelse(data$Binary_T > 0, 0, 1)
data$Binary_C <- ifelse(data$Binary_C > 0, 1, 0)

data$MHS <- data$Binary_N + data$Binary_T + 
  data$Binary_C

table(data$PFS,data$MHS)
data$MHS <- ifelse(data$MHS == "0", 0,
                   ifelse(data$MHS == "1", 1,
                          ifelse(data$MHS == "2", 2,
                                 ifelse(data$MHS == "3", 2, NA)))) 

write.csv(data,file = "E:/data.csv")

data$MHS <- as.factor(data$MHS)

table(data$Tumor_stage,data$MHS)


fit <- survfit(Surv(PFS_time, PFS) ~ MHS, data = data)

ggsurvplot(
  fit,
  data = data,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  surv.median.line = "hv",
  palette = c("#7BC8F6", "#FFC0CB","#FA8072"), 
  ggtheme = theme_minimal(), 
  legend.title = "Group",
  legend.labs = c("Low", "Moderate","High"),
  xlab = "Time", 
  ylab = "PFS" 
)
####OS
fit <- survfit(Surv(OS_time, OS) ~ MHS, data = data)

ggsurvplot(
  fit,
  data = data,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  surv.median.line = "hv",
  palette = c("#7BC8F6", "#FFC0CB","#FA8072"), 
  ggtheme = theme_minimal(), 
  legend.title = "Group",
  legend.labs = c("Low", "Moderate","High"),
  xlab = "Time", 
  ylab = "OS" 
)