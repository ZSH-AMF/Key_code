#################machine learning#RSF
###############################################
########################################################### 
data_set <- data_train[ , c(87, 94,8,61,67)]
# 
cv_rfsrc <- function(data_set, formula, ntree = 100, nodesize, mtry, nsplit = 10, seed = 1000) {
  set.seed(seed)  # 
  folds <- createFolds(data_set$PFS, k = nsplit, returnTrain = FALSE)
  
  c_index <- numeric(nsplit)
  models <- vector("list", nsplit)  # 
  
  for (i in seq_along(folds)) {
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_len(nrow(data_set)), test_idx)
    
    train_data <- data_set[train_idx, ]
    test_data <- data_set[test_idx, ]
    
    model <- rfsrc(formula, data = train_data, ntree = ntree, nodesize = nodesize, 
                   mtry = mtry, splitrule = "logrank", importance = TRUE, proximity = TRUE, forest = TRUE, seed = seed)
    
    models[[i]] <- model  # 
    
    predictions <- predict(model, newdata = test_data)
    
    # 
    surv_obj <- Surv(test_data$PFS_time, test_data$PFS)
    
    # 
    c_index[i] <- survConcordance(surv_obj ~ predictions$predicted)$concordance
  }
  
  # 
  c_index_mean <- mean(c_index)
  c_index_sd <- sd(c_index)
  
  # 
  return(list(mean = c_index_mean, sd = c_index_sd, models = models))
}

# 
data_set <- data_train[,c(87, 94,8,61,67)]
formula <- Surv(PFS_time, PFS) ~ .

# 
nodesize_values <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
mtry_values <- c(2, 3, 4, 5, 6, 7, 8)

results <- expand.grid(nodesize = nodesize_values, mtry = mtry_values)
results$mean_c_index <- numeric(nrow(results))
results$sd_c_index <- numeric(nrow(results))
best_models <- list()  # 

# 
for (i in seq_len(nrow(results))) {
  nodesize <- results$nodesize[i]
  mtry <- results$mtry[i]
  cv_results <- cv_rfsrc(data_set, formula, ntree = 100, nodesize = nodesize, mtry = mtry, nsplit = 10, seed = 1000)
  results$mean_c_index[i] <- cv_results$mean
  results$sd_c_index[i] <- cv_results$sd
  best_models[[i]] <- cv_results$models  # 
}

# 
print(results)

# 
best_idx <- which.max(results$mean_c_index)
best_nodesize <- results$nodesize[best_idx]
best_mtry <- results$mtry[best_idx]
print(paste("Best nodesize:", best_nodesize))
print(paste("Best mtry:", best_mtry))

# 
best_model_com <- best_models[[best_idx]][[1]]  # 

# 
predictions <- predict(best_model_com, newdata = data_test)

# 
surv_obj <- Surv(data_test$PFS_time, data_test$PFS)
c_index_total <- survConcordance(surv_obj ~ predictions$predicted)$concordance

# 
print(paste("C-index for the test dataset:", c_index_total))

# 
predictions <- predict(best_model_com, newdata = data)
# 
surv_obj <- Surv(data$PFS_time, data$PFS)
c_index_total <- survConcordance(surv_obj ~ predictions$predicted)$concordance
print(paste("C-index for the test dataset:", c_index_total))

####
best_model_com <- best_models[[best_idx]][[1]]

# 
bootstrap_c_index <- function(data_train, indices) {
  d <- data_train[indices, ]
  predictions <- predict(best_model_com, newdata = d)
  surv_obj <- Surv(d$PFS_time, d$PFS)
  return(survConcordance(surv_obj ~ predictions$predicted)$concordance)
}

# 
set.seed(123)
boot_results <- boot(data_train, bootstrap_c_index, R = 1000)

# 
c_index_mean <- mean(boot_results$t)
c_index_ci <- boot.ci(boot_results, type = "perc")$percent[4:5]

# 
print(paste("Mean C-index:", c_index_mean))
print(paste("95% CI:", c_index_ci[1], "-", c_index_ci[2]))

####
best_model_com <- best_models[[best_idx]][[1]]

# 
bootstrap_c_index <- function(data_test, indices) {
  d <- data_test[indices, ]
  predictions <- predict(best_model_com, newdata = d)
  surv_obj <- Surv(d$PFS_time, d$PFS)
  return(survConcordance(surv_obj ~ predictions$predicted)$concordance)
}

# 
set.seed(123)
boot_results <- boot(data_test, bootstrap_c_index, R = 1000)

# 
c_index_mean <- mean(boot_results$t)
c_index_ci <- boot.ci(boot_results, type = "perc")$percent[4:5]

# 
print(paste("Mean C-index:", c_index_mean))
print(paste("95% CI:", c_index_ci[1], "-", c_index_ci[2]))

####
best_model_com <- best_models[[best_idx]][[1]]

# 
bootstrap_c_index <- function(data, indices) {
  d <- data[indices, ]
  predictions <- predict(best_model_com, newdata = d)
  surv_obj <- Surv(d$PFS_time, d$PFS)
  return(survConcordance(surv_obj ~ predictions$predicted)$concordance)
}

# 
set.seed(123)
boot_results <- boot(data, bootstrap_c_index, R = 1000)

# 
c_index_mean <- mean(boot_results$t)
c_index_ci <- boot.ci(boot_results, type = "perc")$percent[4:5]

# 
print(paste("Mean C-index:", c_index_mean))
print(paste("95% CI:", c_index_ci[1], "-", c_index_ci[2]))
###############################################
data_set <- data_train[ , c(87, 94,266)]
# 
cv_rfsrc <- function(data_set, formula, ntree = 100, nodesize, mtry, nsplit = 10, seed = 1000) {
  set.seed(seed)  # 
  folds <- createFolds(data_set$PFS, k = nsplit, returnTrain = FALSE)
  
  c_index <- numeric(nsplit)
  models <- vector("list", nsplit)  # 
  
  for (i in seq_along(folds)) {
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_len(nrow(data_set)), test_idx)
    
    train_data <- data_set[train_idx, ]
    test_data <- data_set[test_idx, ]
    
    model <- rfsrc(formula, data = train_data, ntree = ntree, nodesize = nodesize, 
                   mtry = mtry, splitrule = "logrank", importance = TRUE, proximity = TRUE, forest = TRUE, seed = seed)
    
    models[[i]] <- model  # 
    
    predictions <- predict(model, newdata = test_data)
    
    # 
    surv_obj <- Surv(test_data$PFS_time, test_data$PFS)
    
    # 
    c_index[i] <- survConcordance(surv_obj ~ predictions$predicted)$concordance
  }
  
  # 
  c_index_mean <- mean(c_index)
  c_index_sd <- sd(c_index)
  
  # 
  return(list(mean = c_index_mean, sd = c_index_sd, models = models))
}

# 
data_set <- data_train[,c(87, 94,266)]
formula <- Surv(PFS_time, PFS) ~ .

# 
nodesize_values <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
mtry_values <- c(2, 3, 4, 5, 6, 7, 8)

results <- expand.grid(nodesize = nodesize_values, mtry = mtry_values)
results$mean_c_index <- numeric(nrow(results))
results$sd_c_index <- numeric(nrow(results))
best_models <- list()  # 

# 
for (i in seq_len(nrow(results))) {
  nodesize <- results$nodesize[i]
  mtry <- results$mtry[i]
  cv_results <- cv_rfsrc(data_set, formula, ntree = 100, nodesize = nodesize, mtry = mtry, nsplit = 10, seed = 1000)
  results$mean_c_index[i] <- cv_results$mean
  results$sd_c_index[i] <- cv_results$sd
  best_models[[i]] <- cv_results$models  # 
}

# 
print(results)

# 
best_idx <- which.max(results$mean_c_index)
best_nodesize <- results$nodesize[best_idx]
best_mtry <- results$mtry[best_idx]
print(paste("Best nodesize:", best_nodesize))
print(paste("Best mtry:", best_mtry))

# 
best_model_com <- best_models[[best_idx]][[1]]  # 

# 
predictions <- predict(best_model_com, newdata = data_test)

# 
surv_obj <- Surv(data_test$PFS_time, data_test$PFS)
c_index_total <- survConcordance(surv_obj ~ predictions$predicted)$concordance

# 
print(paste("C-index for the test dataset:", c_index_total))

# 
predictions <- predict(best_model_com, newdata = data)
# 
surv_obj <- Surv(data$PFS_time, data$PFS)
c_index_total <- survConcordance(surv_obj ~ predictions$predicted)$concordance
print(paste("C-index for the test dataset:", c_index_total))

####
best_model_com <- best_models[[best_idx]][[1]]

# 
bootstrap_c_index <- function(data_train, indices) {
  d <- data_train[indices, ]
  predictions <- predict(best_model_com, newdata = d)
  surv_obj <- Surv(d$PFS_time, d$PFS)
  return(survConcordance(surv_obj ~ predictions$predicted)$concordance)
}

# 
set.seed(123)
boot_results <- boot(data_train, bootstrap_c_index, R = 1000)

# 
c_index_mean <- mean(boot_results$t)
c_index_ci <- boot.ci(boot_results, type = "perc")$percent[4:5]

# 
print(paste("Mean C-index:", c_index_mean))
print(paste("95% CI:", c_index_ci[1], "-", c_index_ci[2]))

####
best_model_com <- best_models[[best_idx]][[1]]

# 
bootstrap_c_index <- function(data_test, indices) {
  d <- data_test[indices, ]
  predictions <- predict(best_model_com, newdata = d)
  surv_obj <- Surv(d$PFS_time, d$PFS)
  return(survConcordance(surv_obj ~ predictions$predicted)$concordance)
}

# 
set.seed(123)
boot_results <- boot(data_test, bootstrap_c_index, R = 1000)

# 
c_index_mean <- mean(boot_results$t)
c_index_ci <- boot.ci(boot_results, type = "perc")$percent[4:5]

# 
print(paste("Mean C-index:", c_index_mean))
print(paste("95% CI:", c_index_ci[1], "-", c_index_ci[2]))

####
best_model_com <- best_models[[best_idx]][[1]]

# 
bootstrap_c_index <- function(data, indices) {
  d <- data[indices, ]
  predictions <- predict(best_model_com, newdata = d)
  surv_obj <- Surv(d$PFS_time, d$PFS)
  return(survConcordance(surv_obj ~ predictions$predicted)$concordance)
}

# 
set.seed(123)
boot_results <- boot(data, bootstrap_c_index, R = 1000)

# 
c_index_mean <- mean(boot_results$t)
c_index_ci <- boot.ci(boot_results, type = "perc")$percent[4:5]

# 
print(paste("Mean C-index:", c_index_mean))
print(paste("95% CI:", c_index_ci[1], "-", c_index_ci[2]))
#################################
data_set <- data_train[ , c(87, 94,8,61,67,266)]
# 
cv_rfsrc <- function(data_set, formula, ntree = 100, nodesize, mtry, nsplit = 10, seed = 1000) {
  set.seed(seed)  # 
  folds <- createFolds(data_set$PFS, k = nsplit, returnTrain = FALSE)
  
  c_index <- numeric(nsplit)
  models <- vector("list", nsplit)  # 
  
  for (i in seq_along(folds)) {
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_len(nrow(data_set)), test_idx)
    
    train_data <- data_set[train_idx, ]
    test_data <- data_set[test_idx, ]
    
    model <- rfsrc(formula, data = train_data, ntree = ntree, nodesize = nodesize, 
                   mtry = mtry, splitrule = "logrank", importance = TRUE, proximity = TRUE, forest = TRUE, seed = seed)
    
    models[[i]] <- model  # 
    
    predictions <- predict(model, newdata = test_data)
    
    # 
    surv_obj <- Surv(test_data$PFS_time, test_data$PFS)
    
    # 
    c_index[i] <- survConcordance(surv_obj ~ predictions$predicted)$concordance
  }
  
  # 
  c_index_mean <- mean(c_index)
  c_index_sd <- sd(c_index)
  
  # 
  return(list(mean = c_index_mean, sd = c_index_sd, models = models))
}

# 
data_set <- data_train[,c(87, 94,8,61,67,266)]
formula <- Surv(PFS_time, PFS) ~ .

# 
nodesize_values <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
mtry_values <- c(2, 3, 4, 5, 6, 7, 8)

results <- expand.grid(nodesize = nodesize_values, mtry = mtry_values)
results$mean_c_index <- numeric(nrow(results))
results$sd_c_index <- numeric(nrow(results))
best_models <- list()  # 

# 
for (i in seq_len(nrow(results))) {
  nodesize <- results$nodesize[i]
  mtry <- results$mtry[i]
  cv_results <- cv_rfsrc(data_set, formula, ntree = 100, nodesize = nodesize, mtry = mtry, nsplit = 10, seed = 1000)
  results$mean_c_index[i] <- cv_results$mean
  results$sd_c_index[i] <- cv_results$sd
  best_models[[i]] <- cv_results$models  # 
}

# 
print(results)

# 
best_idx <- which.max(results$mean_c_index)
best_nodesize <- results$nodesize[best_idx]
best_mtry <- results$mtry[best_idx]
print(paste("Best nodesize:", best_nodesize))
print(paste("Best mtry:", best_mtry))

# 
best_model_com <- best_models[[best_idx]][[1]]  # 

# 
predictions <- predict(best_model_com, newdata = data_test)

# 
surv_obj <- Surv(data_test$PFS_time, data_test$PFS)
c_index_total <- survConcordance(surv_obj ~ predictions$predicted)$concordance

# 
print(paste("C-index for the test dataset:", c_index_total))

# 
predictions <- predict(best_model_com, newdata = data)
# 
surv_obj <- Surv(data$PFS_time, data$PFS)
c_index_total <- survConcordance(surv_obj ~ predictions$predicted)$concordance
print(paste("C-index for the test dataset:", c_index_total))

####
best_model_com <- best_models[[best_idx]][[1]]

# 
bootstrap_c_index <- function(data_train, indices) {
  d <- data_train[indices, ]
  predictions <- predict(best_model_com, newdata = d)
  surv_obj <- Surv(d$PFS_time, d$PFS)
  return(survConcordance(surv_obj ~ predictions$predicted)$concordance)
}

# 
set.seed(123)
boot_results <- boot(data_train, bootstrap_c_index, R = 1000)

# 
c_index_mean <- mean(boot_results$t)
c_index_ci <- boot.ci(boot_results, type = "perc")$percent[4:5]

# 
print(paste("Mean C-index:", c_index_mean))
print(paste("95% CI:", c_index_ci[1], "-", c_index_ci[2]))

####
best_model_com <- best_models[[best_idx]][[1]]

# 
bootstrap_c_index <- function(data_test, indices) {
  d <- data_test[indices, ]
  predictions <- predict(best_model_com, newdata = d)
  surv_obj <- Surv(d$PFS_time, d$PFS)
  return(survConcordance(surv_obj ~ predictions$predicted)$concordance)
}

# 
set.seed(123)
boot_results <- boot(data_test, bootstrap_c_index, R = 1000)

# 
c_index_mean <- mean(boot_results$t)
c_index_ci <- boot.ci(boot_results, type = "perc")$percent[4:5]

# 
print(paste("Mean C-index:", c_index_mean))
print(paste("95% CI:", c_index_ci[1], "-", c_index_ci[2]))

####
best_model_com <- best_models[[best_idx]][[1]]

# 
bootstrap_c_index <- function(data, indices) {
  d <- data[indices, ]
  predictions <- predict(best_model_com, newdata = d)
  surv_obj <- Surv(d$PFS_time, d$PFS)
  return(survConcordance(surv_obj ~ predictions$predicted)$concordance)
}

# 
set.seed(123)
boot_results <- boot(data, bootstrap_c_index, R = 1000)

# 
c_index_mean <- mean(boot_results$t)
c_index_ci <- boot.ci(boot_results, type = "perc")$percent[4:5]

# 
print(paste("Mean C-index:", c_index_mean))
print(paste("95% CI:", c_index_ci[1], "-", c_index_ci[2]))