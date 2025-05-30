###########################
## Libraries
###########################

library(survival)
library(refund)
library(caret)
library(MASS)
library(grpreg)
library(mgcv)
library(randomForestSRC)
library(ggplot2)
library(scales)
library(patchwork)

############################
## The DataSet
############################

setwd("D:/Final_Year_Project")
load("NHANES_11_14_survPA_TD_Final.RData") 
library(dplyr)

time_mort<-baseline$permth_exm/12 #age in month to year-unit
baseline$time_mort<-time_mort
baseline$Age<-baseline$RIDAGEYR
baseline$Race<-baseline$RIDRETH3

df <- baseline %>%
  select(SEQN, Age, Race, BMI, sex, Mobility, mortstat, 
         diabetes.y, poverty_level, Asthma, 
         Arthritis, heart_failure, coronary_heart_disease, angina, stroke, 
         thyroid, bronchitis, cancer, time_mort, actmat) %>%
  filter(!is.na(mortstat) &
           !is.na(time_mort) &
           !is.na(BMI))

c_index_df <- read.csv(file = 'C_index.csv')
names(df)

####################################
## A Visualization of Survival Data
####################################

# Select relevant data
df_fig <- df %>% 
  select(SEQN, time = time_mort, event = mortstat, age = Age, BMI, 
         race = Race, gender = sex) %>%
  filter(SEQN %in% c(64759, 67293, 73561, 75630, 68808, 82820))

# Custom Colors & Theme
bg_color <- "white"
title_color <- "#333333"
text_color <- "#222222"
event_color <- "red"
censored_color <- "#0072B2"

# Set layout
layout(matrix(c(1:21), 7, 3, byrow = FALSE), widths = c(1, 2, 4))
par(bg = bg_color, family = "serif", col = text_color)

# Title: SEQN
par(mar = c(0, 0, 0, 0))
plot.new()
text(0.5, 0.5, "SEQN", cex = 1.9, font = 2, col = title_color)

# SEQN values
par(mar = c(1, 0, 1, 1))
for(i in 1:6){
  plot.new()
  text(0.5, 0.5, df_fig$SEQN[i], cex = 1.6, font = 2, col = "#444444")
}

# Title: Predictors
par(mar = c(0, 0, 0, 0))
plot.new()
text(0.5, 0.5, "Predictors", cex = 2.2, font = 2, col = title_color)

# Predictor Info
par(mar = c(1, 0, 1, 1))
for(i in 1:6){
  plot.new()
  text(0.5, 0.6, paste0("Age: ", df_fig$age[i], 
                        ", Gender: ", df_fig$gender[i], 
                        "\nBMI: ", df_fig$BMI[i], 
                        ", Race: ", df_fig$race[i]), 
       cex = 0.9, font = 1, col = "#555555")
}

# Title: Time to Event
par(mar = c(0, 0, 0, 0))
plot.new()
text(0.5, 0.5, "Time to Event", cex = 2.0, font = 2, col = title_color)

# Survival Data Visualization
par(mar = c(0.2, 0.5, 0.2, 0.5))
for(j in 1:6){
  plot(0, 0, type = "n", xlim = c(0, 10), ylim = c(-1, 1), xaxt = "n", 
       yaxt = "n", bty = "n")
  
  # Gridlines for better visibility
  abline(h = 0, col = "gray90", lwd = 1)
  abline(v = seq(0, 10, 2), col = "gray90", lwd = 0.5, lty = 2)
  
  # Timeline
  segments(0, 0, df_fig$time[j], 0, lwd = 4, col = ifelse(df_fig$event[j] == 1, 
                                                          event_color, 
                                                          censored_color))
  
  # Event Marker: Make Deceased (Cross) More Visible
  if(df_fig$event[j] == 1) {
    # White border for better visibility
    points(df_fig$time[j], 0, pch = 4, col = "white", cex = 2, lwd = 4)  
    # Bold red cross
    points(df_fig$time[j], 0, pch = 4, col = event_color, cex = 2, lwd = 4)  
  } else {
    points(df_fig$time[j], 0, pch = 16, col = censored_color, cex = 2.5)
  }
  
  # **Adjusted Time Label Placement** (Shifted left & up)
  text(df_fig$time[j] + 0.3, 0.5, paste0(format(round(df_fig$time[j], 2), 
                                                nsmall = 2), " years"), 
       cex = 1.4, col = ifelse(df_fig$event[j] == 1, event_color, 
                               censored_color), font = 2)
  
  # Legend (only once)
  if(j == 1){
    legend("topright", legend = c("Censored", "Deceased"), 
           col = c(censored_color, event_color), pch = c(16, 4), pt.cex = 2, 
           lwd = 3, cex = 1.2, bty = "n", y.intersp=1.4)
  }
}

##############################
## FPCA & Data Preprocessing
##############################

tlen <- 10
nt   <- floor(1440/tlen)
tind <- seq(0,1,len=nt)
# create a list of indices for binning into tlen minute windows
inx_col_ls <- split(1:1440, rep(1:nt,each=tlen))
act_mat_bin <- sapply(inx_col_ls, function(x) rowMeans(df$actmat[,x,drop=FALSE],
                                                       na.rm=TRUE))
Y_S<-act_mat_bin
df$Actsm<-Y_S

cat(paste(" Dimensions of actmat:", paste(dim(df$actmat), collapse = " x "),
          "\n","Dimensions of Actsm:", paste(dim(df$Actsm), collapse = " x "), 
          "\n"))

fpca_result <- fpca.sc(Y = df$Actsm, pve = 0.99) 
scoremat<-fpca_result$scores
df$scoremat<-scoremat

tbtick <- seq(0, 1440, by = tlen)
binmid<-tbtick+5
binmid<-binmid[-145]
efmat <- fpca_result$efunctions
explained_variance <- fpca_result$evalues / sum(fpca_result$evalues)  
lwd_values <- 1 + 5 * explained_variance

# Plot all eigenfunctions with scaled line thickness
matplot(binmid / 60, efmat, type = "l", lty = 1, col = rainbow(ncol(efmat)),
        lwd = lwd_values, # Use scaled line widths
        xlab = "Time (Hours)", ylab = "Eigenfunctions", 
        main = "Eigenfunctions Over Time")

# Add a legend with line thickness explanation
legend(x = "topright", inset = c(0.1, 0.05),
       legend = paste("PC", 1:ncol(efmat)), 
       col = rainbow(ncol(efmat)), lty = 1, lwd = lwd_values, cex = 0.5, 
       ncol = 2)

# Set up a 2x2 panel layout
par(mfrow = c(3, 3))

# Define colors for consistency
colors <- rainbow(9)

# Compute common y-axis limits
ylim_range <- range(efmat[, 1:9])

# Loop through the first 4 eigenfunctions and plot them separately
for (i in 1:9) {
  plot(binmid / 60, efmat[, i], type = "l", lty = 1, col = colors[i],
       lwd = 2 + 5 * explained_variance[i], 
       xlab = "Time (Hours)", ylab = paste("PC", i),
       main = paste("Eigenfunction", i), ylim = ylim_range)
}

dfrfs <- subset(df, select = -c(actmat,Actsm,scoremat))

# Automatically add columns for all principal components
for (i in 1:ncol(scoremat)) {
  dfrfs[[paste0("PC", i)]] <- scoremat[, i]
}
# Columns to convert to numeric
cols_to_convert <- c("Race", "sex", "Mobility")

# Convert specified columns to numeric
dfrfs[cols_to_convert] <- lapply(dfrfs[cols_to_convert], as.factor)

##############################
## Cox PH Model
##############################

cox_model <- coxph(Surv(time_mort, mortstat) ~ Age + Race + BMI + sex + 
                     Mobility + diabetes.y + poverty_level + 
                     Asthma + Arthritis + heart_failure + 
                     coronary_heart_disease + angina + stroke + thyroid + 
                     bronchitis + cancer  + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + 
                     PC7 + PC8 + PC9, data = dfrfs)

cox_model_subset = stepAIC(cox_model, direction = 'both')
summary(cox_model_subset)

##############################
## Penalized Cox
##############################

df_pencox <- subset(dfrfs, select = -c(SEQN))
race_dummies <- model.matrix(~ Race - 1, data = df_pencox)  # One-hot encoding

binary_vars <- c("sex", "Mobility", "diabetes.y", "Asthma", "Arthritis", 
                 "heart_failure", "coronary_heart_disease", "angina", "stroke",
                 "thyroid", "bronchitis", "cancer")

df_pencox[binary_vars] <- lapply(df_pencox[binary_vars], 
                                 function(x) as.numeric(as.factor(x)) - 1)

continuous_vars <- c("Age", "BMI", "poverty_level")  # Excluding PC1 to PC9
pc_vars <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9")

X <- cbind(df_pencox[, continuous_vars], df_pencox[, pc_vars], 
           df_pencox[, binary_vars], race_dummies)

groups <- c(1:length(continuous_vars),  # Continuous variables: separate groups
            rep(length(continuous_vars) + 1, length(pc_vars)),  
            (length(continuous_vars) + length(pc_vars) + 1):
              (length(continuous_vars) + length(pc_vars) + length(binary_vars)),
            rep(length(continuous_vars) + length(pc_vars) + length(binary_vars) 
                + 1, 
                ncol(race_dummies)))  # Race: one group

y <- Surv(df_pencox$time_mort, df_pencox$mortstat)  # Survival outcome

fit <- grpsurv(X, y, group = groups, penalty = "grLasso")

cvfit <- cv.grpsurv(X, y, group = groups, penalty = "grLasso")

# Best lambda value
best_lambda <- cvfit$lambda.min

# Coefficients at the optimal lambda
coef(fit, lambda = best_lambda)

###############################
## GAM Cox
###############################

cox_gam <- mgcv::gam(time_mort ~ s(Age,k=20,bs="cr") + s(BMI,k=20,bs="cr") + 
                       Race + sex + Mobility + diabetes.y + 
                       s(poverty_level,k=20,bs="cr") + 
                       Asthma + Arthritis + heart_failure + 
                       coronary_heart_disease + angina + stroke + thyroid + 
                       bronchitis + cancer + s(PC1,k=20,bs="cr") + 
                       s(PC2,k=20,bs="cr") + s(PC3,k=20,bs="cr") + 
                       s(PC4,k=20,bs="cr") + s(PC5,k=20,bs="cr") + 
                       s(PC6,k=20,bs="cr") + s(PC7,k=20,bs="cr") + 
                       s(PC8,k=20,bs="cr") + s(PC9,k=20,bs="cr"), 
                     method = "REML", data = dfrfs, weights = mortstat, 
                     family = cox.ph)

summary(cox_gam)

################################
## RSF
################################

rsf_tuned <- rfsrc(Surv(time_mort,mortstat) ~ Age + Race + BMI + sex + 
                     Mobility + diabetes.y + poverty_level + 
                     Asthma + Arthritis + heart_failure + 
                     coronary_heart_disease + angina + stroke + thyroid + 
                     bronchitis + cancer  + PC1 + PC2 + PC3 + PC4 + PC5 + 
                     PC6 + PC7 + PC8 + PC9,
                   data = dfrfs, ntree = 500, mtry = 3, 
                   nodesize = 15, nsplit = 10, importance = TRUE)

print(rsf_tuned)
plot(rsf_tuned)

#################################
## 100 Split C-Index Comparison
#################################

# Set path to your folder containing the CSV files
data_dir <- "D:/Final_Year_Project/splits_final"

##############################
## Cox PH Model
##############################

calculate_C_index_cox <- function(n_splits = 100, data_path = "./splits", 
                                  train_prefix = "train_split_", 
                                  test_prefix = "test_split_",
                                  factor_cols = c("Race", "sex", "Mobility", 
                                                  "diabetes.y", "Asthma", 
                                                  "Arthritis", 
                                                  "heart_failure", 
                                                  "coronary_heart_disease", 
                                                  "angina", "stroke", 
                                                  "thyroid", "bronchitis", 
                                                  "cancer")) {
  c_index_values <- numeric(n_splits)
  
  for (i in 1:n_splits) {
    train_file <- file.path(data_path, paste0(train_prefix, i, ".csv"))
    test_file  <- file.path(data_path, paste0(test_prefix, i, ".csv"))
    
    df_train <- read.csv(train_file)
    df_test  <- read.csv(test_file)
    
    # Convert specified columns to factors in both train and test sets
    for (col in factor_cols) {
      df_train[[col]] <- as.factor(df_train[[col]])
      df_test[[col]]  <- as.factor(df_test[[col]])
    }
    
    # Fit Cox model
    cox_model <- coxph(Surv(time_mort, mortstat) ~ Age + Race + BMI + sex + 
                         Mobility + diabetes.y + poverty_level + 
                         Asthma + Arthritis + heart_failure + 
                         coronary_heart_disease + angina + stroke + thyroid + 
                         bronchitis + cancer + PC1 + PC2 + PC3 + PC4 + PC5 + 
                         PC6 + PC7 + PC8 + PC9, data = df_train)
    
    # Compute risk scores
    risk_cox <- rowSums(predict(cox_model, df_test, type = "terms"))
    
    # Compute C-index
    c_index_values[i] <- cal_c(marker = risk_cox, Stime = df_test$time_mort, 
                               status = df_test$mortstat)
  }
  
  return(c_index_values)
}

c_index_cox <- calculate_C_index_cox(data_path = data_dir)
mean_c_index_cox <- mean(c_index_cox)

##############################
## Penalized Cox
##############################

calculate_C_index_pencox <- function(n_splits = 100,
                                     data_path = "./splits",
                                     train_prefix = "train_split_",
                                     test_prefix = "test_split_") {
  c_index_values <- numeric(n_splits)
  
  for (i in 1:n_splits) {
    train_file <- file.path(data_path, paste0(train_prefix, i, ".csv"))
    test_file  <- file.path(data_path, paste0(test_prefix, i, ".csv"))
    
    df_train <- read.csv(train_file)
    df_test  <- read.csv(test_file)
    
    # One-hot encode Race (drop intercept)
    race_dummies_train <- model.matrix(~ Race - 1, data = df_train)
    race_dummies_test  <- model.matrix(~ Race - 1, data = df_test)
    
    # Binary variables
    binary_vars <- c("sex", "Mobility", "diabetes.y", "Asthma", "Arthritis", 
                     "heart_failure", "coronary_heart_disease", "angina", 
                     "stroke","thyroid", "bronchitis", "cancer")
    
    df_train[binary_vars] <- lapply(df_train[binary_vars], function(x) 
      as.numeric(as.factor(x)) - 1)
    df_test[binary_vars]  <- lapply(df_test[binary_vars], function(x) 
      as.numeric(as.factor(x)) - 1)
    
    # Continuous and PC features
    continuous_vars <- c("Age", "BMI", "poverty_level")
    pc_vars <- paste0("PC", 1:9)
    
    # Design matrices
    X_train <- cbind(df_train[, continuous_vars],
                     df_train[, pc_vars],
                     df_train[, binary_vars],
                     race_dummies_train)
    
    X_test <- cbind(df_test[, continuous_vars],
                    df_test[, pc_vars],
                    df_test[, binary_vars],
                    race_dummies_test)
    
    # Group structure
    groups <- c(
      1:length(continuous_vars),  # Continuous: separate groups
      rep(length(continuous_vars) + 1, length(pc_vars)),  # PCs: one group
      (length(continuous_vars) + length(pc_vars) + 1):
        (length(continuous_vars) + length(pc_vars) + length(binary_vars)),  
      rep(length(continuous_vars) + length(pc_vars) + length(binary_vars) + 1, 
          ncol(race_dummies_train))  # Race: one group
    )
    
    # Survival outcomes
    y_train <- Surv(df_train$time_mort, df_train$mortstat)
    y_test  <- Surv(df_test$time_mort, df_test$mortstat)
    
    # Fit penalized Cox model
    fit <- grpsurv(X_train, y_train, group = groups, penalty = "grLasso")
    cvfit <- cv.grpsurv(X_train, y_train, group = groups, penalty = "grLasso")
    best_lambda <- cvfit$lambda.min
    
    coefs <- as.numeric(coef(fit, lambda = best_lambda))
    risk_scores <- as.matrix(X_test) %*% coefs
    
    # Compute C-index
    c_index_values[i] <- cal_c(marker = risk_scores,
                               Stime = df_test$time_mort,
                               status = df_test$mortstat)
  }
  
  return(c_index_values)
}

c_index_pencox <- calculate_C_index_pencox(data_path = data_dir)
mean_c_index_pencox <- mean(c_index_pencox)

###############################
## GAM Cox
###############################

calculate_C_index_gamcox <- function(n_splits = 100, data_path = "./splits", 
                                     train_prefix = "train_split_", 
                                     test_prefix = "test_split_",
                                     factor_cols = c("Race", "sex", "Mobility", 
                                                     "diabetes.y", "Asthma", 
                                                     "Arthritis", 
                                                     "heart_failure", 
                                                     "coronary_heart_disease", 
                                                     "angina", "stroke", 
                                                     "thyroid", "bronchitis", 
                                                     "cancer")) {
  c_index_values <- numeric(n_splits)
  
  for (i in 1:n_splits) {
    # Construct full file paths
    train_file <- file.path(data_path, paste0(train_prefix, i, ".csv"))
    test_file  <- file.path(data_path, paste0(test_prefix, i, ".csv"))
    
    # Load train and test data
    df_train <- read.csv(train_file)
    df_test  <- read.csv(test_file)
    
    # Convert specified columns to factors
    for (col in factor_cols) {
      df_train[[col]] <- as.factor(df_train[[col]])
      df_test[[col]]  <- as.factor(df_test[[col]])
    }
    
    # Fit GAM Cox model
    cox_gam <- mgcv::gam(time_mort ~ s(Age, k = 20, bs = "cr") + 
                           s(BMI, k = 20, bs = "cr") + 
                           Race + sex + Mobility + diabetes.y + 
                           s(poverty_level, k = 20, bs = "cr") + 
                           Asthma + Arthritis + heart_failure + 
                           coronary_heart_disease + angina + stroke + 
                           thyroid + bronchitis + cancer + 
                           s(PC1, k = 20, bs = "cr") + 
                           s(PC2, k = 20, bs = "cr") + 
                           s(PC3, k = 20, bs = "cr") + 
                           s(PC4, k = 20, bs = "cr") + 
                           s(PC5, k = 20, bs = "cr") + 
                           s(PC6, k = 20, bs = "cr") + 
                           s(PC7, k = 20, bs = "cr") + 
                           s(PC8, k = 20, bs = "cr") + 
                           s(PC9, k = 20, bs = "cr"),
                         method = "REML", data = df_train, 
                         weights = mortstat, family = cox.ph)
    
    # Predict risk scores
    risk_gamcox <- rowSums(predict(cox_gam, df_test, type = "terms"))
    
    # Compute C-index
    c_index_values[i] <- cal_c(marker = risk_gamcox, Stime = df_test$time_mort, 
                               status = df_test$mortstat)
  }
  
  return(c_index_values)
}

c_index_gamcox <- calculate_C_index_gamcox(data_path = data_dir)
mean_c_index_gamcox <- mean(c_index_gamcox)

################################
## RSF
################################

calculate_C_index_rsf <- function(n_splits = 100, data_path = "./splits", 
                                  train_prefix = "train_split_", 
                                  test_prefix = "test_split_",
                                  factor_cols = c("Race", "sex", "Mobility", 
                                                  "diabetes.y", "Asthma", 
                                                  "Arthritis", 
                                                  "heart_failure", 
                                                  "coronary_heart_disease", 
                                                  "angina", "stroke", 
                                                  "thyroid", "bronchitis", 
                                                  "cancer")) {
  c_index_values <- numeric(n_splits)
  
  for (i in 1:n_splits) {
    # File paths
    train_file <- file.path(data_path, paste0(train_prefix, i, ".csv"))
    test_file  <- file.path(data_path, paste0(test_prefix, i, ".csv"))
    
    # Read train and test data
    df_train <- read.csv(train_file)
    df_test  <- read.csv(test_file)
    
    # Convert specified columns to factors
    for (col in factor_cols) {
      df_train[[col]] <- as.factor(df_train[[col]])
      df_test[[col]]  <- as.factor(df_test[[col]])
    }
    
    # Fit RSF model
    rsf_tuned <- rfsrc(Surv(time_mort, mortstat) ~ Age + Race + BMI + sex + 
                         Mobility + diabetes.y + poverty_level + 
                         Asthma + Arthritis + heart_failure + 
                         coronary_heart_disease + angina + stroke + thyroid + 
                         bronchitis + cancer + PC1 + PC2 + PC3 + PC4 + PC5 + 
                         PC6 + PC7 + PC8 + PC9,
                       data = df_train, ntree = 500, mtry = 3, 
                       nodesize = 15, nsplit = 10, importance = TRUE)
    
    # Determine median death time index for C-index computation
    dis_time <- sort(unique(df_train$time_mort[df_train$mortstat == 1]))
    med_index <- median(1:length(dis_time))
    
    # Predict survival
    pred <- predict(rsf_tuned, newdata = df_test, times = dis_time, se = FALSE)
    surv_probs <- pred$survival
    
    # Compute C-index
    surv_obj <- Surv(df_test$time_mort, df_test$mortstat)
    c_index_values[i] <- Cindex(surv_obj, predicted = surv_probs[, med_index])
  }
  
  return(c_index_values)
}

c_index_rsf <- calculate_C_index_rsf(data_path = data_dir)
mean_c_index_rsf <- mean(c_index_rsf)

results_df <- data.frame(
  c_index_cox = c_index_cox,
  c_index_rsf = c_index_rsf,
  c_index_gamcox = c_index_gamcox,
  c_index_pencox = c_index_pencox
)

#########################################
## Comparison Plot over 100 splits
#########################################

# Load data
c_index_df <- read.csv("c_index_all_collected.csv")
c_index_df$Split <- 1:nrow(c_index_df)

# Compute means
mean_c_index_deepsurv <- mean(c_index_df$c_index_deepsurv)
mean_c_index_cox <- mean(c_index_df$c_index_cox)
mean_c_index_rsf <- mean(c_index_df$c_index_rsf)
mean_c_index_gamcox <- mean(c_index_df$c_index_gamcox)
mean_c_index_pencox <- mean(c_index_df$c_index_pencox)
mean_c_index_bilstm <- mean(c_index_df$c_index_bilstmdeepsurv)

# Prepare colors
default_colors <- hue_pal()(6)
model_colors <- c(
  "DeepSurv Model" = default_colors[1],
  "BiLSTM DeepSurv" = default_colors[2],
  "Cox Model" = default_colors[3],
  "Random Survival Forest" = default_colors[4],
  "GAM Cox Model" = default_colors[5],
  "Penalized Cox Model" = default_colors[6]
)

# Main plot
main_plot <- ggplot() +
  geom_line(data = c_index_df, aes(x = Split, y = c_index_deepsurv, 
                                   color = "DeepSurv Model")) +
  geom_line(data = c_index_df, aes(x = Split, y = c_index_bilstmdeepsurv, 
                                   color = "BiLSTM DeepSurv")) +
  geom_line(data = c_index_df, aes(x = Split, y = c_index_cox, 
                                   color = "Cox Model")) +
  geom_line(data = c_index_df, aes(x = Split, y = c_index_rsf, 
                                   color = "Random Survival Forest")) +
  geom_line(data = c_index_df, aes(x = Split, y = c_index_gamcox, 
                                   color = "GAM Cox Model")) +
  geom_line(data = c_index_df, aes(x = Split, y = c_index_pencox, 
                                   color = "Penalized Cox Model")) +
  
  geom_point(data = c_index_df, aes(x = Split, y = c_index_deepsurv, 
                                    color = "DeepSurv Model")) +
  geom_point(data = c_index_df, aes(x = Split, y = c_index_bilstmdeepsurv, 
                                    color = "BiLSTM DeepSurv")) +
  geom_point(data = c_index_df, aes(x = Split, y = c_index_cox, 
                                    color = "Cox Model")) +
  geom_point(data = c_index_df, aes(x = Split, y = c_index_rsf, 
                                    color = "Random Survival Forest")) +
  geom_point(data = c_index_df, aes(x = Split, y = c_index_gamcox, 
                                    color = "GAM Cox Model")) +
  geom_point(data = c_index_df, aes(x = Split, y = c_index_pencox, 
                                    color = "Penalized Cox Model")) +
  
  geom_hline(aes(yintercept = mean_c_index_deepsurv, 
                 color = "DeepSurv Model"), linetype = "dashed", size = 1) +
  geom_hline(aes(yintercept = mean_c_index_bilstm, 
                 color = "BiLSTM DeepSurv"), linetype = "dashed", size = 1) +
  geom_hline(aes(yintercept = mean_c_index_cox, 
                 color = "Cox Model"), linetype = "dashed", size = 1) +
  geom_hline(aes(yintercept = mean_c_index_rsf, 
                 color = "Random Survival Forest"), linetype = "dashed", 
             size = 1) +
  geom_hline(aes(yintercept = mean_c_index_gamcox, 
                 color = "GAM Cox Model"), linetype = "dashed", size = 1) +
  geom_hline(aes(yintercept = mean_c_index_pencox, 
                 color = "Penalized Cox Model"), linetype = "dashed", 
             size = 1) +
  
  scale_color_manual(values = model_colors) +
  labs(
    x = "Split Number",
    y = "C-index",
    color = "Model"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    legend.position = "bottom"
  )

# Title plot
title_plot <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, 
           label = "Comparison of C-index across 100 Splits", size = 6, 
           fontface = "bold") +
  theme_void()

# Bottom caption with updated means
mean_label <- paste(
  "Mean C-index | DeepSurv:", round(mean_c_index_deepsurv, 5),
  "| BiLSTM DeepSurv:", round(mean_c_index_bilstm, 5),
  "| Cox:", round(mean_c_index_cox, 5),
  "| RSF:", round(mean_c_index_rsf, 5),
  "| GAM Cox:", round(mean_c_index_gamcox, 5),
  "| Penalized Cox:", round(mean_c_index_pencox, 5)
)

bottom_caption <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = mean_label, size = 4, 
           hjust = 0.5) +
  theme_void()

# Combine all
final_plot <- title_plot / main_plot / bottom_caption + 
  plot_layout(heights = c(0.08, 1, 0.1))
print(final_plot)
@