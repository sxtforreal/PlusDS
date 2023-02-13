library(glmnet)
library(readxl)
library(tidyverse)
library(dplyr)
library(knitr)
library(ggplot2)
library(corrplot)
library(caret)
library(psych)
library(GPArotation)
library(limma)
library(Biobase)
library(survival)
library(survminer)

### Import data
covariates <- read.csv(file = '/Users/sunxiaotan/Desktop/+DS/Data+ Project/dataset/Raw data.csv')
outcomes <- read.csv(file = '/Users/sunxiaotan/Desktop/+DS/Data+ Project/dataset/shhs-cvd-summary-dataset-0.13.0.csv')
int.var <- read_xlsx('/Users/sunxiaotan/Desktop/+DS/Data+ Project/Interested variables.xlsx')

### Data curation
# Remove 'alphad1' because it contains 0 only
dat <- dat[ ,-which(names(dat) %in% c("alphad1"))]

# Create dummy variables for nominal categorical variables with more than 2 levels -- 'mstat', 'race'

# Create dummy variables for 'mstat'
# 'mstat' 1:married 2:widowed 3:divorced/seperated 4: never married
dat$married[dat$mstat == 1] = 1
dat$married[is.na(dat$married)] <- 0
dat$widowed[dat$mstat == 2] = 1
dat$widowed[is.na(dat$widowed)] <- 0
dat$seperated[dat$mstat == 3] = 1
dat$seperated[is.na(dat$seperated)] <- 0
dat$never_married[dat$mstat == 4] = 1
dat$never_married[is.na(dat$never_married)] <- 0

# Create dummy variables for 'race'
# 'race' 1:white 2:black 3:other 
dat$white[dat$race == 1] = 1
dat$white[is.na(dat$white)] <- 0
dat$black[dat$race == 2] = 1
dat$black[is.na(dat$black)] <- 0
dat$other[dat$race == 3] = 1
dat$other[is.na(dat$other)] <- 0

# Remove original 'mstat' and 'race'
dat <- dat[, -which(names(dat) %in% c("mstat","race"))]

# Replace level == '8' with NA
dat$angina15 <- factor(dat$angina15)
levels(dat$angina15)[levels(dat$angina15)=="8"]<-NA
dat$mi15 <- factor(dat$mi15)
levels(dat$mi15)[levels(dat$mi15)=="8"]<-NA
dat$stroke15 <- factor(dat$stroke15)
levels(dat$stroke15)[levels(dat$stroke15)=="8"]<-NA
dat$hf15 <- factor(dat$hf15)
levels(dat$hf15)[levels(dat$hf15)=="8"]<-NA
dat$cabg15 <- factor(dat$cabg15)
levels(dat$cabg15)[levels(dat$cabg15)=="8"]<-NA
dat$ca15 <- factor(dat$ca15)
levels(dat$ca15)[levels(dat$ca15)=="8"]<-NA
dat$othrcs15 <- factor(dat$othrcs15)
levels(dat$othrcs15)[levels(dat$othrcs15)=="8"]<-NA
dat$pacem15 <- factor(dat$pacem15)
levels(dat$pacem15)[levels(dat$pacem15)=="8"]<-NA
dat$sa15 <- factor(dat$sa15)
levels(dat$sa15)[levels(dat$sa15)=="8"]<-NA
dat$emphys15 <- factor(dat$emphys15)
levels(dat$emphys15)[levels(dat$emphys15)=="8"]<-NA
dat$crbron15 <- factor(dat$crbron15)
levels(dat$crbron15)[levels(dat$crbron15)=="8"]<-NA
dat$copd15 <- factor(dat$copd15)
levels(dat$copd15)[levels(dat$copd15)=="8"]<-NA
dat$asthma15 <- factor(dat$asthma15)
levels(dat$asthma15)[levels(dat$asthma15)=="8"]<-NA
dat$asth1215 <- factor(dat$asth1215)
levels(dat$asth1215)[levels(dat$asth1215)=="8"]<-NA
dat$cough315 <- factor(dat$cough315)
levels(dat$cough315)[levels(dat$cough315)=="8"]<-NA
dat$phlegm15 <- factor(dat$phlegm15)
levels(dat$phlegm15)[levels(dat$phlegm15)=="8"]<-NA
dat$runny15 <- factor(dat$runny15)
levels(dat$runny15)[levels(dat$runny15)=="8"]<-NA
dat$sinus15 <- factor(dat$sinus15)
levels(dat$sinus15)[levels(dat$sinus15)=="8"]<-NA
dat$nitro15 <- factor(dat$nitro15)
levels(dat$nitro15)[levels(dat$nitro15)=="8"]<-NA

# Number of missing values
na_count <- map(dat[, !names(dat) %in% c('nsrrid','any_chd','any_cvd','date_to_event')], ~sum(is.na(.)))
na_count <- data.frame(t(data.frame(na_count)))

# Percentage of missing values
na_pct <- na_count/nrow(dat)
names(na_pct) <- 'na_pct'

# Create list of 'binary variables under 40% NA rate', 'non-binary variables under 40% NA rate', and 'variables above 40% NA rate'
# Imputation methods for variables: 0%-40% replace with mean or mode, 40%-95% replace with binary indicator
na_pct$category = NA
na_pct$category[na_pct$na_pct < 0.4] = 'impute'
na_pct$category[na_pct$na_pct >= 0.4] = 'binary_indicator'
over40 <- row.names(na_pct[na_pct[,'category'] == 'binary_indicator',])

# Create level counts table
levels_count <- map(dat[, !names(dat) %in% c('nsrrid','any_chd','any_cvd','date_to_event')], ~length(table(.)))
levels_count <- data.frame(t(data.frame(levels_count)))
names(levels_count) <- 'count'
binary_below40 <- row.names(subset(levels_count, count == 2))
binary_below40 <- binary_below40[! binary_below40 %in% c('prev_hx_mi', 'prev_hx_stroke')]
other_below40 <- row.names(subset(levels_count, count > 2))
other_below40 <- other_below40[! other_below40 %in% c('ankbp', 'armbp')]

# Create get mode function
getmode <- function(x) {
   tbl <- sort(table(x))
   if(var(tbl) == 0){NA
     }else{as.numeric(names(tbl)[which(tbl==max(tbl))])}
}
# Replace missing values in each binary variable with their mode respectively
for(i in binary_below40){
  dat[is.na(dat[,i]), i] <- getmode(dat[,i])
}
# Replace missing values in each non-binary variable with their mean respectively
for(i in other_below40){
  dat[is.na(dat[,i]), i] <- mean(dat[,i], na.rm = TRUE)
  }
# Replace with binary indicators
for(i in over40){
  for(j in 1:nrow(dat)){
    if(is.na(dat[j, i]) == T){
      dat[j,i] <- 0
    } else{
      dat[j,i] <- 1
    }
  }
}

# Turn factor variables into numerical variables
dat$angina15 <- as.numeric(as.character(dat$angina15))
dat$mi15 <- as.numeric(as.character(dat$mi15))
dat$stroke15 <- as.numeric(as.character(dat$stroke15))
dat$hf15 <- as.numeric(as.character(dat$hf15))
dat$cabg15 <- as.numeric(as.character(dat$cabg15))
dat$ca15 <- as.numeric(as.character(dat$ca15))
dat$othrcs15 <- as.numeric(as.character(dat$othrcs15))
dat$pacem15 <- as.numeric(as.character(dat$pacem15))
dat$sa15 <- as.numeric(as.character(dat$sa15))
dat$emphys15 <- as.numeric(as.character(dat$emphys15))
dat$crbron15 <- as.numeric(as.character(dat$crbron15))
dat$copd15 <- as.numeric(as.character(dat$copd15))
dat$asthma15 <- as.numeric(as.character(dat$asthma15))
dat$asth1215 <- as.numeric(as.character(dat$asth1215))
dat$cough315 <- as.numeric(as.character(dat$cough315))
dat$phlegm15 <- as.numeric(as.character(dat$phlegm15))
dat$runny15 <- as.numeric(as.character(dat$runny15))
dat$sinus15 <- as.numeric(as.character(dat$sinus15))
dat$nitro15 <- as.numeric(as.character(dat$nitro15))

# Check correlation matrix
correlations_before <- cor(dat[, !names(dat) %in% c('nsrrid','any_chd','any_cvd','date_to_event')])

# Remove one variable from the pair with a correlation of 0.9 or above. This aims to solve co-linearity, prevent correlation matrix become singular
suggest_removal <- findCorrelation(correlations_before, cutoff = 0.9, verbose = FALSE, names = TRUE, exact = ncol(correlations_before) < 100)
dat <- dat[, -which(names(dat) %in% suggest_removal)]

# Re-order columns
dat <- dat %>% relocate(any_chd, .after = last_col())
dat <- dat %>% relocate(any_cvd, .after = last_col())
dat <- dat %>% relocate(date_to_event, .after = last_col())

### Exploratory Data Analysis
#Quantify the univariate associations between features and outcomes

# chd
design_chd <- model.matrix(~ any_chd, data = dat[,2:86])
exp <- t(as.matrix(dat[,2:85]))
fit_chd <- lmFit(exp, design_chd)
fit_chd <- eBayes(fit_chd)
tT_chd <- topTable(fit_chd, number = 84, adjust.method = "BH")
tT_chd <- subset(tT_chd, select = c("P.Value", "adj.P.Val", "B"))
chd_sig <- tT_chd[tT_chd$adj.P.Val < 0.05,]
# Limma suggests 58 features that have significant univariate association with any_chd
print(rownames(chd_sig))

design_cvd <- model.matrix(~ any_cvd, data = dat[,c(2:85,87)])
fit_cvd <- lmFit(exp, design_cvd)
fit_cvd <- eBayes(fit_cvd)
tT_cvd <- topTable(fit_cvd, number = 84, adjust.method = "BH")
tT_cvd <- subset(tT_cvd, select = c("P.Value", "adj.P.Val", "B"))
cvd_sig <- tT_cvd[tT_cvd$adj.P.Val < 0.05,]
# Limma suggests 59 features that have significant univariate association with any_cvd
print(rownames(cvd_sig))

# Dimension reduction by factor analysis

# Round correlation matrix to 2 decimal places in order to solve singularity
cor_dat <- round(cor(dat[, !names(dat) %in% c('nsrrid','any_chd','any_cvd','date_to_event')]), 2)
corrplot(cor_dat)

# Parallel Analysis
# Extract important eigenvalues through comparison with a random matrix with same dimension. If the eigenvalue from our data is greater than the average eigenvalues from a random matrix, we keep this eigenvalue and mark it as significant. Note, eigenvalues correspond to PCs.
fa.parallel(cor_dat, n.obs = 5042, fa = "fa", n.iter = 100, main = "Parallel Analysis Scree Plot")
# Parallel analysis suggests we use 28 factors to represent the whole dataset

# Extract factors
Factor.Result <- fa(cor_dat, nfactors = 28, rotate = "promax", fm = "pa")
# PA -- Component loading
Factor.Result

# Create a table for variables
# For each row, includes the name of the variable, its optimal factor selection, the corresponding value in that factor and the corresponding absolute value
PA_dat <- as.data.frame(colnames(dat[, !names(dat) %in% c('nsrrid','any_chd','any_cvd','date_to_event')]))
PA_dat$opt_factor <- NA
PA_dat$Values <- NA
names(PA_dat) <- c('Variables', 'Opt_factor', 'Values')

for(i in 1:nrow(PA_dat)){
  a <- max(abs(Factor.Result[["loadings"]][i,]))
  b <- which(abs(Factor.Result[["loadings"]][i,]) == a)
  PA_dat[i,2] <- b
  PA_dat[i,3] <- Factor.Result[["loadings"]][i,][b]
}
PA_dat$abs_values <- abs(PA_dat$Values)

# Check number of components in each factor
table(PA_dat$Opt_factor)

# Create a table for factors and their components
PA_Comp <- data.frame(PA = paste('PA',1:28, sep = ''))
PA_Comp$Components <- NA
for(i in 1:28){
  PA_Comp[i,2] <- toString(paste(PA_dat[PA_dat$Opt_factor == i,]$Variables, sep = ','))
}

# Compute score matrix (Score matrix = X matrix * Loadings matrix)
loadings <- as.matrix(Factor.Result$loadings)
score_matrix <- as.matrix(dat[, !names(dat) %in% c('nsrrid','any_chd','any_cvd','date_to_event')])%*%loadings
score_matrix <- as.data.frame(score_matrix)

# Join two datasets, combine 28 PAs with 2 outcomes
Reg_dat <- cbind(score_matrix, oc)

### Lasso - Cox Proportion Hazards Model

# Remove rows with 0 date_to_event
# Split data into training set and validation
cox_dat <- dat %>% filter(date_to_event != 0)
cox_train <- sample(1:nrow(cox_dat), nrow(cox_dat)*0.75)
cox_train_dat <- cox_dat[cox_train,]
cox_test <- (-cox_train)
cox_test_dat <- cox_dat[cox_test,]
cox_exp_train <- cox_train_dat[,2:85]
cox_exp_test <- cox_test_dat[,2:85]
set.seed(1)

# CHD
cox_y_chd_train <- Surv(time = cox_train_dat$date_to_event, event = cox_train_dat$any_chd)
cox_y_chd_test <- Surv(time = cox_test_dat$date_to_event, event = cox_test_dat$any_chd)
# Cross-validation with c-index
cvfit_chd <- cv.glmnet(data.matrix(cox_exp_train), cox_y_chd_train, family = "cox", type.measure = "C", alpha = 1)
coef(cvfit_chd, cvfit_chd$lambda.min)
a <- do.call(rbind.data.frame, as.list(coef(cvfit_chd, cvfit_chd$lambda.min)))
# Predict relative risk of chd
pred_chd <- predict(cvfit_chd, newx = as.matrix(cox_exp_test), s = "lambda.min", type = "response")

# C-index
Cindex(pred_chd, cox_y_chd_test)

# CVD
cox_y_cvd_train <- Surv(time = cox_train_dat$date_to_event, event = cox_train_dat$any_cvd)
cox_y_cvd_test <- Surv(time = cox_test_dat$date_to_event, event = cox_test_dat$any_cvd)
# Cross-validation with c-index
cvfit_cvd <- cv.glmnet(data.matrix(cox_exp_train), cox_y_cvd_train, family = "cox", type.measure = "C", alpha = 1)
coef(cvfit_cvd, cvfit_cvd$lambda.min)
b <- do.call(rbind.data.frame, as.list(coef(cvfit_cvd, cvfit_cvd$lambda.min)))
# Predict relative risk of cvd
pred_cvd <- predict(cvfit_cvd, newx = as.matrix(cox_exp_test), s = "lambda.min", type = "response")

# C-index
Cindex(pred_cvd, cox_y_cvd_test)

### XGBoost-AFT model
library(xgboost)

# Associate ranged labels with the data matrix(chd)
for (i in 1:nrow(outcomes)) {
  if (outcomes$any_chd[i] == 0){
    outcomes$chd_lower_bound[i] <- outcomes$date_to_event[i]
    outcomes$chd_upper_bound[i] <- +Inf
  }
  if (outcomes$any_chd[i] == 1){
    outcomes$chd_lower_bound[i] <- outcomes$date_to_event[i]
    outcomes$chd_upper_bound[i] <- outcomes$date_to_event[i]
  }
}

# Associate ranged labels with the data matrix(cvd)
for (i in 1:nrow(outcomes)) {
  if (outcomes$any_cvd[i] == 0){
    outcomes$cvd_lower_bound[i] <- outcomes$date_to_event[i]
    outcomes$cvd_upper_bound[i] <- +Inf
  }
  if (outcomes$any_cvd[i] == 1){
    outcomes$cvd_lower_bound[i] <- outcomes$date_to_event[i]
    outcomes$cvd_upper_bound[i] <- outcomes$date_to_event[i]
  }
}

dat <- inner_join(dat, outcomes[,c("nsrrid", "chd_lower_bound", "chd_upper_bound", "cvd_lower_bound", "cvd_upper_bound")], by = 'nsrrid')

training <- sample(1:nrow(dat), nrow(dat)*0.75)
training_dat <- dat[training,]
testing <- (-training)
testing_dat <- dat[testing,]

X_train <- as.matrix(training_dat[,2:85])
X_test <- as.matrix(testing_dat[,2:85])
dtrain <- xgb.DMatrix(X_train)
dtest <- xgb.DMatrix(X_test)

# CHD
# Training
y_lower_bound_train_chd <- training_dat$chd_lower_bound
y_upper_bound_train_chd <- training_dat$chd_upper_bound
setinfo(dtrain, 'label_lower_bound', y_lower_bound_train_chd)
setinfo(dtrain, 'label_upper_bound', y_upper_bound_train_chd)

y_lower_bound_test_chd <- testing_dat$chd_lower_bound
y_upper_bound_test_chd <- testing_dat$chd_upper_bound
setinfo(dtest, 'label_lower_bound', y_lower_bound_test_chd)
setinfo(dtest, 'label_upper_bound', y_upper_bound_test_chd)

params <- list(objective='survival:aft',
               eval_metric='aft-nloglik',
               aft_loss_distribution='normal',
               aft_loss_distribution_scale=1.20,
               tree_method='hist',
               learning_rate=0.3,
               max_depth=4)
watchlist <- list(train = dtrain, test = dtest)
bst <- xgb.train(params, dtrain, nrounds=5, watchlist)


# View feature importance
importance_matrix <- xgb.importance(model = bst)
print(importance_matrix)
xgb.plot.importance(importance_matrix = importance_matrix)

# View the trees from a model
xgb.dump(bst, with_stats = TRUE)
xgb.plot.tree(model = bst)

# CVD
# Training
y_lower_bound_train <- training_dat$cvd_lower_bound
y_upper_bound_train <- training_dat$cvd_upper_bound
setinfo(dtrain, 'label_lower_bound', y_lower_bound_train)
setinfo(dtrain, 'label_upper_bound', y_upper_bound_train)

y_lower_bound_test <- testing_dat$cvd_lower_bound
y_upper_bound_test <- testing_dat$cvd_upper_bound
setinfo(dtest, 'label_lower_bound', y_lower_bound_test)
setinfo(dtest, 'label_upper_bound', y_upper_bound_test)

params <- list(objective='survival:aft',
               eval_metric='aft-nloglik',
               aft_loss_distribution='normal',
               aft_loss_distribution_scale=1.20,
               tree_method='hist',
               learning_rate=0.3,
               max_depth=4)
watchlist <- list(train = dtrain, test = dtest)
bst <- xgb.train(params, dtrain, nrounds=5, watchlist)


# View feature importance
importance_matrix <- xgb.importance(model = bst)
print(importance_matrix)
xgb.plot.importance(importance_matrix = importance_matrix)

# View the trees from a model
xgb.dump(bst, with_stats = TRUE)
xgb.plot.tree(model = bst)
