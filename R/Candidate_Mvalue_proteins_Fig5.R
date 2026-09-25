#### All cohort M-value prediction ####
library("tidyverse")

Discovery_proteome_data <- readRDS(file = "data/processed_datasets.rds")
Discovery_clinical_data <- readRDS(file = "data/Clinical_data_filt")
Discovery_sample_data <- readRDS(file = "data/Sample_ID_filt.rds")
Discovery_proteome_data <- Discovery_proteome_data$Olink_discovery_cohort
rownames(Discovery_proteome_data) <- Discovery_proteome_data$ID
Discovery_proteome_data <- Discovery_proteome_data[,-1]

Discovery_proteome_data <- Discovery_proteome_data[,Discovery_sample_data$Clamp == "Pre"]
Discovery_clinical_data <- Discovery_clinical_data[Discovery_sample_data$Clamp == "Pre",]
Discovery_proteome_data <- Discovery_proteome_data[,-c(14)]
Discovery_clinical_data <- Discovery_clinical_data[-c(14),]

length(which(Discovery_clinical_data$Group == "T2D" & Discovery_clinical_data$Gender == "Male"))
length(which(Discovery_clinical_data$Group == "T2D" & Discovery_clinical_data$Gender == "Female"))
length(which(Discovery_clinical_data$Group == "NGT" & Discovery_clinical_data$Gender == "Male"))
length(which(Discovery_clinical_data$Group == "NGT" & Discovery_clinical_data$Gender == "Female"))

length(which(Discovery_clinical_data$Gender == "Male"))
sd(Discovery_clinical_data$M..µmol..kg.min..)/sqrt(length(Discovery_clinical_data$M..µmol..kg.min..))
mean(Discovery_clinical_data$M..µmol..kg.min..)
median(Discovery_clinical_data$M..µmol..kg.min..)

HIIT_proteome_data <- readRDS(file = "data/processed_datasets_validation.rds")
HIIT_clinical_data <- readRDS(file = "data/Val_Clinical_data_Olink.rds")
HIIT_proteome_data <- HIIT_proteome_data$Olink_validation_cohort
HIIT_proteome_data <- HIIT_proteome_data[,HIIT_clinical_data$Condition == "Pre"]

HIIT_clinical_data <- HIIT_clinical_data[HIIT_clinical_data$Condition == "Pre",]
HIIT_clinical_data$Gender <- "Male"

length(which(HIIT_clinical_data$Group == "3"))
length(which(HIIT_clinical_data$Group == "2"))
length(which(HIIT_clinical_data$Group == "1"))

sd(HIIT_clinical_data$M.value)/sqrt(length(HIIT_clinical_data$M.value))
mean(HIIT_clinical_data$M.value)
median(HIIT_clinical_data$M.value)



IDA_proteome_data <- readRDS(file = "Olink_proteome_validation_cohort.rds")
IDA_proteome_data <- IDA_proteome_data[, !colnames(IDA_proteome_data) %in% c("PARTICIPANT_V1", "PARTICIPANT_V2", "PARTICIPANT_V3")]
IDA_clinical_data <- readRDS(file = "Clinical_data_validation_cohort.rds")
IDA_clinical_data <- IDA_clinical_data[!IDA_clinical_data$ID %in% c("PARTICIPANT_V1","PARTICIPANT_V2"),]
IDA_clinical_data$Gender <- as.character(IDA_clinical_data$sex)

IDA_additional_clinical_data <- read.delim("data-raw/Clinical_data_validation_cohort.txt", dec = ",")
IDA_additional_clinical_data <- IDA_additional_clinical_data[!IDA_additional_clinical_data$PatientID %in% c("PARTICIPANT_V1","PARTICIPANT_V2"),]
IDA_additional_clinical_data$Gender <- as.character(IDA_additional_clinical_data$Gender)
IDA_additional_clinical_data <- IDA_additional_clinical_data[IDA_additional_clinical_data$PatientID %in% IDA_clinical_data$ID,]
IDA_additional_clinical_data <- IDA_additional_clinical_data[
  match(IDA_clinical_data$ID, IDA_additional_clinical_data$PatientID),
]
IDA_additional_clinical_data$ID <- IDA_additional_clinical_data$PatientID
IDA_combined_clinical_data <- inner_join(IDA_clinical_data[,c(3,12)],IDA_additional_clinical_data, by = "ID")
IDA_combined_clinical_data <- IDA_combined_clinical_data[,!colnames(IDA_combined_clinical_data) %in% c("PatientID","X")]

IDA_Lipid_data <- read.delim("data-raw/Clinical_data_validation_cohort_lipids.txt",dec=",")
IDA_Lipid_data <- IDA_Lipid_data[!is.na(IDA_Lipid_data$HDL),]
IDA_combined_clinical_data <- cbind(IDA_combined_clinical_data,IDA_Lipid_data[match(IDA_combined_clinical_data$ID,IDA_Lipid_data$PatientID),c(1,8:17)])
IDA_combined_clinical_data<- IDA_combined_clinical_data[!is.na(IDA_combined_clinical_data$PatientID),]

length(common_proteins)
# 1. Identify common proteins (rownames must match exactly)
common_proteins <- Reduce(intersect, list(
  rownames(Discovery_proteome_data),
  rownames(HIIT_proteome_data),
  rownames(IDA_proteome_data)
))

# 2. Subset each dataset to the shared proteins
Discovery_sub <- Discovery_proteome_data[common_proteins, ]
HIIT_sub      <- HIIT_proteome_data[common_proteins, ]
IDA_sub       <- IDA_proteome_data[common_proteins, ]

# 3. Combine all datasets column-wise
proteome_combined <- cbind(Discovery_sub, HIIT_sub, IDA_sub)

library(limma)

# 1. Create a batch vector corresponding to each cohort
batch <- c(
  rep("Discovery", ncol(Discovery_sub)),
  rep("HIIT", ncol(HIIT_sub)),
  rep("IDA", ncol(IDA_sub))
)

# 2. Apply batch correction
proteome_corrected <- removeBatchEffect(proteome_combined, batch = batch)

#### Clinical data ####
m_values <- c(
  Discovery_clinical_data$`M..µmol..kg.min..`,
  HIIT_clinical_data$M.value,
  IDA_combined_clinical_data$M.value
)

Gender <- c(
  Discovery_clinical_data$Gender,
  HIIT_clinical_data$Gender,
  IDA_combined_clinical_data$Gender
)

Age <- c(
  Discovery_clinical_data$Age,
  HIIT_clinical_data$Age,
  IDA_combined_clinical_data$Age
)

BMI <- c(
  Discovery_clinical_data$BMI..kg.m2.,
  HIIT_clinical_data$BMI1,
  IDA_combined_clinical_data$BMI
)

HbA1c <- c(
  Discovery_clinical_data$HbA1c..mmol.mol.,
  HIIT_clinical_data$HbA1c1,
  IDA_combined_clinical_data$HbA1c
)

FM <- c(
  Discovery_clinical_data$Body.fat....,
  HIIT_clinical_data$Fat_per1,
  IDA_combined_clinical_data$Fat.mass
)

FPG <- c(
  Discovery_clinical_data$FP.Glucose..mmol.L.,
  HIIT_clinical_data$F_gluc1,
  IDA_combined_clinical_data$Fast.gluco
)

TG <- c(
    Discovery_clinical_data$P.TG..mmol.L.,
    HIIT_clinical_data$TG1,
    IDA_combined_clinical_data$Triglyc
)

HDL <- c(
    Discovery_clinical_data$HDL.Chol..mmol.L.,
    HIIT_clinical_data$HDL1,
    IDA_combined_clinical_data$HDL
)

TG_HDL_ratio <- c(
    Discovery_clinical_data$P.TG..mmol.L./Discovery_clinical_data$HDL.Chol..mmol.L.,
    HIIT_clinical_data$TG1/HIIT_clinical_data$HDL1,
    IDA_combined_clinical_data$Triglyc/IDA_combined_clinical_data$HDL
)

M_value_z_by_cohort <- c(
    scale(Discovery_clinical_data$M..µmol..kg.min..),
    scale(HIIT_clinical_data$M.value),
    scale(IDA_combined_clinical_data$M.value)
)

# Optional: Create a corresponding cohort label
cohort <- c(
  rep("Discovery", length(Discovery_clinical_data$`M..µmol..kg.min..`)),
  rep("HIIT", length(HIIT_clinical_data$M.value)),
  rep("IDA", length(IDA_combined_clinical_data$M.value))
)

# Create a dataframe
combined_clinical <- data.frame(
  M_value = m_values,
  Cohort = cohort,
  Gender = Gender,
  Age = Age,
  BMI = BMI,
  HbA1c = HbA1c,
  TG = TG,
  HDL = HDL,
  TG_HDL_ratio,
  M_value_z_by_cohort
)


# Example M-value vectors from each cohort
m_discovery <- Discovery_clinical_data$`M..µmol..kg.min..`
m_hiit <- HIIT_clinical_data$M.value
m_ida <- IDA_combined_clinical_data$M.value

# Combine M-values and cohort labels
combined_df <- data.frame(
  M_value = c(m_discovery, m_hiit, m_ida),
  Cohort = rep(c("Discovery", "HIIT", "IDA"),
               times = c(length(m_discovery), length(m_hiit), length(m_ida)))
)

# Create rank variable
combined_df$Rank <- rank(combined_df$M_value, ties.method = "first")

# Plot
library(ggplot2)
ggplot(combined_df, aes(x = reorder(factor(Rank), Rank), y = M_value, fill = Cohort)) +
  geom_bar(stat = "identity", width = 1, alpha=0.9) +
  scale_fill_manual(values = c("Discovery" = "darkgrey", "HIIT" = "darkgrey", "IDA" = "#4564a3")) +
  labs(x = "Ranked Individuals", y = "M-value (μmol/kg/min)", title = "Ranked M-values by Cohort") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "top"
  )



library(ppcor)

# Initialize result holder
partial_results <- data.frame(Protein = rownames(proteome_corrected),
                              rho = NA, pval = NA, adj_pval = NA)

# Loop through each protein
for (i in seq_len(nrow(proteome_corrected))) {
  protein_values <- as.numeric(proteome_corrected[i, ])
  m_values <- combined_clinical$M_value_z_by_cohort
  gender <- factor(combined_clinical$Gender)
  Age <- as.numeric(Age)
  BMI <- as.numeric(BMI)


  # Remove missing values
  complete_idx <- complete.cases(protein_values, m_values, gender)

  if (sum(complete_idx) >= 4) {  # Require at least 4 data points
    test <- pcor.test(
      x = protein_values[complete_idx],
      y = m_values[complete_idx],
      z = data.frame(Gender = as.numeric(gender[complete_idx]),Age, BMI),
      method = "spearman"
    )

    partial_results$rho[i] <- test$estimate
    partial_results$pval[i] <- test$p.value
  }
}

# Adjust p-values for multiple testing
partial_results$adj_pval <- p.adjust(partial_results$pval, method = "BH")

#### M-value prediction ####
# Filter by Olink 3K UKBB
Olink_3K_UKBB <- read.delim("data-raw/olink_assay_1013.txt")
proteome_corrected_UKBB <- proteome_corrected[word(rownames(proteome_corrected),1,sep="_") %in% Olink_3K_UKBB$Assay,]

# Identify significant proteins
significant_proteins <- partial_results$Protein[partial_results$adj_pval < 0.05]

# Subset the proteomics data matrix to only those proteins
proteome_significant <- proteome_corrected_UKBB[rownames(proteome_corrected_UKBB) %in% significant_proteins, ]

proteome_scaled <- t(scale(t(proteome_significant)))
rownames(combined_clinical) <- colnames(proteome_scaled)






















#### Mean-centered Mvalue ####

#### LOCO IDA cohort M-val prediction

# Safe LOCO forward selection
loco_forward_selection <- function(clinical_data, protein_data, target = "M_value_z_by_cohort",
                                   method = c("rf", "lasso", "svr", "knn"),
                                   max_features = 20, cv_folds = 5, seed = 123) {

    set.seed(seed)
    method <- match.arg(method)
    clinical_data$Gender <- ifelse(clinical_data$Gender == "Male", 1, 0)

    common_ids <- intersect(colnames(protein_data), rownames(clinical_data))
    protein_data <- protein_data[, common_ids]
    clinical_data <- clinical_data[common_ids, ]

    cohorts <- unique(clinical_data$Cohort)
    results <- list()

    for (test_cohort in cohorts) {
        message(sprintf("Testing on %s", test_cohort))

        train_ids <- rownames(clinical_data)[clinical_data$Cohort != test_cohort]
        test_ids  <- rownames(clinical_data)[clinical_data$Cohort == test_cohort]

        y_train <- clinical_data[train_ids, target]
        y_test  <- clinical_data[test_ids, target]

        clinical_vars <- clinical_data[, c("Age", "Gender", "BMI", "HbA1c","TG_HDL_ratio")]
        pre_proc <- caret::preProcess(clinical_vars, method = c("center", "scale"))
        X_scaled <- predict(pre_proc, clinical_vars)

        X_train <- X_scaled[train_ids, , drop = FALSE]
        X_test  <- X_scaled[test_ids, , drop = FALSE]

        selected <- c()
        remaining <- rownames(protein_data)
        r2_progress <- c()

        for (step in 1:max_features) {
            best_cv_r2 <- -Inf
            best_prot <- NULL

            for (prot in setdiff(remaining, selected)) {
                temp_train <- X_train
                temp_train[[prot]] <- as.numeric(protein_data[prot, train_ids])

                folds <- caret::createFolds(y_train, k = cv_folds, returnTrain = TRUE)
                fold_r2s <- c()

                for (fold in folds) {
                    X_tr <- temp_train[fold, , drop = FALSE]
                    y_tr <- y_train[fold]
                    X_val <- temp_train[-fold, , drop = FALSE]
                    y_val <- y_train[-fold]

                    pred <- switch(method,
                                   rf = {
                                       fit <- ranger::ranger(y = y_tr, x = X_tr)
                                       predict(fit, data = X_val)$predictions
                                   },
                                   lasso = {
                                       fit <- glmnet::cv.glmnet(as.matrix(X_tr), y_tr, alpha = 1)
                                       predict(fit, newx = as.matrix(X_val), s = "lambda.min")
                                   },
                                   svr = {
                                       fit <- svm(x = X_tr, y = y_tr, kernel = "radial")
                                       predict(fit, newdata = X_val)
                                   },
                                   knn = {
                                       FNN::knn.reg(train = X_tr, test = X_val, y = y_tr, k = 5)$pred
                                   },
                                   xgb = {
                                       dtrain <- xgb.DMatrix(data = as.matrix(X_tr), label = y_tr)
                                       dval <- xgb.DMatrix(data = as.matrix(X_val))
                                       model <- xgboost(data = dtrain, nrounds = 100, objective = "reg:squarederror", verbose = 0)
                                       predict(model, dval)
                                   }
                    )

                    fold_r2s <- c(fold_r2s, r2_score(y_val, pred))
                }

                mean_r2 <- mean(fold_r2s, na.rm = TRUE)
                if (mean_r2 > best_cv_r2) {
                    best_cv_r2 <- mean_r2
                    best_prot <- prot
                }
            }

            selected <- c(selected, best_prot)
            r2_progress <- c(r2_progress, best_cv_r2)
            X_train[[best_prot]] <- as.numeric(protein_data[best_prot, train_ids])
            X_test[[best_prot]]  <- as.numeric(protein_data[best_prot, test_ids])
        }

        best_k <- which.max(r2_progress)
        final_prots <- selected[1:best_k]

        final_train <- X_scaled[train_ids, , drop = FALSE]
        final_test  <- X_scaled[test_ids, , drop = FALSE]
        for (prot in final_prots) {
            final_train[[prot]] <- as.numeric(protein_data[prot, train_ids])
            final_test[[prot]]  <- as.numeric(protein_data[prot, test_ids])
        }

        final_pred <- switch(method,
                             rf = {
                                 fit <- ranger::ranger(y = y_train, x = final_train)
                                 predict(fit, data = final_test)$predictions
                             },
                             lasso = {
                                 fit <- glmnet::cv.glmnet(as.matrix(final_train), y_train, alpha = 1)
                                 predict(fit, newx = as.matrix(final_test), s = "lambda.min")
                             },
                             svr = {
                                 fit <- svm(x = final_train, y = y_train, kernel = "radial")
                                 predict(fit, newdata = final_test)
                             },
                             knn = {
                                 FNN::knn.reg(train = final_train, test = final_test, y = y_train, k = 5)$pred
                             },
                             xgb = {
                                 dtrain <- xgb.DMatrix(data = as.matrix(final_train), label = y_train)
                                 dtest <- xgb.DMatrix(data = as.matrix(final_test))
                                 model <- xgboost(data = dtrain, nrounds = 100, objective = "reg:squarederror", verbose = 0)
                                 predict(model, dtest)
                             }
        )

        scores <- correlation_scores(y_test, final_pred)

        results[[test_cohort]] <- list(
            r2 = r2_score(y_test, final_pred),
            pearson2 = scores$pearson2,
            spearman2 = scores$spearman2,
            predicted = final_pred,
            observed = y_test,
            selected_proteins = final_prots
        )

        message(sprintf("Cohort %s — Best R²: %.4f using %d features", test_cohort, r2_score(y_test, final_pred), best_k))
    }

    return(results)
}

result_rf_z <- loco_forward_selection(
    clinical_data = combined_clinical,
    protein_data = proteome_significant,
    method = "rf"
)

result_lasso_z <- loco_forward_selection(
    clinical_data = combined_clinical,
    protein_data = proteome_significant,
    method = "lasso"
)

result_svr_z <- loco_forward_selection(
    clinical_data = combined_clinical,
    protein_data = proteome_significant,
    method = "svr"
)

result_knn_z <- loco_forward_selection(
    clinical_data = combined_clinical,
    protein_data = proteome_significant,
    method = "knn"
)

save(result_rf_z, result_lasso_z, result_svr_z, result_knn_z,
     file = "data/ML_Mvalue_prediction_data/loco_forward_selection_results.RData")



#### Protein only IDA cohort prediction ####
# Safe LOCO forward selection
loco_forward_selection_proteinONLY <- function(clinical_data, protein_data, target = "M_value_z_by_cohort",
                                               method = c("rf", "lasso", "svr", "knn"),
                                               max_features = 20, cv_folds = 5, seed = 123) {

    set.seed(seed)
    method <- match.arg(method)
    clinical_data$Gender <- ifelse(clinical_data$Gender == "Male", 1, 0)

    common_ids <- intersect(colnames(protein_data), rownames(clinical_data))
    protein_data <- protein_data[, common_ids]
    clinical_data <- clinical_data[common_ids, ]

    cohorts <- unique(clinical_data$Cohort)
    results <- list()

    for (test_cohort in cohorts) {
        message(sprintf("Testing on %s", test_cohort))

        train_ids <- rownames(clinical_data)[clinical_data$Cohort != test_cohort]
        test_ids  <- rownames(clinical_data)[clinical_data$Cohort == test_cohort]

        y_train <- clinical_data[train_ids, target]
        y_test  <- clinical_data[test_ids, target]

        # No clinical features: start with empty data frame
        X_train <- data.frame(row.names = train_ids)
        X_test  <- data.frame(row.names = test_ids)

        selected <- c()
        remaining <- rownames(protein_data)
        r2_progress <- c()

        for (step in 1:max_features) {
            best_cv_r2 <- -Inf
            best_prot <- NULL

            for (prot in setdiff(remaining, selected)) {
                temp_train <- X_train
                temp_train[[prot]] <- as.numeric(protein_data[prot, train_ids])

                folds <- caret::createFolds(y_train, k = cv_folds, returnTrain = TRUE)
                fold_r2s <- c()

                for (fold in folds) {
                    X_tr <- temp_train[fold, , drop = FALSE]
                    y_tr <- y_train[fold]
                    X_val <- temp_train[-fold, , drop = FALSE]
                    y_val <- y_train[-fold]

                    pred <- switch(method,
                                   rf = {
                                       fit <- ranger::ranger(y = y_tr, x = X_tr)
                                       predict(fit, data = X_val)$predictions
                                   },
                                   lasso = {
                                       X_tr_mat <- as.matrix(X_tr)
                                       X_val_mat <- as.matrix(X_val)
                                       if (ncol(X_tr_mat) < 2) {
                                           X_tr_mat <- cbind(X_tr_mat, dummy = 0)
                                           X_val_mat <- cbind(X_val_mat, dummy = 0)
                                       }
                                       fit <- glmnet::cv.glmnet(X_tr_mat, y_tr, alpha = 1)
                                       predict(fit, newx = X_val_mat, s = "lambda.min")
                                   },
                                   svr = {
                                       fit <- svm(x = X_tr, y = y_tr, kernel = "radial")
                                       predict(fit, newdata = X_val)
                                   },
                                   knn = {
                                       FNN::knn.reg(train = X_tr, test = X_val, y = y_tr, k = 5)$pred
                                   },
                                   xgb = {
                                       dtrain <- xgb.DMatrix(data = as.matrix(X_tr), label = y_tr)
                                       dval <- xgb.DMatrix(data = as.matrix(X_val))
                                       model <- xgboost(data = dtrain, nrounds = 100, objective = "reg:squarederror", verbose = 0)
                                       predict(model, dval)
                                   }
                    )

                    fold_r2s <- c(fold_r2s, r2_score(y_val, pred))
                }

                mean_r2 <- mean(fold_r2s, na.rm = TRUE)
                if (mean_r2 > best_cv_r2) {
                    best_cv_r2 <- mean_r2
                    best_prot <- prot
                }
            }

            selected <- c(selected, best_prot)
            r2_progress <- c(r2_progress, best_cv_r2)
            X_train[[best_prot]] <- as.numeric(protein_data[best_prot, train_ids])
            X_test[[best_prot]]  <- as.numeric(protein_data[best_prot, test_ids])
        }

        best_k <- which.max(r2_progress)
        final_prots <- selected[1:best_k]

        final_train <- data.frame(row.names = train_ids)
        final_test  <- data.frame(row.names = test_ids)
        for (prot in final_prots) {
            final_train[[prot]] <- as.numeric(protein_data[prot, train_ids])
            final_test[[prot]]  <- as.numeric(protein_data[prot, test_ids])
        }

        final_pred <- switch(method,
                             rf = {
                                 fit <- ranger::ranger(y = y_train, x = final_train)
                                 predict(fit, data = final_test)$predictions
                             },
                             lasso = {
                                 final_train_mat <- as.matrix(final_train)
                                 final_test_mat <- as.matrix(final_test)
                                 if (ncol(final_train_mat) < 2) {
                                     final_train_mat <- cbind(final_train_mat, dummy = 0)
                                     final_test_mat <- cbind(final_test_mat, dummy = 0)
                                 }
                                 fit <- glmnet::cv.glmnet(final_train_mat, y_train, alpha = 1)
                                 predict(fit, newx = final_test_mat, s = "lambda.min")
                             },
                             svr = {
                                 fit <- svm(x = final_train, y = y_train, kernel = "radial")
                                 predict(fit, newdata = final_test)
                             },
                             knn = {
                                 FNN::knn.reg(train = final_train, test = final_test, y = y_train, k = 5)$pred
                             },
                             xgb = {
                                 dtrain <- xgb.DMatrix(data = as.matrix(final_train), label = y_train)
                                 dtest <- xgb.DMatrix(data = as.matrix(final_test))
                                 model <- xgboost(data = dtrain, nrounds = 100, objective = "reg:squarederror", verbose = 0)
                                 predict(model, dtest)
                             }
        )

        scores <- correlation_scores(y_test, final_pred)

        results[[test_cohort]] <- list(
            r2 = r2_score(y_test, final_pred),
            pearson2 = scores$pearson2,
            spearman2 = scores$spearman2,
            predicted = final_pred,
            observed = y_test,
            selected_proteins = final_prots
        )

        message(sprintf("Cohort %s — Best R²: %.4f using %d features", test_cohort, r2_score(y_test, final_pred), best_k))
    }

    return(results)
}


result_rf_LOCO_protein_Z <- loco_forward_selection_proteinONLY(
    clinical_data = combined_clinical,
    protein_data = proteome_significant,
    method = "rf"
)

result_lasso_LOCO_protein_Z <- loco_forward_selection_proteinONLY(
    clinical_data = combined_clinical,
    protein_data = proteome_significant,
    method = "lasso"
)

result_svr_LOCO_protein_Z <- loco_forward_selection_proteinONLY(
    clinical_data = combined_clinical,
    protein_data = proteome_significant,
    method = "svr"
)

result_knn_LOCO_protein_Z <- loco_forward_selection_proteinONLY(
    clinical_data = combined_clinical,
    protein_data = proteome_significant,
    method = "knn"
)


save(result_rf_LOCO_protein_Z, result_lasso_LOCO_protein_Z, result_svr_LOCO_protein_Z, result_knn_LOCO_protein_Z,
     file = "data/ML_Mvalue_prediction_data/loco_forward_selection_results_PROTEINs.RData")


#### Clinical only IDA cohort prediction ####

loco_forward_selection_clinical_only <- function(clinical_data, protein_data, target = "M_value_z_by_cohort",
                                                 method = c("rf", "lasso", "svr", "knn"),
                                                 cv_folds = 5, seed = 123) {

    set.seed(seed)
    method <- match.arg(method)
    clinical_data$Gender <- ifelse(clinical_data$Gender == "Male", 1, 0)

    common_ids <- intersect(colnames(protein_data), rownames(clinical_data))
    clinical_data <- clinical_data[common_ids, ]

    cohorts <- unique(clinical_data$Cohort)
    results <- list()

    for (test_cohort in cohorts) {
        message(sprintf("Testing on %s", test_cohort))

        train_ids <- rownames(clinical_data)[clinical_data$Cohort != test_cohort]
        test_ids  <- rownames(clinical_data)[clinical_data$Cohort == test_cohort]

        y_train <- clinical_data[train_ids, target]
        y_test  <- clinical_data[test_ids, target]

        clinical_vars <- clinical_data[, c("Age", "Gender", "BMI", "HbA1c", "TG_HDL_ratio")]
        pre_proc <- caret::preProcess(clinical_vars, method = c("center", "scale"))
        X_scaled <- predict(pre_proc, clinical_vars)

        X_train <- X_scaled[train_ids, , drop = FALSE]
        X_test  <- X_scaled[test_ids, , drop = FALSE]

        final_pred <- switch(method,
                             rf = {
                                 fit <- ranger::ranger(y = y_train, x = X_train)
                                 predict(fit, data = X_test)$predictions
                             },
                             lasso = {
                                 fit <- glmnet::cv.glmnet(as.matrix(X_train), y_train, alpha = 1)
                                 predict(fit, newx = as.matrix(X_test), s = "lambda.min")
                             },
                             svr = {
                                 fit <- svm(x = X_train, y = y_train, kernel = "radial")
                                 predict(fit, newdata = X_test)
                             },
                             knn = {
                                 FNN::knn.reg(train = X_train, test = X_test, y = y_train, k = 5)$pred
                             },
                             xgb = {
                                 dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
                                 dtest <- xgb.DMatrix(data = as.matrix(X_test))
                                 model <- xgboost(data = dtrain, nrounds = 100, objective = "reg:squarederror", verbose = 0)
                                 predict(model, dtest)
                             }
        )

        scores <- correlation_scores(y_test, final_pred)

        results[[test_cohort]] <- list(
            r2 = r2_score(y_test, final_pred),
            pearson2 = scores$pearson2,
            spearman2 = scores$spearman2,
            predicted = final_pred,
            observed = y_test,
            selected_proteins = character(0)  # none selected
        )

        message(sprintf("Cohort %s — R² (clinical only): %.4f", test_cohort, r2_score(y_test, final_pred)))
    }

    return(results)
}


result_rf_LOCO_clinical_z <- loco_forward_selection_clinical_only(
    clinical_data = combined_clinical,
    protein_data = proteome_significant,
    method = "rf"
)

result_knn_LOCO_clinical_z <- loco_forward_selection_clinical_only(
    clinical_data = combined_clinical,
    protein_data = proteome_significant,
    method = "knn"
)

result_svr_LOCO_clinical_z <- loco_forward_selection_clinical_only(
    clinical_data = combined_clinical,
    protein_data = proteome_significant,
    method = "svr"
)

result_lasso_LOCO_clinical_z <- loco_forward_selection_clinical_only(
    clinical_data = combined_clinical,
    protein_data = proteome_significant,
    method = "lasso"
)


save(result_rf_LOCO_clinical_z, result_knn_LOCO_clinical_z, result_svr_LOCO_clinical_z, result_lasso_LOCO_clinical_z,
     file = "data/ML_Mvalue_prediction_data/loco_forward_selection_results_CLINICAL.RData")
















#### ALL DATA Mval prediction ####

detach("package:MOFA2", unload = TRUE)
library(glmnet)
library(caret)
library(xgboost)
library(ranger)
library(e1071)  # for SVR
library(FNN)    # for KNN
library(dplyr)
library(tidyr)
library(ggplot2)

# Helper: Compute R²
r2_score <- function(y_true, y_pred) {
    1 - mean((y_true - y_pred)^2) / var(y_true)
}

# Helper: Compute Pearson² and Spearman²
correlation_scores <- function(y_true, y_pred) {
    list(
        pearson2 = cor(y_true, y_pred, method = "pearson")^2,
        spearman2 = cor(y_true, y_pred, method = "spearman")^2
    )
}

#### Train model on all data ####
cv_forward_selection <- function(clinical_data, protein_data, target = "M_value_z_by_cohort",
                                 method = c("rf", "lasso", "svr", "knn", "xgb"),
                                 max_features = 20, cv_folds = 10, seed = 123) {

    set.seed(seed)
    method <- match.arg(method)
    clinical_data$Gender <- ifelse(clinical_data$Gender == "Male", 1, 0)

    common_ids <- intersect(colnames(protein_data), rownames(clinical_data))
    protein_data <- protein_data[, common_ids]
    clinical_data <- clinical_data[common_ids, ]
    target_vector <- clinical_data[[target]]
    names(target_vector) <- rownames(clinical_data)

    # Standardize clinical features
    clinical_vars <- clinical_data[, c("Age", "Gender", "BMI", "HbA1c", "TG_HDL_ratio")]
    pre_proc <- caret::preProcess(clinical_vars, method = c("center", "scale"))
    X_all <- predict(pre_proc, clinical_vars)

    selected <- c()
    remaining <- rownames(protein_data)
    r2_progress <- c()

    for (step in 1:max_features) {
        best_cv_r2 <- -Inf
        best_prot <- NULL

        for (prot in setdiff(remaining, selected)) {
            temp_X <- X_all
            temp_X[[prot]] <- as.numeric(protein_data[prot, rownames(clinical_data)])

            folds <- caret::createFolds(target_vector, k = cv_folds, returnTrain = TRUE)
            fold_r2s <- c()

            for (fold in folds) {
                X_train <- temp_X[fold, , drop = FALSE]
                y_train <- target_vector[fold]
                X_val   <- temp_X[-fold, , drop = FALSE]
                y_val   <- target_vector[-fold]

                pred <- switch(method,
                               rf = {
                                   fit <- ranger::ranger(y = y_train, x = X_train)
                                   predict(fit, data = X_val)$predictions
                               },
                               lasso = {
                                   fit <- glmnet::cv.glmnet(as.matrix(X_train), y_train, alpha = 1)
                                   predict(fit, newx = as.matrix(X_val), s = "lambda.min")
                               },
                               svr = {
                                   fit <- svm(x = X_train, y = y_train, kernel = "radial")
                                   predict(fit, newdata = X_val)
                               },
                               knn = {
                                   FNN::knn.reg(train = X_train, test = X_val, y = y_train, k = 5)$pred
                               },
                               xgb = {
                                   dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
                                   dval <- xgb.DMatrix(data = as.matrix(X_val))
                                   model <- xgboost(data = dtrain, nrounds = 100, objective = "reg:squarederror", verbose = 0)
                                   predict(model, dval)
                               }
                )

                fold_r2s <- c(fold_r2s, r2_score(y_val, pred))
            }

            mean_r2 <- mean(fold_r2s, na.rm = TRUE)
            if (mean_r2 > best_cv_r2) {
                best_cv_r2 <- mean_r2
                best_prot <- prot
            }
        }

        selected <- c(selected, best_prot)
        r2_progress <- c(r2_progress, best_cv_r2)
        X_all[[best_prot]] <- as.numeric(protein_data[best_prot, rownames(clinical_data)])
    }

    # Select best number of features
    best_k <- which.max(r2_progress)
    final_prots <- selected[1:best_k]
    message(sprintf("Selected %d proteins with best CV R² = %.4f", best_k, r2_progress[best_k]))

    # Evaluate final model using CV again (on full data)
    folds <- caret::createFolds(target_vector, k = cv_folds, returnTrain = TRUE)
    all_pred <- rep(NA, length(target_vector))
    names(all_pred) <- rownames(clinical_data)

    for (fold in folds) {
        X_train <- X_all[fold, , drop = FALSE]
        X_val   <- X_all[-fold, , drop = FALSE]
        y_train <- target_vector[fold]

        for (prot in final_prots) {
            X_train[[prot]] <- as.numeric(protein_data[prot, rownames(X_train)])
            X_val[[prot]]   <- as.numeric(protein_data[prot, rownames(X_val)])
        }

        y_pred <- switch(method,
                         rf = {
                             fit <- ranger::ranger(y = y_train, x = X_train)
                             predict(fit, data = X_val)$predictions
                         },
                         lasso = {
                             fit <- glmnet::cv.glmnet(as.matrix(X_train), y_train, alpha = 1)
                             predict(fit, newx = as.matrix(X_val), s = "lambda.min")
                         },
                         svr = {
                             fit <- svm(x = X_train, y = y_train, kernel = "radial")
                             predict(fit, newdata = X_val)
                         },
                         knn = {
                             FNN::knn.reg(train = X_train, test = X_val, y = y_train, k = 5)$pred
                         },
                         xgb = {
                             dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
                             dval <- xgb.DMatrix(data = as.matrix(X_val))
                             model <- xgboost(data = dtrain, nrounds = 100, objective = "reg:squarederror", verbose = 0)
                             predict(model, dval)
                         }
        )

        all_pred[rownames(X_val)] <- y_pred
    }

    y_true <- target_vector[match(names(all_pred), names(target_vector))]
    scores <- correlation_scores(y_true, all_pred)

    #### EXTRA FOR LASSO EXPORT ####
    final_model <- NULL
    scaling_params <- NULL
    active_features <- NULL

    if (method == "lasso") {
        # Rebuild full matrix with final features (clinical + proteins)
        for (prot in final_prots) {
            X_all[[prot]] <- as.numeric(protein_data[prot, rownames(clinical_data)])
        }

        full_X <- X_all[, c("Age", "Gender", "BMI", "HbA1c", "TG_HDL_ratio", final_prots)]
        final_model <- glmnet::cv.glmnet(as.matrix(full_X), target_vector, alpha = 1)
        scaling_params <- list(center = pre_proc$mean, scale = pre_proc$std)
        active_features <- colnames(full_X)
    }

    #### FINAL RETURN ####
    return(list(
        r2 = r2_score(y_true, all_pred),
        pearson2 = scores$pearson2,
        spearman2 = scores$spearman2,
        predicted = all_pred,
        observed = y_true,
        selected_proteins = final_prots,
        final_model = final_model,
        scaling_params = scaling_params,
        active_features = active_features
    ))
}



result_knn_cv_z <- cv_forward_selection(
    clinical_data = combined_clinical,
    protein_data = proteome_significant,
    method = "knn",
    max_features = 20,
    cv_folds = 10
)

result_svr_cv_z <- cv_forward_selection(
    clinical_data = combined_clinical,
    protein_data = proteome_significant,
    method = "svr",
    max_features = 20,
    cv_folds = 10
)

result_rf_cv_z <- cv_forward_selection(
    clinical_data = combined_clinical,
    protein_data = proteome_significant,
    method = "rf",
    max_features = 20,
    cv_folds = 10
)

result_lasso_cv_z <- cv_forward_selection(
    clinical_data = combined_clinical,
    protein_data = proteome_significant,
    method = "lasso",
    max_features = 20,
    cv_folds = 10
)
names(result_lasso_cv_z$final_model)
coefs <- coef(result_lasso_cv_z$final_model, s = "lambda.min")
coefs_df <- as.data.frame(as.matrix(coefs))
coefs_df$Feature <- rownames(coefs_df)
colnames(coefs_df)[1] <- "Coefficient"

# Filter to non-zero coefficients only
nonzero_coefs <- subset(coefs_df, Coefficient != 0)

result_lasso_cv_z$active_features

result_LASSO_final_model <- result_lasso_cv_z$final_model
result_LASSO_scaling_params <- result_lasso_cv_z$scaling_params
result_LASSO_final_active_features <- result_lasso_cv_z$active_features
result_LASSO_final_active_features[6:23] <- paste0(word(result_LASSO_final_active_features[6:23],1,sep="_"),"_",word(result_LASSO_final_active_features[6:23],3,sep="_"))

save(
    result_LASSO_final_model,
    result_LASSO_scaling_params,
    result_LASSO_final_active_features,
    file = "M_value_prediction_model.RData"
)

save(result_knn_cv_z, result_svr_cv_z, result_rf_cv_z, result_lasso_cv_z,
     file = "data/ML_Mvalue_prediction_data/ALL_CV_forward_selection_results.RData")

## which are active features ##
# Get the non-zero coefficient set for deployment
coefs <- coef(result_lasso_cv_z$final_model, s = "lambda.min")
coefs_df <- as.data.frame(as.matrix(coefs))
coefs_df$Feature <- rownames(coefs_df)
colnames(coefs_df)[1] <- "Coefficient"

nonzero <- subset(coefs_df, Coefficient != 0)
# separate intercept and features
intercept <- nonzero$Coefficient[nonzero$Feature == "(Intercept)"]
nonzero_features <- nonzero$Feature[nonzero$Feature != "(Intercept)"]

# These are the *actual* features used by the UKBB LASSO model
nonzero_features
intercept





#### ALL data Proteins only Mean-scored M value prediction ####
cv_forward_selection_protein_only <- function(clinical_data, protein_data, target = "M_value_z_by_cohort",
                                              method = c("rf", "lasso", "svr", "knn"),
                                              max_features = 20, cv_folds = 10, seed = 123) {

    set.seed(seed)
    method <- match.arg(method)

    common_ids <- intersect(colnames(protein_data), rownames(clinical_data))
    protein_data <- protein_data[, common_ids]
    clinical_data <- clinical_data[common_ids, ]
    target_vector <- clinical_data[[target]]
    names(target_vector) <- rownames(clinical_data)

    X_all <- data.frame(row.names = rownames(clinical_data))  # Start with empty data frame

    selected <- c()
    remaining <- rownames(protein_data)
    r2_progress <- c()

    for (step in 1:max_features) {
        best_cv_r2 <- -Inf
        best_prot <- NULL

        for (prot in setdiff(remaining, selected)) {
            temp_X <- X_all
            temp_X[[prot]] <- as.numeric(protein_data[prot, rownames(clinical_data)])

            folds <- caret::createFolds(target_vector, k = cv_folds, returnTrain = TRUE)
            fold_r2s <- c()

            for (fold in folds) {
                X_train <- temp_X[fold, , drop = FALSE]
                y_train <- target_vector[fold]
                X_val   <- temp_X[-fold, , drop = FALSE]
                y_val   <- target_vector[-fold]

                pred <- switch(method,
                               rf = {
                                   fit <- ranger::ranger(y = y_train, x = X_train)
                                   predict(fit, data = X_val)$predictions
                               },
                               lasso = {
                                   x_train_mat <- as.matrix(X_train)
                                   x_val_mat <- as.matrix(X_val)

                                   # glmnet needs 2+ columns; handle single-feature edge case
                                   if (ncol(x_train_mat) == 1) {
                                       # Add dummy column with 0s to avoid glmnet error
                                       x_train_mat <- cbind(x_train_mat, dummy = 0)
                                       x_val_mat <- cbind(x_val_mat, dummy = 0)

                                       fit <- glmnet::glmnet(x_train_mat, y_train, alpha = 1)
                                       predict(fit, newx = x_val_mat, s = 0.01)
                                   } else {
                                       fit <- glmnet::cv.glmnet(x_train_mat, y_train, alpha = 1)
                                       predict(fit, newx = x_val_mat, s = "lambda.min")
                                   }
                               }
                               ,
                               svr = {
                                   fit <- e1071::svm(x = X_train, y = y_train, kernel = "radial")
                                   predict(fit, newdata = X_val)
                               },
                               knn = {
                                   FNN::knn.reg(train = X_train, test = X_val, y = y_train, k = 5)$pred
                               },
                               xgb = {
                                   dtrain <- xgboost::xgb.DMatrix(data = as.matrix(X_train), label = y_train)
                                   dval <- xgboost::xgb.DMatrix(data = as.matrix(X_val))
                                   model <- xgboost::xgboost(data = dtrain, nrounds = 100, objective = "reg:squarederror", verbose = 0)
                                   predict(model, dval)
                               }
                )

                fold_r2s <- c(fold_r2s, r2_score(y_val, pred))
            }

            mean_r2 <- mean(fold_r2s, na.rm = TRUE)
            if (mean_r2 > best_cv_r2) {
                best_cv_r2 <- mean_r2
                best_prot <- prot
            }
        }

        selected <- c(selected, best_prot)
        r2_progress <- c(r2_progress, best_cv_r2)
        X_all[[best_prot]] <- as.numeric(protein_data[best_prot, rownames(clinical_data)])
    }

    best_k <- which.max(r2_progress)
    final_prots <- selected[1:best_k]
    message(sprintf("Selected %d proteins with best CV R² = %.4f", best_k, r2_progress[best_k]))

    folds <- caret::createFolds(target_vector, k = cv_folds, returnTrain = TRUE)
    all_pred <- rep(NA, length(target_vector))
    names(all_pred) <- rownames(clinical_data)

    for (fold in folds) {
        X_train <- data.frame(row.names = rownames(clinical_data)[fold])
        X_val   <- data.frame(row.names = rownames(clinical_data)[-fold])

        for (prot in final_prots) {
            X_train[[prot]] <- as.numeric(protein_data[prot, rownames(X_train)])
            X_val[[prot]]   <- as.numeric(protein_data[prot, rownames(X_val)])
        }

        y_train <- target_vector[fold]
        y_val   <- target_vector[-fold]

        y_pred <- switch(method,
                         rf = {
                             fit <- ranger::ranger(y = y_train, x = X_train)
                             predict(fit, data = X_val)$predictions
                         },
                         lasso = {
                             fit <- glmnet::cv.glmnet(as.matrix(X_train), y_train, alpha = 1)
                             predict(fit, newx = as.matrix(X_val), s = "lambda.min")
                         },
                         svr = {
                             fit <- e1071::svm(x = X_train, y = y_train, kernel = "radial")
                             predict(fit, newdata = X_val)
                         },
                         knn = {
                             FNN::knn.reg(train = X_train, test = X_val, y = y_train, k = 5)$pred
                         },
                         xgb = {
                             dtrain <- xgboost::xgb.DMatrix(data = as.matrix(X_train), label = y_train)
                             dval <- xgboost::xgb.DMatrix(data = as.matrix(X_val))
                             model <- xgboost::xgboost(data = dtrain, nrounds = 100, objective = "reg:squarederror", verbose = 0)
                             predict(model, dval)
                         }
        )

        all_pred[rownames(X_val)] <- y_pred
    }

    y_true <- target_vector[match(names(all_pred), names(target_vector))]
    scores <- correlation_scores(y_true, all_pred)

    return(list(
        r2 = r2_score(y_true, all_pred),
        pearson2 = scores$pearson2,
        spearman2 = scores$spearman2,
        predicted = all_pred,
        observed = y_true,
        selected_proteins = final_prots
    ))
}

result_knn_cv_noclin_z <- cv_forward_selection_protein_only(
    clinical_data = combined_clinical,
    protein_data = proteome_significant,
    method = "knn",
    max_features = 20,
    cv_folds = 10
)
result_knn_cv_noclin_z$r2

result_svr_cv_noclin_z <- cv_forward_selection_protein_only(
    clinical_data = combined_clinical,
    protein_data = proteome_significant,
    method = "svr",
    max_features = 20,
    cv_folds = 10
)

result_rf_cv_noclin_z <- cv_forward_selection_protein_only(
    clinical_data = combined_clinical,
    protein_data = proteome_significant,
    method = "rf",
    max_features = 20,
    cv_folds = 10
)

result_lasso_cv_noclin_z <- cv_forward_selection_protein_only(
    clinical_data = combined_clinical,
    protein_data = proteome_significant,
    method = "lasso",
    max_features = 20,
    cv_folds = 10
)

save(result_knn_cv_noclin_z, result_svr_cv_noclin_z, result_rf_cv_noclin_z, result_lasso_cv_noclin_z,
     file = "data/ML_Mvalue_prediction_data/ALL_CV_forward_selection_results_PROTEINS.RData")



#### Clinical features only - Mean scaled M-value ####
cv_model_clinical_only <- function(clinical_data, target = "M_value_z_by_cohort",
                                   method = c("rf", "lasso", "svr", "knn"),
                                   cv_folds = 10, seed = 123) {

    set.seed(seed)
    method <- match.arg(method)

    # Convert gender to numeric if necessary
    if (is.character(clinical_data$Gender) || is.factor(clinical_data$Gender)) {
        clinical_data$Gender <- ifelse(clinical_data$Gender == "Male", 1, 0)
    }

    # Extract predictors and target
    clinical_vars <- clinical_data[, c("Age", "Gender", "BMI", "HbA1c", "TG_HDL_ratio")]
    target_vector <- clinical_data[[target]]

    # Standardize predictors (across all data)
    pre_proc <- caret::preProcess(clinical_vars, method = c("center", "scale"))
    X_all <- predict(pre_proc, clinical_vars)

    # Cross-validation setup
    folds <- caret::createFolds(target_vector, k = cv_folds, returnTrain = FALSE)
    all_pred <- rep(NA, length(target_vector))
    names(all_pred) <- rownames(clinical_data)

    for (fold_idx in seq_along(folds)) {
        test_idx <- folds[[fold_idx]]
        train_idx <- setdiff(seq_along(target_vector), test_idx)

        X_train <- X_all[train_idx, , drop = FALSE]
        X_val   <- X_all[test_idx, , drop = FALSE]
        y_train <- target_vector[train_idx]
        y_val   <- target_vector[test_idx]

        y_pred <- switch(method,
                         rf = {
                             fit <- ranger::ranger(y = y_train, x = X_train)
                             predict(fit, data = X_val)$predictions
                         },
                         lasso = {
                             fit <- glmnet::cv.glmnet(as.matrix(X_train), y_train, alpha = 1)
                             predict(fit, newx = as.matrix(X_val), s = "lambda.min")
                         },
                         svr = {
                             fit <- e1071::svm(x = X_train, y = y_train, kernel = "radial")
                             predict(fit, newdata = X_val)
                         },
                         knn = {
                             FNN::knn.reg(train = X_train, test = X_val, y = y_train, k = 5)$pred
                         },
                         xgb = {
                             dtrain <- xgboost::xgb.DMatrix(data = as.matrix(X_train), label = y_train)
                             dval <- xgboost::xgb.DMatrix(data = as.matrix(X_val))
                             model <- xgboost::xgboost(data = dtrain, nrounds = 100, objective = "reg:squarederror", verbose = 0)
                             predict(model, dval)
                         }
        )

        all_pred[test_idx] <- y_pred
    }

    # Final performance evaluation
    valid_idx <- !is.na(all_pred)
    y_true <- target_vector[valid_idx]
    y_pred <- all_pred[valid_idx]

    scores <- correlation_scores(y_true, y_pred)

    return(list(
        r2 = r2_score(y_true, y_pred),
        pearson2 = scores$pearson2,
        spearman2 = scores$spearman2,
        predicted = y_pred,
        observed = y_true,
        used_features = colnames(X_all)
    ))

}



result_knn_clinical_z <- cv_model_clinical_only(
    clinical_data = combined_clinical,
    target = "M_value_z_by_cohort",
    method = "knn"
)

result_rf_clinical_z <- cv_model_clinical_only(
    clinical_data = combined_clinical,
    target = "M_value_z_by_cohort",
    method = "rf"
)

result_lasso_clinical_z <- cv_model_clinical_only(
    clinical_data = combined_clinical,
    target = "M_value_z_by_cohort",
    method = "lasso"
)

result_svr_clinical_z <- cv_model_clinical_only(
    clinical_data = combined_clinical,
    target = "M_value_z_by_cohort",
    method = "svr"
)


save(result_knn_clinical_z, result_rf_clinical_z, result_lasso_clinical_z, result_svr_clinical_z,
     file = "data/ML_Mvalue_prediction_data/ALL_CV_forward_selection_results_CLINICAL.RData")




#### Summary plot of all ML performance ####
# leave-out-one-cohort analysis
library(ggplot2)

# Combine all R² values into a single data frame
r2_summary <- data.frame(
    Model = rep(c("RF", "SVR", "KNN", "LASSO"), 3),
    R2 = c(
        result_rf_z$IDA$r2,
        result_svr_z$IDA$r2,
        result_knn_z$IDA$r2,
        result_lasso_z$IDA$r2,

        result_rf_LOCO_protein_Z$IDA$r2,
        result_svr_LOCO_protein_Z$IDA$r2,
        result_knn_LOCO_protein_Z$IDA$r2,
        result_lasso_LOCO_protein_Z$IDA$r2,

        result_rf_LOCO_clinical_z$IDA$r2,
        result_svr_LOCO_clinical_z$IDA$r2,
        result_knn_LOCO_clinical_z$IDA$r2,
        result_lasso_LOCO_clinical_z$IDA$r2
    ),
    Input = rep(c("Clinical_Protein", "Protein", "Clinical"), each = 4)
)

unique(r2_summary$Input)
class(r2_summary$Input)

r2_summary$Input <- factor(r2_summary$Input, levels=c("Clinical","Protein", "Clinical_Protein"))

ggplot(r2_summary, aes(x = Model, y = R2, fill = Input)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.75), width = 0.6) +
    scale_fill_manual(
        values = c(
            "Clinical" = "#e8eefa",
            "Protein" = "#8da0cb",
            "Clinical_Protein" = "#4564a3"
        ),
        labels = c("Clinical", "Protein", "Clinical + Protein")
    ) +
    labs(
        title = expression("M-value prediction (Validation cohort)"),
        x = "Model",
        y = expression(R^2),
        fill = "Input Features"
    ) +
    coord_cartesian(ylim = c(0, 0.8)) +
    theme_minimal(base_size = 14) +
    theme(
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = c(0.87, 0.95),
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2)
    )


#all cross validation
# Combine all R² values into a single data frame
r2_summary_CV <- data.frame(
    Model = rep(c("RF", "SVR", "KNN", "LASSO"), 3),
    R2 = c(
        result_rf_cv_z$r2,
        result_svr_cv_z$r2,
        result_knn_cv_z$r2,
        result_lasso_cv_z$r2,

        result_rf_cv_noclin_z$r2,
        result_svr_cv_noclin_z$r2,
        result_knn_cv_noclin_z$r2,
        result_lasso_cv_noclin_z$r2,

        result_rf_clinical_z$r2,
        result_svr_clinical_z$r2,
        result_knn_clinical_z$r2,
        result_lasso_clinical_z$r2
    ),
    Input = rep(c("Clinical_Protein", "Protein", "Clinical"), each = 4)
)

r2_summary_CV$Input <- factor(r2_summary_CV$Input, levels=c("Clinical","Protein", "Clinical_Protein"))

ggplot(r2_summary_CV, aes(x = Model, y = R2, fill = Input)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.75), width = 0.6) +
    scale_fill_manual(
        values = c(
            "Clinical" = "#e8eefa",
            "Protein" = "#8da0cb",
            "Clinical_Protein" = "#4564a3")) +
    labs(
        title = "M-value prediction (All data)",
        x = "Model",
        y = expression(R^2),
        fill = "Input Features"
    ) + coord_cartesian(ylim = c(0, 0.8)) +
    theme_minimal(base_size = 14) +
    theme(
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = c(0.95, 1),
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.2)
    )


library(ggplot2)

# Data
df_ida <- data.frame(
    Observed = result_lasso_z$IDA$observed,
    Predicted = as.numeric(result_lasso_z$IDA$predicted)
)

# Compute R²
r2_ida <- round(r2_score(df_ida$Observed, df_ida$Predicted), 2)

# Plot
ggplot(df_ida, aes(x = Observed, y = Predicted)) +
    geom_point(alpha = 0.6, size = 3.5, color = "darkblue", fill="#e8eefa", shape=21) +
    geom_smooth(method = "lm", se = TRUE, color = "#4564a3", fill = "#4564a3", alpha = 0.1, linetype="dashed") +
    annotate("text", x = min(df_ida$Observed), y = max(df_ida$Predicted),
             label = paste0("R² = ", r2_ida), hjust = 0, vjust = 1.2, size = 5) +
    theme_minimal(base_size = 14) +
    theme(
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(
        title = "LASSO model prediction (validation cohort)",
        x = "Observed M-value (Z-scored)",
        y = "Predicted M-value (Z-scored)"
    )


# Data
df_all <- data.frame(
    Observed = result_lasso_cv_z$observed,
    Predicted = as.numeric(result_lasso_cv_z$predicted)
)

# Compute R²
r2_all <- round(r2_score(df_all$Observed, df_all$Predicted), 2)

# Plot
ggplot(df_all, aes(x = Observed, y = Predicted)) +
    geom_point(alpha = 0.6, size = 3.5, color = "darkblue", fill="#e8eefa", shape=21) +
    geom_smooth(method = "lm", se = TRUE, color = "#4564a3", fill = "#4564a3", alpha = 0.1, linetype="dashed") +
    annotate("text", x = min(df_all$Observed), y = max(df_all$Predicted),
             label = paste0("R² = ", r2_all), hjust = 0, vjust = 1.2, size = 5) +
    theme_minimal(base_size = 14) +
    theme(
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(
        title = "LASSO model prediction (All data)",
        x = "Observed M-value (Z-scored)",
        y = "Predicted M-value (Z-scored)"
    )
