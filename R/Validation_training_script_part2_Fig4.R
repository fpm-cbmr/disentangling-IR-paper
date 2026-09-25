#### HIIT- discovery cohort validation M-value correlations ####
library(dplyr)

# Load clinical and expression data
Val_Clinical_data <- readRDS("data/Val_Clinical_data.rds")
datasets <- readRDS("data/processed_datasets_validation.rds")

# Subset only workflows you want (MS-based and/or Olink)
workflow_names <- names(processed_datasets_validation)

# Clinical variables to correlate
clinical_vars <- c("Age", "BMI1", "Fat_per1", "VO2max_kg1", "HbA1c1",
                   "F_gluc1", "F_ins1", "F_Cpep1", "GIR1", "RD_ba1",
                   "RD_cl1", "M.value")

# Spearman correlation function
run_correlation_analysis <- function(expr_mat_subset, clinical_data_subset, label = "Pre") {
  cor_results <- list()

  for (var in clinical_vars) {
    clinical_vector <- clinical_data_subset[[var]]

    cor_out <- apply(expr_mat_subset, 1, function(protein_expr) {
      test <- suppressWarnings(cor.test(protein_expr, clinical_vector, method = "spearman", use = "pairwise.complete.obs"))
      c(estimate = test$estimate, pval = test$p.value)
    })

    cor_df <- as.data.frame(t(cor_out))
    cor_df$adj.pval <- p.adjust(cor_df$pval, method = "BH")
    colnames(cor_df) <- paste0(label, "_", var, c("_rho", "_pval", "_adj.pval"))

    cor_results[[var]] <- cor_df
  }

  cor_matrix <- do.call(cbind, cor_results)
  cor_matrix$Protein <- rownames(expr_mat_subset)
  cor_matrix <- cor_matrix[, c("Protein", setdiff(names(cor_matrix), "Protein"))]
  return(cor_matrix)
}

# Subset sample indices
pre_idx <- Val_Clinical_data$Condition == "Pre"
post_idx <- Val_Clinical_data$Condition == "Post"

# Storage for all workflows
workflow_correlation_results <- list()

for (wf in workflow_names) {
  message("Running correlation for workflow: ", wf)

  # Extract expression matrix
  if (wf == "Olink_validation_cohort") {
    expr_mat <- processed_datasets_validation[[wf]]
  } else {
    expr_mat <- processed_datasets_validation[[wf]]$scaled
  }

  # Match to clinical data
  expr_pre <- expr_mat[, pre_idx]
  expr_post <- expr_mat[, post_idx]

  clin_pre <- Val_Clinical_data[pre_idx, ]
  clin_post <- Val_Clinical_data[post_idx, ]

  # Run analysis
  cor_pre <- run_correlation_analysis(expr_pre, clin_pre, label = "Pre")
  cor_post <- run_correlation_analysis(expr_post, clin_post, label = "Post")

  # Merge and store
  cor_combined <- left_join(cor_pre, cor_post, by = "Protein")
  workflow_correlation_results[[wf]] <- cor_combined
}



#### HIIT-OLINK validation M-value correlations ####
Val_Clinical_data_Olink <- readRDS("data/Val_Clinical_data_Olink.rds")
clinical_vars <- c("Age", "BMI1", "Fat_per1", "VO2max_kg1", "HbA1c1",
                   "F_gluc1", "F_ins1", "F_Cpep1", "GIR1", "RD_ba1",
                   "RD_cl1", "M.value")


# Function to perform correlation analysis
run_correlation_analysis <- function(expr_mat_subset, clinical_data_subset, label = "Pre") {
  cor_results <- list()

  for (var in clinical_vars) {
    clinical_vector <- clinical_data_subset[[var]]

    cor_out <- apply(expr_mat_subset, 1, function(protein_expr) {
      test <- cor.test(protein_expr, clinical_vector, method = "spearman", use = "pairwise.complete.obs")
      c(estimate = test$estimate, pval = test$p.value)
    })

    cor_df <- as.data.frame(t(cor_out))
    cor_df$adj.pval <- p.adjust(cor_df$pval, method = "BH")
    colnames(cor_df) <- paste0(label, "_", var, c("_rho", "_pval", "_adj.pval"))

    cor_results[[var]] <- cor_df
  }

  # Combine all into one table
  cor_matrix <- do.call(cbind, cor_results)
  cor_matrix$Protein <- rownames(expr_mat_subset)
  cor_matrix <- cor_matrix[, c("Protein", setdiff(names(cor_matrix), "Protein"))]
  return(cor_matrix)
}

# Subset for Pre and Post samples
pre_idx <- Val_Clinical_data_Olink$Condition == "Pre"
post_idx <- Val_Clinical_data_Olink$Condition == "Post"

expr_pre <- processed_datasets_validation$Olink_validation_cohort[, pre_idx]
expr_post <- processed_datasets_validation$Olink_validation_cohort[, post_idx]

clin_pre <- Val_Clinical_data_Olink[pre_idx, ]
clin_post <- Val_Clinical_data_Olink[post_idx, ]

# Run the analyses
cor_pre <- run_correlation_analysis(expr_pre, clin_pre, label = "Pre")
cor_post <- run_correlation_analysis(expr_post, clin_post, label = "Post")

# Combine Olink pre and post correlation results by Protein
cor_combined_olink <- left_join(cor_pre, cor_post, by = "Protein")

cor_combined_olink$Protein # change ID

cor_combined_olink$Protein <- paste0(word(cor_combined_olink$Protein,3,sep="_"),"_",word(cor_combined_olink$Protein,1,sep="_"))

# Store it in the same list as used for MS workflows
workflow_correlation_results$Olink_validation_cohort <- cor_combined_olink

workflow_correlation_results$Neat_validation_cohort$Protein <- paste0(str_extract(workflow_correlation_results$Neat_validation_cohort$Protein, "^[^;_]+"),"_",
    str_extract(workflow_correlation_results$Neat_validation_cohort$Protein, "(?<=_)[^;]+"))
workflow_correlation_results$PCA_validation_cohort$Protein <- paste0(str_extract(workflow_correlation_results$PCA_validation_cohort$Protein, "^[^;_]+"),"_",
                                                                     str_extract(workflow_correlation_results$PCA_validation_cohort$Protein, "(?<=_)[^;]+"))
workflow_correlation_results$MagNet_validation_cohort$Protein <- paste0(str_extract(workflow_correlation_results$MagNet_validation_cohort$Protein, "^[^;_]+"),"_",
                                                                        str_extract(workflow_correlation_results$MagNet_validation_cohort$Protein, "(?<=_)[^;]+"))
workflow_correlation_results$Depleted_validation_cohort$Protein <- paste0(str_extract(workflow_correlation_results$Depleted_validation_cohort$Protein, "^[^;_]+"),"_",
                                                                          str_extract(workflow_correlation_results$Depleted_validation_cohort$Protein, "(?<=_)[^;]+"))


# Step 1: Combine all workflows into one data.frame with workflow label
combined_df <- bind_rows(
  lapply(names(workflow_correlation_results), function(wf) {
    df <- workflow_correlation_results[[wf]]
    df$Workflow <- wf
    df
  })
)

combined_df$Workflow <- word(combined_df$Workflow,1,sep="_")

# Step 2: Filter down to most significant M-value Pre-pval per unique Protein
# (You can also use adj.pval instead if preferred)
final_df <- combined_df %>%
  filter(!is.na(M.value.Pre_M.value_pval)) %>%
  group_by(Protein) %>%
  slice_min(order_by = M.value.Pre_M.value_pval, with_ties = FALSE) %>%
  ungroup()

Discovery_cohort_clin_cor <- read.delim(file="All_unique_significant_results.txt")
Discovery_cohort_clin_cor$Protein <- paste0(Discovery_cohort_clin_cor$UniProt_ID,"_",Discovery_cohort_clin_cor$ID)


Cohort_HIIT_combined <- left_join(
  Discovery_cohort_clin_cor,
  combined_df,
  by = c("Protein", "Workflow")
)


library(ggplot2)


ggplot(Cohort_HIIT_combined, aes(x = M.value.Pre_M.value_rho, y = M.value_estimate_pre, color=Workflow)) +
  geom_point(size = 2, alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Rho Insulin sensitivity association (Validation cohort)",
       y = "Rho Insulin sensitivity association (Discovery cohort)") +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),              # remove gridlines
    panel.background = element_blank(),        # ensure no background
    axis.line = element_line(color = "black"), # add black axis lines
    axis.text = element_text(color = "black"), # black axis text
    axis.title = element_text(color = "black"),# black axis titles
    plot.title = element_text(color = "black", face = "bold", hjust = 0.5)
  ) + scale_color_manual(values = c("darkred","darkgreen", "gray","orange", "darkblue"))


cor.test(Cohort_HIIT_combined$M.value.Pre_M.value_rho,Cohort_HIIT_combined$M.value_estimate_pre, method="spearman")

#### Training-induced changes in insulin sensitivity-assocaited proteins ####
dataset_names <- c(
  "Overall_training_neat",
  "Overall_training_PCA",
  "Overall_training_MagNet",
  "Overall_training_Depleted",
  "Overall_training_Olink",

  "Groupwise_training_neat",
  "Groupwise_training_PCA",
  "Groupwise_training_Magnet",
  "Groupwise_training_Depleted",
  "Groupwise_training_Olink",

  "ANOVA_training_neat",
  "ANOVA_training_PCA",
  "ANOVA_training_Magnet",
  "ANOVA_training_Depleted",
  "ANOVA_training_Olink",

  "Group_diff_neat",
  "Group_diff_PCA",
  "Group_diff_MagNet",
  "Group_diff_Depleted",
  "Group_diff_Olink",

  "Group_ANOVA_neat",
  "Group_ANOVA_PCA",
  "Group_ANOVA_MagNet",
  "Group_ANOVA_Depleted",
  "Group_ANOVA_Olink"
)

datasets <- list()

for (name in dataset_names) {
  file_rds <- file.path("data/validation_LMM_data/", paste0(name, ".rds"))
  datasets[[name]] <- readRDS(file_rds)
}

datasets$Overall_training_Olink$Protein <- sapply(datasets$Overall_training_Olink$Protein, function(x) {
  parts <- unlist(strsplit(x, "_"))
  parts <- parts[!grepl("^\\d+$", parts)]  # Remove purely numeric segments
  paste0(tail(parts, 1), "_", parts[1])
})

datasets$Overall_training_neat$Protein <- paste0(str_extract(datasets$Overall_training_neat$Protein, "^[^;_]+"),"_",
                                                                      str_extract(datasets$Overall_training_neat$Protein, "(?<=_)[^;]+"))
datasets$Overall_training_PCA$Protein <- paste0(str_extract(datasets$Overall_training_PCA$Protein, "^[^;_]+"),"_",
                                                                     str_extract(datasets$Overall_training_PCA$Protein, "(?<=_)[^;]+"))
datasets$Overall_training_MagNet$Protein <- paste0(str_extract(datasets$Overall_training_MagNet$Protein, "^[^;_]+"),"_",
                                                                        str_extract(datasets$Overall_training_MagNet$Protein, "(?<=_)[^;]+"))
datasets$Overall_training_Depleted$Protein <- paste0(str_extract(datasets$Overall_training_Depleted$Protein, "^[^;_]+"),"_",
                                                                          str_extract(datasets$Overall_training_Depleted$Protein, "(?<=_)[^;]+"))

# Add BH correction and a Workflow column to each dataset
datasets$Overall_training_Olink <- datasets$Overall_training_Olink %>%
  mutate(adj.P_value = p.adjust(P_value, method = "BH"),
         Workflow = "Olink")

datasets$Overall_training_neat <- datasets$Overall_training_neat %>%
  mutate(adj.P_value = p.adjust(P_value, method = "BH"),
         Workflow = "Neat")

datasets$Overall_training_PCA <- datasets$Overall_training_PCA %>%
  mutate(adj.P_value = p.adjust(P_value, method = "BH"),
         Workflow = "PCA")

datasets$Overall_training_MagNet <- datasets$Overall_training_MagNet %>%
  mutate(adj.P_value = p.adjust(P_value, method = "BH"),
         Workflow = "MagNet")

datasets$Overall_training_Depleted <- datasets$Overall_training_Depleted %>%
  mutate(adj.P_value = p.adjust(P_value, method = "BH"),
         Workflow = "Depleted")

# Combine all into one long data frame
all_combined <- bind_rows(
  datasets$Overall_training_Olink,
  datasets$Overall_training_neat,
  datasets$Overall_training_PCA,
  datasets$Overall_training_MagNet,
  datasets$Overall_training_Depleted
)

# Keep only the row with the lowest raw P_value per protein
final_combined <- all_combined %>%
  group_by(Protein) %>%
  slice_min(order_by = P_value, with_ties = FALSE) %>%
  ungroup()


library(ggplot2)
library(ggrepel)
final_combined$label <- word(final_combined$Protein,2,sep="_")
#export in 6x6 PDF
ggplot(final_combined[final_combined$adj.P_value < 0.05,], aes(x = Estimate, y = -log10(P_value), color = Workflow)) +
  geom_point(alpha = 0.4, size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    x = "Effect Size (Standardized Estimate)",
    y = "-log10(P-value)",
    title = "Training-induced protein changes",
    color = "Workflow"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    plot.title = element_text(color = "black", face = "bold", hjust = 0.5)
  ) +
  scale_color_manual(values = c("darkred","grey","orange")) + ylim(3,15) + geom_text_repel(aes(label = label))





#### HIIT M-value associated proteins responding to training ####
HIIT_Mval_and_training <- inner_join(final_combined,final_df, by=c("Protein","Workflow"))
HIIT_Mval_and_training <- inner_join(all_combined,combined_df, by=c("Protein","Workflow"))

Discovery_cohort_clin_cor <- read.delim(file="All_unique_significant_results.txt")
Discovery_cohort_clin_cor$Protein <- paste0(Discovery_cohort_clin_cor$UniProt_ID,"_",Discovery_cohort_clin_cor$ID)

Mvalue_discovery_Training_response <- inner_join(Discovery_cohort_clin_cor,all_combined, by=c("Protein","Workflow"))



ggplot(Mvalue_discovery_Training_response[Mvalue_discovery_Training_response$M.value_adj.p.value_pre < 0.05,],
       aes(M.value_estimate_pre, Estimate, color = abs(Estimate))) +
  geom_point(alpha = 0.4, size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    x = "Rho insulin sensitivity association (Discovery cohort)",
    y = "Effect Size (Standardized Estimate)",
    title = ""
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(color = "black", face = "bold", hjust = 0.5)
  ) +
  scale_color_gradient(low = "grey", high = "darkred") +
  guides(color = "none")


library(dplyr)
library(stringr)

protein_df <- final_combined %>%
    mutate(
        UniProt_ID = str_extract(Protein, "^[^_]+"),        # before underscore
        GeneID = str_extract(Protein, "(?<=_).*"),         # after underscore
        Is_Significant = ifelse(P_value < 0.01, "Yes", "No"),
        Direction = case_when(
            Is_Significant == "Yes" & Estimate > 0 ~ "Up",
            Is_Significant == "Yes" & Estimate < 0 ~ "Down",
            TRUE ~ NA_character_
        )
    ) %>%
    dplyr::select(GeneID, UniProt_ID, Is_Significant, Direction)

