library(lme4)
library(emmeans)
library(car)
library(dplyr)
library(broom.mixed)
library(tidyr)
library(purrr)
library(data.table)
library(readr)
library(pheatmap)
library(patchwork)
library(ggrepel)
library(lmerTest)

#### LMM code to get overall training effect and group-specific responses ####
#Interception was to be "Post" sample by default why I take "-" of the estimate afterwards.
processed_datasets_validation <- readRDS("data/processed_datasets_validation.rds")
Val_Clinical_data <- readRDS("data/Val_Clinical_data.rds")

# Ensure Training and Group are factors
Val_Clinical_data$Training <- factor(Val_Clinical_data$Condition, levels = c("Pre", "Post"))
Val_Clinical_data$Group <- factor(Val_Clinical_data$Group, levels = c("1", "2", "3"))  # 1 = Lean, 2 = Obese, 3 = Obese+T2D

# Standardize (Z-score) protein levels for all workflows
#for (workflow in names(processed_datasets_validation)) {
#  processed_datasets_validation[[workflow]]$scaled <- t(apply(
#    processed_datasets_validation[[workflow]]$scaled, 1, scale
#  ))
#}

for (workflow in setdiff(names(processed_datasets_validation), "Olink_validation_cohort")) {
    processed_datasets_validation[[workflow]]$scaled <- t(apply(
        processed_datasets_validation[[workflow]]$scaled, 1, scale
    ))
}



# Keep M-value in original units
Val_Clinical_data$M.value <- Val_Clinical_data$M.value

# Initialize lists to store results
overall_training_results <- list()
groupwise_training_results <- list()
anova_results <- list()

# Loop through all workflows
for (workflow in names(processed_datasets_validation)) {

  cat("\nRunning analysis for:", workflow, "\n")

  # Get standardized protein dataset for this workflow
  protein_data <- processed_datasets_validation[[workflow]]$scaled

  for (protein in rownames(protein_data)) {

    # Prepare dataset for this protein
    valid_data <- data.frame(
      ID = Val_Clinical_data$ID,
      Training = Val_Clinical_data$Training,
      Group = Val_Clinical_data$Group,
      Protein_Level = as.numeric(protein_data[protein, ])
    ) %>% na.omit()  # Remove NAs

    # Ensure enough paired samples
    if (nrow(valid_data) < 5) {
      overall_training_results[[workflow]][[protein]] <- data.frame(Protein = protein, Estimate = 0, P_value = 1)
      groupwise_training_results[[workflow]][[protein]] <- data.frame(Protein = protein, Group = NA, Estimate = 0, P_value = 1)
      anova_results[[workflow]][[protein]] <- data.frame(Protein = protein, Anova_P_value = 1)
      next
    }

    # Try fitting the model
    model <- tryCatch(
      lmer(Protein_Level ~ Training * Group + (1 | ID), data = valid_data),
      error = function(e) return(NULL)
    )

    if (is.null(model)) {
      overall_training_results[[workflow]][[protein]] <- data.frame(Protein = protein, Estimate = 0, P_value = 1)
      groupwise_training_results[[workflow]][[protein]] <- data.frame(Protein = protein, Group = NA, Estimate = 0, P_value = 1)
      anova_results[[workflow]][[protein]] <- data.frame(Protein = protein, Anova_P_value = 1)
      next
    }

    # 1️⃣ Extract Overall Training Effect (using contrasts)
    overall_emm <- tryCatch(
      emmeans(model, ~ Training),
      error = function(e) return(NULL)
    )

    if (!is.null(overall_emm)) {
      contrast_overall <- tryCatch(
        as.data.frame(contrast(overall_emm, "pairwise")),  # Convert to dataframe
        error = function(e) return(NULL)
      )

      if (!is.null(contrast_overall)) {
        overall_training_results[[workflow]][[protein]] <- data.frame(
          Protein = protein,
          Estimate = contrast_overall$estimate[1],  # Training effect (Post vs. Pre)
          P_value = contrast_overall$p.value[1]  # p-value for Training effect
        )
      } else {
        overall_training_results[[workflow]][[protein]] <- data.frame(Protein = protein, Estimate = 0, P_value = 1)
      }
    } else {
      overall_training_results[[workflow]][[protein]] <- data.frame(Protein = protein, Estimate = 0, P_value = 1)
    }

    # 2️⃣ Extract Group-Specific Training Effects
    groupwise_emm <- tryCatch(
      emmeans(model, pairwise ~ Training | Group),
      error = function(e) return(NULL)
    )

    if (!is.null(groupwise_emm)) {
      groupwise_contrasts <- tryCatch(
        as.data.frame(groupwise_emm$contrasts),  # Convert to dataframe
        error = function(e) return(NULL)
      )

      if (!is.null(groupwise_contrasts)) {
        groupwise_contrasts$Protein <- protein
        groupwise_training_results[[workflow]][[protein]] <- groupwise_contrasts
      } else {
        groupwise_training_results[[workflow]][[protein]] <- data.frame(Protein = protein, Group = NA, Estimate = 0, P_value = 1)
      }
    } else {
      groupwise_training_results[[workflow]][[protein]] <- data.frame(Protein = protein, Group = NA, Estimate = 0, P_value = 1)
    }

    # 3️⃣ ANOVA to Test if Training Effects Differ Across Groups
    anova_test <- tryCatch(
      Anova(model, type = "III"),
      error = function(e) return(NULL)
    )

    if (!is.null(anova_test)) {
      anova_results[[workflow]][[protein]] <- data.frame(
        Protein = protein,
        Anova_P_value = anova_test["Training:Group", "Pr(>Chisq)"]
      )
    } else {
      anova_results[[workflow]][[protein]] <- data.frame(Protein = protein, Anova_P_value = 1)
    }
    # 4️⃣ Group contrasts from the same model
    group_contrasts <- tryCatch(
        emmeans(model, pairwise ~ Group),
        error = function(e) return(NULL)
    )

    if (!is.null(group_contrasts)) {
        group_contrasts_df <- tryCatch(
            as.data.frame(group_contrasts$contrasts),
            error = function(e) return(NULL)
        )
        if (!is.null(group_contrasts_df)) {
            group_contrasts_df$Protein <- protein
            if (!exists("group_diff_results")) group_diff_results <- list()
            if (is.null(group_diff_results[[workflow]])) {
                group_diff_results[[workflow]] <- list()
            }
            group_diff_results[[workflow]][[protein]] <- group_contrasts_df
        }
    }

    # 5️⃣ ANOVA for main group effect
    group_anova_p <- tryCatch({
        anova_test["Group", "Pr(>Chisq)"]
    }, error = function(e) 1)

    if (!exists("anova_group_results")) anova_group_results <- list()
    if (is.null(anova_group_results[[workflow]])) {
        anova_group_results[[workflow]] <- list()
    }

    anova_group_results[[workflow]][[protein]] <- data.frame(
        Protein = protein,
        Group_Anova_P_value = group_anova_p
    )
  }
}

# Convert lists to dataframes
overall_training_df <- lapply(overall_training_results, bind_rows)
groupwise_training_df <- lapply(groupwise_training_results, bind_rows)
anova_training_df <- lapply(anova_results, bind_rows)
group_diff_df <- lapply(group_diff_results, bind_rows)
anova_group_df <- lapply(anova_group_results, bind_rows)

# correcting sign of estimate
Overall_training_neat <- overall_training_df$Neat
Overall_training_neat$Estimate <- -Overall_training_neat$Estimate
Overall_training_PCA <- overall_training_df$PCA
Overall_training_PCA$Estimate <- -Overall_training_PCA$Estimate
Overall_training_MagNet <- overall_training_df$MagNet
Overall_training_MagNet$Estimate <- -Overall_training_MagNet$Estimate
Overall_training_Depleted <- overall_training_df$Depleted
Overall_training_Depleted$Estimate <- -Overall_training_Depleted$Estimate

groupwise_training_neat <- groupwise_training_df$Neat
groupwise_training_neat$estimate <- -groupwise_training_df$Neat$estimate

groupwise_training_PCA <- groupwise_training_df$PCA
groupwise_training_PCA$estimate <- -groupwise_training_df$PCA$estimate

groupwise_training_Magnet <- groupwise_training_df$MagNet
groupwise_training_Magnet$estimate <- -groupwise_training_df$MagNet$estimate

groupwise_training_Depleted <- groupwise_training_df$Depleted
groupwise_training_Depleted$estimate <- -groupwise_training_df$Depleted$estimate

anova_training_neat <- anova_training_df$Neat
anova_training_PCA <- anova_training_df$PCA
anova_training_Magnet <- anova_training_df$MagNet
anova_training_Depleted <- anova_training_df$Depleted


# Function to process each workflow
process_workflow <- function(df) {
    df %>%
        mutate(
            contrast = case_when(
                contrast == "Group1 - Group2" ~ "Obese_vs_Lean",
                contrast == "Group1 - Group3" ~ "T2D_vs_Lean",
                contrast == "Group2 - Group3" ~ "T2D_vs_Obese",
                TRUE ~ contrast
            ),
            estimate = -estimate  # Flip sign so it's X - Lean
        ) %>%
        pivot_longer(cols = c(estimate, SE, df, t.ratio, p.value), names_to = "stat", values_to = "value") %>%
        unite("contrast_stat", contrast, stat) %>%
        pivot_wider(names_from = contrast_stat, values_from = value)
}

# Apply to all workflows
group_diff_wide <- list(
    Neat = group_diff_df$Neat_validation_cohort,
    PCA = group_diff_df$PCA_validation_cohort,
    MagNet = group_diff_df$MagNet_validation_cohort,
    Depleted = group_diff_df$Depleted_validation_cohort
) %>%
    map(~ group_by(.x, Protein) %>% process_workflow())


# running Olink separately
# Ensure Training and Group are factors
Val_Clinical_data_Olink <- readRDS("data/Val_Clinical_data_Olink.rds")
Val_Clinical_data_Olink$Training <- factor(Val_Clinical_data_Olink$Condition, levels = c("Pre", "Post"))
Val_Clinical_data_Olink$Group <- factor(Val_Clinical_data_Olink$Group, levels = c("1", "2", "3"))

# Standardize protein levels
expr_mat_scaled <- t(apply(processed_datasets_validation$Olink_validation_cohort, 1, scale))
colnames(expr_mat_scaled) <- colnames(processed_datasets_validation$Olink_validation_cohort)

# Initialize lists to store results
overall_training_results <- list()
groupwise_training_results <- list()
anova_results <- list()
group_diff_results_olink <- list()
anova_group_results_olink <- list()

# Loop through proteins
for (protein in rownames(expr_mat_scaled)) {
    valid_data <- data.frame(
        ID = Val_Clinical_data_Olink$ID,
        Training = Val_Clinical_data_Olink$Training,
        Group = Val_Clinical_data_Olink$Group,
        Protein_Level = as.numeric(expr_mat_scaled[protein, ])
    ) %>% na.omit()

    if (nrow(valid_data) < 5) {
        overall_training_results[[protein]] <- data.frame(Protein = protein, Estimate = 0, P_value = 1)
        groupwise_training_results[[protein]] <- data.frame(Protein = protein, Group = NA, Estimate = 0, P_value = 1)
        anova_results[[protein]] <- data.frame(Protein = protein, Anova_P_value = 1)
        next
    }

    model <- tryCatch(
        lmer(Protein_Level ~ Training * Group + (1 | ID), data = valid_data),
        error = function(e) return(NULL)
    )

    if (is.null(model)) {
        overall_training_results[[protein]] <- data.frame(Protein = protein, Estimate = 0, P_value = 1)
        groupwise_training_results[[protein]] <- data.frame(Protein = protein, Group = NA, Estimate = 0, P_value = 1)
        anova_results[[protein]] <- data.frame(Protein = protein, Anova_P_value = 1)
        next
    }

    overall_emm <- tryCatch(emmeans(model, ~ Training), error = function(e) return(NULL))
    if (!is.null(overall_emm)) {
        contrast_overall <- tryCatch(as.data.frame(contrast(overall_emm, "pairwise")), error = function(e) return(NULL))
        if (!is.null(contrast_overall)) {
            overall_training_results[[protein]] <- data.frame(
                Protein = protein,
                Estimate = -contrast_overall$estimate[1],
                P_value = contrast_overall$p.value[1]
            )
        }
    }

    groupwise_emm <- tryCatch(emmeans(model, pairwise ~ Training | Group), error = function(e) return(NULL))
    if (!is.null(groupwise_emm)) {
        groupwise_contrasts <- tryCatch(as.data.frame(groupwise_emm$contrasts), error = function(e) return(NULL))
        if (!is.null(groupwise_contrasts)) {
            groupwise_contrasts$Protein <- protein
            groupwise_contrasts$estimate <- -groupwise_contrasts$estimate
            groupwise_training_results[[protein]] <- groupwise_contrasts
        }
    }

    anova_test <- tryCatch(Anova(model, type = "III"), error = function(e) return(NULL))
    if (!is.null(anova_test)) {
        anova_results[[protein]] <- data.frame(
            Protein = protein,
            Anova_P_value = anova_test["Training:Group", "Pr(>Chisq)"]
        )
    }
    # 4️⃣ Group contrasts
    group_contrasts <- tryCatch(
        emmeans(model, pairwise ~ Group),
        error = function(e) return(NULL)
    )

    if (!is.null(group_contrasts)) {
        group_contrasts_df <- tryCatch(
            as.data.frame(group_contrasts$contrasts),
            error = function(e) return(NULL)
        )
        if (!is.null(group_contrasts_df)) {
            group_contrasts_df$Protein <- protein
            group_diff_results_olink[[protein]] <- group_contrasts_df
        }
    }

    # 5️⃣ ANOVA main effect for Group
    group_anova_p <- tryCatch({
        anova_test["Group", "Pr(>Chisq)"]
    }, error = function(e) 1)

    anova_group_results_olink[[protein]] <- data.frame(
        Protein = protein,
        Group_Anova_P_value = group_anova_p
    )
}

# Convert results to dataframes
Overall_training_Olink <- bind_rows(overall_training_results)
groupwise_training_Olink <- bind_rows(groupwise_training_results)
anova_training_Olink <- bind_rows(anova_results)
Group_diff_Olink <- bind_rows(group_diff_results_olink)
Group_anova_Olink <- bind_rows(anova_group_results_olink)

# First, bind all results into a single data frame
group_diff_Olink_raw <- bind_rows(group_diff_results_olink)

# Then, process and reshape using the same function
group_diff_Olink <- group_diff_Olink_raw %>%
    group_by(Protein) %>%
    process_workflow()


#### Saving the LMM results ####
# Define the output directory
output_dir <- "data/validation_LMM_data"

# Ensure the directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# List of datasets with corresponding filenames
datasets <- list(
  "Overall_training_neat" = Overall_training_neat,
  "Overall_training_PCA" = Overall_training_PCA,
  "Overall_training_MagNet" = Overall_training_MagNet,
  "Overall_training_Depleted" = Overall_training_Depleted,
  "Overall_training_Olink" = Overall_training_Olink,

  "Groupwise_training_neat" = groupwise_training_neat,
  "Groupwise_training_PCA" = groupwise_training_PCA,
  "Groupwise_training_Magnet" = groupwise_training_Magnet,
  "Groupwise_training_Depleted" = groupwise_training_Depleted,
  "Groupwise_training_Olink" = groupwise_training_Olink,

  "ANOVA_training_neat" = anova_training_neat,
  "ANOVA_training_PCA" = anova_training_PCA,
  "ANOVA_training_Magnet" = anova_training_Magnet,
  "ANOVA_training_Depleted" = anova_training_Depleted,
  "ANOVA_training_Olink" = anova_training_Olink,

  "Group_diff_neat" = group_diff_wide$Neat,
  "Group_diff_PCA" = group_diff_wide$PCA,
  "Group_diff_MagNet" = group_diff_wide$MagNet,
  "Group_diff_Depleted" = group_diff_wide$Depleted,
  "Group_diff_Olink" = group_diff_Olink,

  "Group_ANOVA_neat" = anova_group_df$Neat,
  "Group_ANOVA_PCA" = anova_group_df$PCA,
  "Group_ANOVA_MagNet" = anova_group_df$MagNet,
  "Group_ANOVA_Depleted" = anova_group_df$Depleted,
  "Group_ANOVA_Olink" = Group_anova_Olink
)

# Save each dataset as a .txt and .rds file
for (name in names(datasets)) {
  file_txt <- file.path(output_dir, paste0(name, ".txt"))
  file_rds <- file.path(output_dir, paste0(name, ".rds"))

  # Save as tab-separated text file
  write_tsv(datasets[[name]], file_txt)

  # Save as RDS file
  saveRDS(datasets[[name]], file_rds)
}




#### overlapping M-value associated proteins with group-wise specific training responses ####
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

# check.names = FALSE keeps the hyphenated column names (e.g. `M-value_estimate_pre`)
# that this script refers to with backticks further down.
All_unique_significant <- read.delim("data-raw/All_unique_significant_results.txt",
                                    check.names = FALSE)


# List of groupwise training dataframes
groupwise_training_list <- list(
  Neat = datasets$Groupwise_training_neat,
  PCA = datasets$Groupwise_training_PCA,
  MagNet = datasets$Groupwise_training_Magnet,
  Depleted = datasets$Groupwise_training_Depleted,
  Olink = datasets$Groupwise_training_Olink
)

groupwise_training_list$Olink$ID <- groupwise_training_list$Olink$Protein
groupwise_training_list$Olink$Protein <- sapply(groupwise_training_list$Olink$Protein, function(x) {
    parts <- unlist(strsplit(x, "_"))
    parts <- parts[!grepl("^\\d+$", parts)]  # Remove purely numeric segments
    paste0(tail(parts, 1), "_", parts[1])
})

#groupwise_training_list$Olink$Group <- factor(
#    groupwise_training_list$Olink$Group,
#    levels = c("lean", "obese", "t2d"),
#    labels = c("1", "2", "3")
#)

groupwise_training_list$Neat$Protein <- paste0(str_extract(groupwise_training_list$Neat$Protein, "^[^;_]+"),"_",
                                                 str_extract(groupwise_training_list$Neat$Protein, "(?<=_)[^;]+"))
groupwise_training_list$PCA$Protein <- paste0(str_extract(groupwise_training_list$PCA$Protein, "^[^;_]+"),"_",
                                               str_extract(groupwise_training_list$PCA$Protein, "(?<=_)[^;]+"))
groupwise_training_list$MagNet$Protein <- paste0(str_extract(groupwise_training_list$MagNet$Protein, "^[^;_]+"),"_",
                                               str_extract(groupwise_training_list$MagNet$Protein, "(?<=_)[^;]+"))
groupwise_training_list$Depleted$Protein <- paste0(str_extract(groupwise_training_list$Depleted$Protein, "^[^;_]+"),"_",
                                               str_extract(groupwise_training_list$Depleted$Protein, "(?<=_)[^;]+"))


# Extract unique significant proteins from "All_unique_significant"
unique_significant_proteins <- unique(All_unique_significant$ID)  # Gene ID
unique_significant_uniprot <- unique(All_unique_significant$UniProt_ID)  # UniProt Accession

# Initialize list to store results
matching_proteins_list <- list()

# Loop through each workflow dataset and find matches
for (workflow in names(groupwise_training_list)) {

  # Extract the groupwise training dataset
  df <- groupwise_training_list[[workflow]]

  # Extract the protein names from "Protein" column
  df$Gene_ID <- sub(".*_", "", df$Protein)  # Extract gene name
  df$UniProt_ID <- sub("_.*", "", df$Protein)  # Extract Uniprot accession

  # Find proteins that are in All_unique_significant
  matched_proteins <- df %>%
    filter(Gene_ID %in% unique_significant_proteins | UniProt_ID %in% unique_significant_uniprot) %>%
    dplyr::select(Protein, Group, estimate, p.value, Gene_ID, UniProt_ID) %>%  # Include p.value
    mutate(Workflow = workflow)  # Add workflow column

  # Store the results
  matching_proteins_list[[workflow]] <- matched_proteins
}

# Combine all workflow results
combined_matching_proteins <- bind_rows(matching_proteins_list)

# Filter for significant results (P < 0.05)
# Filter for proteins that have at least one significant result in any Group
significant_proteins <- combined_matching_proteins %>%
  group_by(Protein) %>%
  filter(any(p.value < 0.01)) %>%  # Keep if at least one p-value < 0.01
  ungroup()

# Merge with All_unique_significant to include M-value estimates
final_results <- All_unique_significant %>%
  dplyr::select(ID, UniProt_ID, `M-value_estimate_pre`, Workflow) %>%
  left_join(significant_proteins, by = c("ID" = "Gene_ID", "UniProt_ID", "Workflow"))


#alternative approach taking the workflow with most significant P-value #
best_workflow_per_protein <- significant_proteins %>%
    group_by(Gene_ID, UniProt_ID) %>%
    slice_min(order_by = p.value, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(Gene_ID, UniProt_ID, Workflow)

significant_proteins_filtered <- significant_proteins %>%
    inner_join(best_workflow_per_protein,
               by = c("Gene_ID", "UniProt_ID", "Workflow"))

final_results <- significant_proteins_filtered %>%
    left_join(All_unique_significant %>%
                  select(ID, UniProt_ID, `M-value_estimate_pre`),
              by = c("Gene_ID" = "ID", "UniProt_ID"))


# Save results
write.csv(final_results, "Matched_Proteins_Significant.csv", row.names = FALSE)


# Remove NAs and ensure Group is a factor
final_results_filtered <- final_results %>%
  filter(Protein %in% significant_proteins$Protein) %>%  # Keep only significant proteins
  mutate(Group = factor(Group, levels = c("1", "2", "3"), labels = c("Lean", "Obese", "Obese T2D")),  # Rename Groups
         M_value_direction = ifelse(`M-value_estimate_pre` >= 0, "Positive", "Negative"))  # Categorize sign


# Prepare data for clustering (wide format)
heatmap_data <- final_results_filtered %>%
  dplyr::select(Protein, Group, estimate, M_value_direction) %>%
  pivot_wider(names_from = Group, values_from = estimate) %>%
  column_to_rownames(var = "Protein")  # Set proteins as row names for clustering

# Convert to matrix for pheatmap
heatmap_matrix <- as.matrix(heatmap_data %>% dplyr::select(-M_value_direction))  # Remove annotation column
rownames(heatmap_matrix) <- word(rownames(heatmap_matrix),2,sep="_")

# Create row annotation for M-value direction
row_annotation <- data.frame(M_value_Direction = heatmap_data$M_value_direction)
rownames(row_annotation) <- rownames(heatmap_matrix)  # Match annotation to proteins

# Define color palettes
heatmap_colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)  # Blue = ↓, Red = ↑
annotation_colors <- list(
  M_value_Direction = c("Positive" = "#1B9E77", "Negative" = "#D95F02")  # Green for positive, Orange for negative
)

# Transpose the matrix
heatmap_matrix_flipped <- t(heatmap_matrix)

# Transpose annotations too (if needed)
annotation_col <- row_annotation  # annotation_row becomes column annotation

# Generate flipped heatmap
pheatmap(
    heatmap_matrix_flipped,
    color = heatmap_colors,
    cluster_rows = FALSE,          # Rows are now original columns (groups)
    cluster_cols = TRUE,           # Columns are now original proteins
    scale = "none",
    border_color = NA,
    main = "Horizontally Flipped Heatmap",
    fontsize = 12,
    show_rownames = TRUE,
    show_colnames = TRUE,
    treeheight_col = 50,           # Now dendrogram is on columns
    legend = TRUE,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors
)

#### Specific examples of group-specific changes ####
PCA_standardized <- t(apply(processed_datasets_validation$PCA_validation_cohort$scaled, 1, scale))
colnames(PCA_standardized) <- colnames(processed_datasets_validation$PCA_validation_cohort$scaled)
Neat_standardized <- t(apply(processed_datasets_validation$Neat_validation_cohort$scaled, 1, scale))
colnames(Neat_standardized) <- colnames(processed_datasets_validation$Neat_validation_cohort$scaled)
Olink_standardized <- t(apply(processed_datasets_validation$Olink_validation_cohort, 1, scale))
colnames(Olink_standardized) <- colnames(processed_datasets_validation$Olink_validation_cohort)

HRC_expr <- as.data.frame(PCA_standardized["P23327_HRC", , drop = FALSE])
HRC_expr <- as.numeric(HRC_expr)
names(HRC_expr) <- colnames(PCA_standardized)

plot_HRC <- data.frame(
    Expression = HRC_expr,
    Sample = names(HRC_expr)
) %>%
    left_join(Val_Clinical_data, by = c("Sample" = "ID_Cond"))


plot_HRC$Group <- factor(plot_HRC$Group)
plot_HRC$Condition <- factor(plot_HRC$Condition, levels=c("Pre","Post"))


HRC_final_plot <- ggplot(plot_HRC, aes(x = Condition, y = Expression)) +
    # Paired subject lines
    geom_line(aes(group = ID), color = "gray70", alpha = 0.5, linewidth = 0.6) +

    # Boxplots per condition
    geom_boxplot(aes(fill = Condition), alpha = 0.4, outlier.shape = NA, width = 0.5) +

    # Points colored by condition
    geom_point(aes(color = Condition), size = 2) +

    facet_wrap(~ Group) +
    theme_minimal(base_size = 14) +
    theme(
        axis.text = element_text(color = "black"),
        strip.text = element_text(size = 12, face = "bold"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")
    ) +
    labs(
        title = "Protein",
        y = "Standardized Expression (Z-score)",
        x = NULL
    ) +
    scale_fill_manual(values = c("grey", "darkred")) +
    scale_color_manual(values = c("grey", "darkred"))



PON3_expr <- Neat_standardized["Q15166_PON3", , drop = FALSE]
PON3_expr <- as.numeric(PON3_expr)
names(PON3_expr) <- colnames(PCA_standardized)

plot_PON3 <- data.frame(
    Expression = PON3_expr,
    Sample = names(PON3_expr)
) %>%
    left_join(Val_Clinical_data, by = c("Sample" = "ID_Cond"))


plot_PON3$Group <- factor(plot_PON3$Group)
plot_PON3$Condition <- factor(plot_PON3$Condition, levels=c("Pre","Post"))


PON3_final_plot <- ggplot(plot_PON3, aes(x = Condition, y = Expression)) +
    # Paired subject lines
    geom_line(aes(group = ID), color = "gray70", alpha = 0.5, linewidth = 0.6) +

    # Boxplots per condition
    geom_boxplot(aes(fill = Condition), alpha = 0.4, outlier.shape = NA, width = 0.5) +

    # Points colored by condition
    geom_point(aes(color = Condition), size = 2) +

    facet_wrap(~ Group) +
    theme_minimal(base_size = 14) +
    theme(
        axis.text = element_text(color = "black"),
        strip.text = element_text(size = 12, face = "bold"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black")
    ) +
    labs(
        title = "Protein",
        y = "Standardized Expression (Z-score)",
        x = NULL
    ) +
    scale_fill_manual(values = c("grey", "darkred")) +
    scale_color_manual(values = c("grey", "darkred"))



(HRC_final_plot / PON3_final_plot) +
    plot_layout(guides = "collect") & theme(legend.position = "right")




#### main effect of training ####
# List of groupwise training dataframes
main_training_list <- list(
    Neat = datasets$Overall_training_neat,
    PCA = datasets$Overall_training_PCA,
    MagNet = datasets$Overall_training_MagNet,
    Depleted = datasets$Overall_training_Depleted,
    Olink = datasets$Overall_training_Olink
)

main_training_list$Olink$ID <- main_training_list$Olink$Protein
main_training_list$Olink$Protein <- paste0(word(main_training_list$Olink$ID,3,sep="_"),"_",word(main_training_list$Olink$ID,1,sep="_"))
main_training_list$Olink$Group <- factor(
    main_training_list$Olink$Group,
    levels = c("lean", "obese", "t2d"),
    labels = c("1", "2", "3")
)


# Extract unique significant proteins from "All_unique_significant"
unique_significant_proteins <- unique(All_unique_significant$ID)  # Gene ID
unique_significant_uniprot <- unique(All_unique_significant$UniProt_ID)  # UniProt Accession

# Initialize list to store results
matching_proteins_list <- list()

# Loop through each workflow dataset and find matches
for (workflow in names(main_training_list)) {

    # Extract the groupwise training dataset
    df <- main_training_list[[workflow]]

    # Extract the protein names from "Protein" column
    df$Gene_ID <- sub(".*_", "", df$Protein)  # Extract gene name
    df$UniProt_ID <- sub("_.*", "", df$Protein)  # Extract Uniprot accession

    # Find proteins that are in All_unique_significant
    matched_proteins <- df %>%
        filter(Gene_ID %in% unique_significant_proteins | UniProt_ID %in% unique_significant_uniprot) %>%
        dplyr::select(Protein, Estimate, P_value, Gene_ID, UniProt_ID) %>%  # Include p.value
        mutate(Workflow = workflow)  # Add workflow column

    # Store the results
    matching_proteins_list[[workflow]] <- matched_proteins
}

# Combine all workflow results
combined_matching_proteins_training <- bind_rows(matching_proteins_list)

# Filter for significant results (P < 0.01)
# Filter for proteins that have at least one significant result in any Group
significant_proteins_training <- combined_matching_proteins_training %>%
    filter(P_value < 0.01) %>%  # Filter to only significant hits
    group_by(Protein) %>%
    slice_min(P_value, with_ties = FALSE) %>%  # Keep row with smallest p-value per protein
    ungroup()

final_results_training <- All_unique_significant %>%
    dplyr::select(ID, UniProt_ID, `M-value_estimate_pre`, Workflow) %>%
    left_join(significant_proteins_training, by = c("ID" = "Gene_ID", "UniProt_ID", "Workflow"))

# Remove NAs and ensure Group is a factor
final_results_filtered_training <- final_results_training %>%
    filter(Protein %in% significant_proteins$Protein)




#### progession plot Lean, obese, T2D ####
# List of workflows
workflows <- c("neat", "PCA", "MagNet", "Depleted", "Olink")

datasets$Group_ANOVA_Olink$Protein <- paste0(word(datasets$Group_ANOVA_Olink$Protein,3,sep="_"),"_",word(datasets$Group_ANOVA_Olink$Protein,1,sep="_"))
datasets$Group_diff_Olink$Protein <- paste0(word(datasets$Group_diff_Olink$Protein,3,sep="_"),"_",word(datasets$Group_diff_Olink$Protein,1,sep="_"))

# Step 1: BH-correct p-values and find most significant workflow per protein
sig_proteins <- list()

for (w in workflows) {
    anova_df <- datasets[[paste0("Group_ANOVA_", w)]]
    anova_df$adj.P.Val <- p.adjust(anova_df$Group_Anova_P_value, method = "BH")
    sig_hits <- anova_df %>% filter(adj.P.Val < 0.05)

    for (i in seq_len(nrow(sig_hits))) {
        prot <- sig_hits$Protein[i]
        pval <- sig_hits$adj.P.Val[i]

        if (!prot %in% names(sig_proteins) || pval < sig_proteins[[prot]]$adj.P.Val) {
            sig_proteins[[prot]] <- list(
                Workflow = w,
                adj.P.Val = pval
            )
        }
    }
}

# Step 2: Extract group-wise estimates from Group_diff_*
library(dplyr)
library(tidyverse)

protein_contrasts <- lapply(names(sig_proteins), function(prot) {
    w <- sig_proteins[[prot]]$Workflow
    df <- datasets[[paste0("Group_diff_", w)]]

    row <- df %>% filter(Protein == prot)
    if (nrow(row) == 0) return(NULL)

    data.frame(
        Protein = prot,
        Workflow = w,
        Obese_vs_Lean = row$Obese_vs_Lean_estimate,
        T2D_vs_Lean = row$T2D_vs_Lean_estimate,
        T2D_vs_Obese = row$T2D_vs_Obese_estimate
    )
})

# Combine all into a single data frame
protein_contrasts_df <- bind_rows(protein_contrasts)
long_profiles <- protein_contrasts_df %>%
    mutate(Lean = 0) %>%
    dplyr::select(Protein, Lean, Obese_vs_Lean, T2D_vs_Lean) %>%
    dplyr::rename(Obese = Obese_vs_Lean, T2D = T2D_vs_Lean) %>%
    pivot_longer(cols = c("Lean", "Obese", "T2D"), names_to = "Group", values_to = "Estimate") %>%
    mutate(Group = factor(Group, levels = c("Lean", "Obese", "T2D")))

# One row per protein, columns = Lean, Obese, T2D (Lean always 0)
wide_profiles <- protein_contrasts_df %>%
    mutate(Lean = 0) %>%
    dplyr::select(Protein, Lean, Obese_vs_Lean, T2D_vs_Lean) %>%
    dplyr::rename(Obese = Obese_vs_Lean, T2D = T2D_vs_Lean)


# Matrix of just values (no protein names)
mat <- as.matrix(wide_profiles[, c("Lean", "Obese", "T2D")])
rownames(mat) <- wide_profiles$Protein

# Optional: scale by row (center each protein)
mat_scaled <- t(scale(t(mat)))

# K-means clustering (choose k based on inspection or elbow method)
set.seed(123)
wss <- sapply(1:10, function(k) {
    kmeans(mat_scaled, centers = k, nstart = 25)$tot.withinss
})

plot(1:10, wss, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Total within-cluster sum of squares",
     main = "Elbow Method")


k <- 4
km <- kmeans(mat_scaled, centers = k)

# Add cluster assignment
wide_profiles$Cluster <- factor(km$cluster)

long_profiles_clustered <- long_profiles %>%
    left_join(wide_profiles %>% dplyr::select(Protein, Cluster), by = "Protein")


highlight_proteins <- long_profiles_clustered %>%
    group_by(Protein, Cluster) %>%
    summarise(range = max(Estimate) - min(Estimate), .groups = "drop") %>%
    group_by(Cluster) %>%
    slice_max(order_by = range, n = 5) %>%  # adjust n as desired
    pull(Protein)


highlight_proteins <- word(word(highlight_proteins,2,sep="_"),1,sep=";")
long_profiles_clustered$label <- word(word(long_profiles_clustered$Protein,2,sep="_"),1,sep=";")

ggplot(long_profiles_clustered, aes(x = Group, y = Estimate, group = Protein, color = Cluster)) +
    geom_line(alpha = 0.5,linewidth = 0.6) +
    facet_wrap(~ Cluster, ncol = 4, scales = "fixed") +
    geom_text_repel(
        data = long_profiles_clustered %>% filter(label %in% highlight_proteins, Group == "T2D"),
        aes(label =label),
        size = 3, segment.color = "gray60", max.overlaps = Inf
    ) +
    theme_minimal(base_size = 13) +
    labs(title = "Protein patterns across metabolic states",
         x = "Group", y = "Standardized effect size (relative to Lean group)") +
    theme(
        legend.position = "none",
        axis.text = element_text(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.spacing = unit(1.2, "lines"),
        strip.text = element_text(size = 12, face = "bold"),
        panel.grid = element_blank()  # remove all grid lines
    ) + geom_hline(yintercept = 0, linetype="dashed") +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1")


#### summary progression plot ####
library(dplyr)

summary_profiles <- long_profiles_clustered %>%
    group_by(Cluster, Group) %>%
    summarise(
        Median = median(Estimate, na.rm = TRUE),
        SD = sd(Estimate, na.rm = TRUE),
        SEM = SD / sqrt(n())
    ) %>%
    ungroup()


library(ggplot2)

ggplot(summary_profiles, aes(x = Group, y = Median, group = Cluster, color = Cluster)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = Median - SD, ymax = Median + SD, fill = Cluster), alpha = 0.2, color = NA) +
    labs(title = "Protein patterns across metabolic states", y = "Standardized effect size (relative to Lean group)", x = NULL) +
    theme_minimal(base_size = 15) +
    theme(
        panel.grid = element_blank(),
        axis.line = element_blank(),  # Remove default axis lines
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        plot.title = element_text(color = "black", face = "bold", hjust = 0.5),
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black")
    ) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    geom_segment(aes(x = 0.75, xend = 0.75, y = min(summary_profiles$Median - summary_profiles$SD),
                     yend = max(summary_profiles$Median + summary_profiles$SD)),
                 inherit.aes = FALSE, color = "black", linewidth = 0.5) + geom_hline(yintercept = 0, linewidth = 0.75, alpha=0.5, linetype="dashed")




#### Showcase 2 exaples ####
protein_id <- "Q15166_PON3"
protein_id <- "P26572_MGAT1"
protein_id <- "P58166_INHBE"
protein_id <- "Q03167_TGFBR3"

best_workflow <- sig_proteins[[protein_id]]$Workflow

expr_matrix <- processed_datasets_validation[[paste0(best_workflow, "_validation_cohort")]]$scaled

protein_expr <- expr_matrix[protein_id, , drop = FALSE]

plot_df <- data.frame(
    Expression = as.numeric(protein_expr),
    Sample = colnames(protein_expr)
) %>%
    left_join(Val_Clinical_data[, c("ID_Cond", "Group", "Condition")],
              by = c("Sample" = "ID_Cond"))

plot_df$Expression_z <- scale(plot_df$Expression)[, 1]

lean_mean_z <- plot_df %>%
    filter(Group == "1", Condition == "Pre") %>%
    summarise(mean_val = mean(Expression_z, na.rm = TRUE)) %>%
    pull(mean_val)

plot_df$Expression_z_lean_centered <- plot_df$Expression_z - lean_mean_z


plot_df$Group <- factor(plot_df$Group, levels = c("1", "2", "3"), labels = c("Lean", "Obese", "Obese T2D"))

MGAT1_plot <- ggplot(plot_df[plot_df$Condition == "Pre",], aes(x = Group, y = Expression_z_lean_centered)) +

    # Boxplots per group, purple color
    geom_boxplot(fill = "white", color = "purple", width = 0.5, outlier.shape = NA, alpha = 0.2) +

    # Jittered points, also purple
    geom_jitter(color = "purple", size = 3, width = 0.25, height = 0, alpha=0.4, stroke=NA) +

    # Line connecting group means
    stat_summary(fun = mean, geom = "line", aes(group = 1), color = "purple", linewidth = 0.4) +

    theme_minimal(base_size = 14) +
    labs(
        title = protein_id,
        y = "Z-scored Expression (Lean-centered)",
        x = NULL
    ) +
    theme(
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),  # Adds x and y axes
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        plot.title = element_text(color = "black", face = "bold", hjust = 0.5),
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black")
    )


protein_id <- "Q9P2J2_IGSF9"
best_workflow <- sig_proteins[[protein_id]]$Workflow
expr_matrix <- processed_datasets_validation[[paste0(best_workflow, "_validation_cohort")]]

protein_expr <- expr_matrix[grepl("IGSF9",rownames(expr_matrix)), , drop = FALSE]


plot_df <- data.frame(
    Expression = as.numeric(protein_expr),
    Sample = colnames(protein_expr)
) %>%
    left_join(Val_Clinical_data_Olink[, c("match_col", "Group", "Condition")],
              by = c("Sample" = "match_col"))


plot_df$Expression_z <- scale(plot_df$Expression)[, 1]

lean_mean_z <- plot_df %>%
    filter(Group == "1", Condition == "Pre") %>%
    summarise(mean_val = mean(Expression_z, na.rm = TRUE)) %>%
    pull(mean_val)

plot_df$Expression_z_lean_centered <- plot_df$Expression_z - lean_mean_z


plot_df$Group <- factor(plot_df$Group, levels = c("1", "2", "3"), labels = c("Lean", "Obese", "Obese T2D"))

IGSF9_plot <- ggplot(plot_df[plot_df$Condition == "Pre",], aes(x = Group, y = Expression_z_lean_centered)) +

    # Boxplots per group, purple color
    geom_boxplot(fill = "white", color = "darkgreen", width = 0.5, outlier.shape = NA, alpha = 0.2) +

    # Jittered points, also purple
    geom_jitter(color = "darkgreen", size = 3, width = 0.25, height = 0,alpha=0.4, stroke=NA) +

    # Line connecting group means
    stat_summary(fun = mean, geom = "line", aes(group = 1), color = "darkgreen", linewidth = 0.4) +

    theme_minimal(base_size = 14) +
    labs(
        title = protein_id,
        y = "Z-scored Expression (Lean-centered)",
        x = NULL
    ) +
    theme(
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),  # Adds x and y axes
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        plot.title = element_text(color = "black", face = "bold", hjust = 0.5),
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black")
    )



(IGSF9_plot / MGAT1_plot) +
    plot_layout(guides = "collect") & theme(legend.position = "right")
