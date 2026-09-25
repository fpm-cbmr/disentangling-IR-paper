processed_datasets <- readRDS("data/processed_datasets.rds")
Clinical_data_filt <- readRDS("data/Clinical_data_filt")
Sample_ID_filt <- readRDS("data/Sample_ID_filt.rds")
rownames(processed_datasets$Olink_discovery_cohort) <- processed_datasets$Olink_discovery_cohort$ID

processed_datasets <- readRDS("data/processed_datasets.rds")
Clinical_data_filt <- readRDS("data/Clinical_data_filt")
Sample_ID_filt <- readRDS("data/Sample_ID_filt.rds")
rownames(processed_datasets$Olink_discovery_cohort) <- processed_datasets$Olink_discovery_cohort$ID

pre_clinical_data <- Clinical_data_filt[which(word(colnames(processed_datasets$MagNet_discovery_cohort$corrected),3,sep="-") == 1),]
rownames(pre_clinical_data) <- pre_clinical_data$Subject_ID
alternate.cols <- function(m1, m2) {
    cbind(m1, m2)[, order(c(seq(ncol(m1)), seq(ncol(m2))))]
}

library(ppcor)

perform_ir_partial_correlation <- function(
    proteomics_data,
    pre_clinical_data,
    covariate = "Gender",
    mvalue_col = "M..µmol..kg.min..",
    homa_col = "HOMA1.IR",
    insulin_col = "FS.Insulin..mIE.L.",
    ffa_col = "FFA..mmol.L.",
    pre_pattern = "-1$"
) {

  if (is.null(dim(proteomics_data))) {
    stop("proteomics_data is not a 2D object.")
  }

  if (is.null(colnames(proteomics_data))) {
    stop("proteomics_data has no column names.")
  }

  pre_idx <- grepl(pre_pattern, colnames(proteomics_data))

  if (!any(pre_idx)) {
    stop(paste("No Pre samples found with pattern:", pre_pattern))
  }

  pre_plasma_intensities <- proteomics_data[, pre_idx, drop = FALSE]
  colnames(pre_plasma_intensities) <- sub(pre_pattern, "", colnames(pre_plasma_intensities))

  matching_samples <- intersect(colnames(pre_plasma_intensities), rownames(pre_clinical_data))
  if (length(matching_samples) == 0) {
    stop("No matching samples between proteomics data and pre_clinical_data.")
  }

  pre_plasma_intensities <- pre_plasma_intensities[, matching_samples, drop = FALSE]
  pre_clinical_data <- pre_clinical_data[matching_samples, , drop = FALSE]

  required_cols <- c(mvalue_col, homa_col, insulin_col, ffa_col, covariate)
  missing_cols <- setdiff(required_cols, colnames(pre_clinical_data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing clinical columns:", paste(missing_cols, collapse = ", ")))
  }

  pre_clinical_data$AdipoIR <- as.numeric(pre_clinical_data[[insulin_col]]) *
    as.numeric(pre_clinical_data[[ffa_col]])
  pre_clinical_data$M_value <- as.numeric(pre_clinical_data[[mvalue_col]])
  pre_clinical_data$HOMA_IR <- as.numeric(pre_clinical_data[[homa_col]])

  clin_var <- c("M_value", "HOMA_IR", "AdipoIR")

  pre_plasma_intensities <- as.matrix(pre_plasma_intensities)
  mode(pre_plasma_intensities) <- "numeric"

  z_raw <- pre_clinical_data[[covariate]]
  if (is.character(z_raw) || is.factor(z_raw)) {
    z_raw <- as.character(z_raw)
    if (all(na.omit(z_raw) %in% c("Male", "Female"))) {
      z <- ifelse(z_raw == "Male", 1, 0)
    } else {
      stop(paste("Covariate", covariate, "is not numeric and not coded as Male/Female."))
    }
  } else {
    z <- as.numeric(z_raw)
  }

  baseline_p.val_matrix <- matrix(NA, nrow = nrow(pre_plasma_intensities), ncol = length(clin_var))
  baseline_estimate_matrix <- matrix(NA, nrow = nrow(pre_plasma_intensities), ncol = length(clin_var))

  rownames(baseline_p.val_matrix) <- rownames(pre_plasma_intensities)
  rownames(baseline_estimate_matrix) <- rownames(pre_plasma_intensities)

  for (k in seq_along(clin_var)) {
    y <- as.numeric(pre_clinical_data[[clin_var[k]]])

    for (i in seq_len(nrow(pre_plasma_intensities))) {
      x <- as.numeric(pre_plasma_intensities[i, ])
      valid_idx <- complete.cases(x, y, z)

      if (sum(valid_idx) < 5) next

      test_result <- tryCatch({
        ppcor::pcor.test(
          x = x[valid_idx],
          y = y[valid_idx],
          z = z[valid_idx],
          method = "spearman"
        )
      }, error = function(e) NULL)

      if (!is.null(test_result)) {
        baseline_p.val_matrix[i, k] <- test_result$p.value
        baseline_estimate_matrix[i, k] <- test_result$estimate
      }
    }
  }

  adj_p_matrix <- apply(baseline_p.val_matrix, 2, function(p) p.adjust(p, method = "BH"))
  if (is.vector(adj_p_matrix)) {
    adj_p_matrix <- matrix(adj_p_matrix, ncol = 1)
  }

  rownames(adj_p_matrix) <- rownames(pre_plasma_intensities)

  colnames(baseline_p.val_matrix) <- paste0(clin_var, "_p.value")
  colnames(adj_p_matrix) <- paste0(clin_var, "_adj.p.value")
  colnames(baseline_estimate_matrix) <- paste0(clin_var, "_estimate")

  complete_cor_matrix <- cbind(baseline_p.val_matrix, adj_p_matrix, baseline_estimate_matrix)

  organized_columns <- unlist(lapply(clin_var, function(var) {
    c(
      paste0(var, "_p.value"),
      paste0(var, "_adj.p.value"),
      paste0(var, "_estimate")
    )
  }))

  complete_cor_matrix <- complete_cor_matrix[, organized_columns, drop = FALSE]

  result_frame <- data.frame(
    ID = rownames(pre_plasma_intensities),
    complete_cor_matrix,
    row.names = NULL,
    check.names = FALSE
  )

  return(result_frame)
}


PCA_ir_partial_cor_frame <- perform_ir_partial_correlation(
  proteomics_data   = processed_datasets$PCA_discovery_cohort$corrected,
  pre_clinical_data = pre_clinical_data,
  covariate         = "Gender"
)

Neat_ir_partial_cor_frame <- perform_ir_partial_correlation(
  proteomics_data   = processed_datasets$Neat_discovery_cohort$corrected,
  pre_clinical_data = pre_clinical_data,
  covariate         = "Gender"
)

Depleted_ir_partial_cor_frame <- perform_ir_partial_correlation(
  proteomics_data   = processed_datasets$Depleted_discovery_cohort$corrected,
  pre_clinical_data = pre_clinical_data,
  covariate         = "Gender"
)

MagNet_ir_partial_cor_frame <- perform_ir_partial_correlation(
  proteomics_data   = processed_datasets$MagNet_discovery_cohort$corrected,
  pre_clinical_data = pre_clinical_data,
  covariate         = "Gender"
)

Olink_ir_partial_cor_frame <- perform_ir_partial_correlation(
  proteomics_data   = processed_datasets$Olink_discovery_cohort,
  pre_clinical_data = pre_clinical_data,
  covariate         = "Gender"
)


library(dplyr)
library(stringr)

# Reformat Olink IDs: GENE_number_UNIPROT -> UNIPROT_GENE
Olink_ir_partial_cor_frame <- Olink_ir_partial_cor_frame %>%
  mutate(
    ID = str_replace(ID, "^([^_]+)_[^_]+_([^_]+)$", "\\2_\\1")
  )




extract_gene_ms <- function(df) {
  df %>%
    mutate(Gene = sub(".*_", "", ID))
}

extract_gene_olink <- function(df) {
  df %>%
    mutate(Gene = sub("_.*", "", ID))
}

PCA_ir_partial_cor_frame      <- extract_gene_ms(PCA_ir_partial_cor_frame)
Neat_ir_partial_cor_frame     <- extract_gene_ms(Neat_ir_partial_cor_frame)
Depleted_ir_partial_cor_frame <- extract_gene_ms(Depleted_ir_partial_cor_frame)
MagNet_ir_partial_cor_frame   <- extract_gene_ms(MagNet_ir_partial_cor_frame)

Olink_ir_partial_cor_frame    <- extract_gene_olink(Olink_ir_partial_cor_frame)

library(dplyr)

all_ir_cor <- bind_rows(
  PCA_ir_partial_cor_frame      %>% mutate(Workflow = "PCA"),
  Neat_ir_partial_cor_frame     %>% mutate(Workflow = "Neat"),
  Depleted_ir_partial_cor_frame %>% mutate(Workflow = "Depleted"),
  MagNet_ir_partial_cor_frame   %>% mutate(Workflow = "MagNet"),
  Olink_ir_partial_cor_frame    %>% mutate(Workflow = "Olink")
)

pick_best_per_gene <- function(df, p_col) {
  df %>%
    filter(!is.na(.data[[p_col]])) %>%
    group_by(Gene) %>%
    slice_min(order_by = .data[[p_col]], n = 1, with_ties = FALSE) %>%
    ungroup()
}


combined_M_value <- pick_best_per_gene(all_ir_cor, "M_value_p.value") %>%
  arrange(M_value_p.value) %>%
  mutate(Significant = ifelse(M_value_adj.p.value < 0.05, "YES", "NO"))

combined_HOMA_IR <- pick_best_per_gene(all_ir_cor, "HOMA_IR_p.value") %>%
  arrange(HOMA_IR_p.value) %>%
  mutate(Significant = ifelse(HOMA_IR_adj.p.value < 0.05, "YES", "NO"))

combined_AdipoIR <- pick_best_per_gene(all_ir_cor, "AdipoIR_p.value") %>%
  arrange(AdipoIR_p.value) %>%
  mutate(Significant = ifelse(AdipoIR_adj.p.value < 0.05, "YES", "NO"))


extract_uniprot <- function(id_vec) {
  ifelse(
    grepl("^[^_]+_[^_]+_[^_]+$", id_vec),  # Olink pattern (3 parts)
    sub("^[^_]+_[^_]+_([^_]+)$", "\\1", id_vec),  # take last
    sub("^([^_]+)_.*$", "\\1", id_vec)            # take first (MS)
  )
}


annotate_results <- function(df, estimate_col, adj_p_col) {
  df %>%
    mutate(
      UniProt_ID = extract_uniprot(ID),
      GeneID = Gene,
      Significant = ifelse(!is.na(.data[[adj_p_col]]) & .data[[adj_p_col]] < 0.05, "YES", "NO"),
      Direction = case_when(
        Significant == "YES" & .data[[estimate_col]] > 0  ~ "Up",
        Significant == "YES" & .data[[estimate_col]] < 0  ~ "Down",
        TRUE ~ NA_character_
      )
    )
}

combined_M_value <- pick_best_per_gene(all_ir_cor, "M_value_p.value") %>%
  arrange(M_value_p.value) %>%
  annotate_results("M_value_estimate", "M_value_adj.p.value")

combined_HOMA_IR <- pick_best_per_gene(all_ir_cor, "HOMA_IR_p.value") %>%
  arrange(HOMA_IR_p.value) %>%
  annotate_results("HOMA_IR_estimate", "HOMA_IR_adj.p.value")

combined_AdipoIR <- pick_best_per_gene(all_ir_cor, "AdipoIR_p.value") %>%
  arrange(AdipoIR_p.value) %>%
  annotate_results("AdipoIR_estimate", "AdipoIR_adj.p.value")


combined_M_value_subset <- combined_M_value[,13:16]
combined_HOMA_IR_subset <- combined_HOMA_IR[,13:16]
combined_AdipoIR_subset <- combined_AdipoIR[,13:16]


saveRDS(combined_M_value_subset, file = "data/IRindices_M_value_subset.rds")
saveRDS(combined_HOMA_IR_subset, file = "data/IRindices_HOMA_IR_subset.rds")
saveRDS(combined_AdipoIR_subset, file = "data/IRindices_AdipoIR_subset.rds")
