#### Protein clinical association ####
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


#### Pre partial correlation adjusted for sex ####
library(ppcor)

perform_partial_correlation <- function(proteomics_data, pre_clinical_data, clin_var, covariate) {
    # Filter data for Pre samples
    pre_plasma_intensities <- proteomics_data[, grepl("-1$", colnames(proteomics_data))]
    if (ncol(pre_plasma_intensities) == 0) stop("No Pre samples found in proteomics data.")
    colnames(pre_plasma_intensities) <- sub("-1$", "", colnames(pre_plasma_intensities))

    # Align samples between proteomics data and clinical data
    matching_samples <- intersect(colnames(pre_plasma_intensities), rownames(pre_clinical_data))
    pre_plasma_intensities <- pre_plasma_intensities[, matching_samples, drop = FALSE]
    pre_clinical_data <- pre_clinical_data[matching_samples, , drop = FALSE]

    # Ensure numeric data
    pre_plasma_intensities <- as.data.frame(pre_plasma_intensities)
    pre_plasma_intensities[] <- lapply(pre_plasma_intensities, as.numeric)
    pre_clinical_data[, covariate] <- ifelse(pre_clinical_data[, covariate] == "Male", 1, 0)

    # Initialize matrices to store results
    baseline_p.val_matrix <- matrix(nrow = nrow(pre_plasma_intensities), ncol = length(clin_var))
    baseline_estimate_matrix <- matrix(nrow = nrow(pre_plasma_intensities), ncol = length(clin_var))

    for (k in 1:length(clin_var)) {
        if (!clin_var[k] %in% colnames(pre_clinical_data)) stop(paste("Clinical variable", clin_var[k], "not found in pre_clinical_data."))

        for (i in 1:nrow(pre_plasma_intensities)) {
            # Extract the protein data for all samples
            x <- as.numeric(pre_plasma_intensities[i, ])

            # Extract the clinical variable data for all samples
            y <- as.numeric(pre_clinical_data[, clin_var[k]])
            z <- as.numeric(pre_clinical_data[, covariate])

            # Filter valid samples
            valid_idx <- complete.cases(x, y, z)

            # Perform partial Spearman correlation on valid data
            test_result <- tryCatch({
                pcor.test(
                    x = x[valid_idx],
                    y = y[valid_idx],
                    z = z[valid_idx],
                    method = "spearman"
                )
            }, error = function(e) list(p.value = NA, estimate = NA))

            baseline_p.val_matrix[i, k] <- test_result$p.value
            baseline_estimate_matrix[i, k] <- test_result$estimate
        }
    }

    # Adjust p-values
    new_matrix_adjust <- matrix(nrow = nrow(baseline_p.val_matrix), ncol = ncol(baseline_p.val_matrix))
    for (i in 1:ncol(baseline_p.val_matrix)) {
        new_matrix_adjust[, i] <- p.adjust(baseline_p.val_matrix[, i], method = "BH")
    }

    # Set column names
    colnames(baseline_p.val_matrix) <- paste0(clin_var, "_p.value")
    colnames(baseline_estimate_matrix) <- paste0(clin_var, "_estimate")
    colnames(new_matrix_adjust) <- paste0(clin_var, "_adj.p.value")

    # Combine adjusted and unadjusted p-values and estimates
    complete_p.val_matrix <- cbind(baseline_p.val_matrix, new_matrix_adjust)
    complete_cor_matrix <- cbind(complete_p.val_matrix, baseline_estimate_matrix)

    # Reorganize columns: p.value, adj.p.value, estimate for each clinical variable
    organized_columns <- unlist(lapply(clin_var, function(var) {
        c(paste0(var, "_p.value"), paste0(var, "_adj.p.value"), paste0(var, "_estimate"))
    }))
    complete_cor_matrix <- complete_cor_matrix[, organized_columns]

    # Create data frame with results
    result_frame <- data.frame(
        ID = rownames(pre_plasma_intensities),
        complete_cor_matrix,
        row.names = NULL
    )

    # Return the result
    return(result_frame)
}


MagNet_partial_cor_frame <- perform_partial_correlation(
    processed_datasets$MagNet_discovery_cohort$corrected,
    pre_clinical_data,
    clin_var,
    "Gender"
)


Neat_partial_cor_frame <- perform_partial_correlation(
    processed_datasets$Neat_discovery_cohort$corrected,
    pre_clinical_data,
    clin_var,
    "Gender"
)


PCA_partial_cor_frame <- perform_partial_correlation(
    processed_datasets$PCA_discovery_cohort$corrected,
    pre_clinical_data,
    clin_var,
    "Gender"
)

Depleted_partial_cor_frame <- perform_partial_correlation(
    processed_datasets$Depleted_discovery_cohort$corrected,
    pre_clinical_data,
    clin_var,
    "Gender"
)

Olink_partial_cor_frame <- perform_partial_correlation(
    processed_datasets$Olink_discovery_cohort[,-1],
    pre_clinical_data,
    clin_var,
    "Gender"
)


#### Post partial correlation adjusted for sex ####
post_clinical_data <- Clinical_data_filt[which(word(colnames(processed_datasets$MagNet_discovery_cohort$corrected), 3, sep = "-") == "2"), ]
rownames(post_clinical_data) <- post_clinical_data$Subject_ID

perform_partial_correlation_post <- function(proteomics_data, post_clinical_data, clin_var, covariate) {
    # Filter data for Post samples
    post_plasma_intensities <- proteomics_data[, grepl("-2$", colnames(proteomics_data))]
    if (ncol(post_plasma_intensities) == 0) stop("No Post samples found in proteomics data.")
    colnames(post_plasma_intensities) <- sub("-2$", "", colnames(post_plasma_intensities))

    # Align samples between proteomics data and clinical data
    matching_samples <- intersect(colnames(post_plasma_intensities), rownames(post_clinical_data))
    post_plasma_intensities <- post_plasma_intensities[, matching_samples, drop = FALSE]
    post_clinical_data <- post_clinical_data[matching_samples, , drop = FALSE]

    # Ensure numeric data
    post_plasma_intensities <- as.data.frame(post_plasma_intensities)
    post_plasma_intensities[] <- lapply(post_plasma_intensities, as.numeric)
    post_clinical_data[, covariate] <- ifelse(post_clinical_data[, covariate] == "Male", 1, 0)

    # Initialize matrices to store results
    baseline_p.val_matrix <- matrix(nrow = nrow(post_plasma_intensities), ncol = length(clin_var))
    baseline_estimate_matrix <- matrix(nrow = nrow(post_plasma_intensities), ncol = length(clin_var))

    for (k in 1:length(clin_var)) {
        if (!clin_var[k] %in% colnames(post_clinical_data)) stop(paste("Clinical variable", clin_var[k], "not found in post_clinical_data."))

        for (i in 1:nrow(post_plasma_intensities)) {
            # Extract the protein data for all samples
            x <- as.numeric(post_plasma_intensities[i, ])

            # Extract the clinical variable data for all samples
            y <- as.numeric(post_clinical_data[, clin_var[k]])
            z <- as.numeric(post_clinical_data[, covariate])

            # Filter valid samples
            valid_idx <- complete.cases(x, y, z)

            # Perform partial Spearman correlation on valid data
            test_result <- tryCatch({
                pcor.test(
                    x = x[valid_idx],
                    y = y[valid_idx],
                    z = z[valid_idx],
                    method = "spearman"
                )
            }, error = function(e) list(p.value = NA, estimate = NA))

            baseline_p.val_matrix[i, k] <- test_result$p.value
            baseline_estimate_matrix[i, k] <- test_result$estimate
        }
    }

    # Adjust p-values
    new_matrix_adjust <- matrix(nrow = nrow(baseline_p.val_matrix), ncol = ncol(baseline_p.val_matrix))
    for (i in 1:ncol(baseline_p.val_matrix)) {
        new_matrix_adjust[, i] <- p.adjust(baseline_p.val_matrix[, i], method = "BH")
    }

    # Set column names
    colnames(baseline_p.val_matrix) <- paste0(clin_var, "_p.value")
    colnames(baseline_estimate_matrix) <- paste0(clin_var, "_estimate")
    colnames(new_matrix_adjust) <- paste0(clin_var, "_adj.p.value")

    # Combine adjusted and unadjusted p-values and estimates
    complete_p.val_matrix <- cbind(baseline_p.val_matrix, new_matrix_adjust)
    complete_cor_matrix <- cbind(complete_p.val_matrix, baseline_estimate_matrix)

    # Reorganize columns: p.value, adj.p.value, estimate for each clinical variable
    organized_columns <- unlist(lapply(clin_var, function(var) {
        c(paste0(var, "_p.value"), paste0(var, "_adj.p.value"), paste0(var, "_estimate"))
    }))
    complete_cor_matrix <- complete_cor_matrix[, organized_columns]

    # Create data frame with results
    result_frame <- data.frame(
        ID = rownames(post_plasma_intensities),
        complete_cor_matrix,
        row.names = NULL
    )

    # Return the result
    return(result_frame)
}


# Run for each workflow
MagNet_partial_cor_frame_post <- perform_partial_correlation_post(
    processed_datasets$MagNet_discovery_cohort$corrected,
    post_clinical_data,
    clin_var,
    "Gender"
)

Neat_partial_cor_frame_post <- perform_partial_correlation_post(
    processed_datasets$Neat_discovery_cohort$corrected,
    post_clinical_data,
    clin_var,
    "Gender"
)

PCA_partial_cor_frame_post <- perform_partial_correlation_post(
    processed_datasets$PCA_discovery_cohort$corrected,
    post_clinical_data,
    clin_var,
    "Gender"
)

Depleted_partial_cor_frame_post <- perform_partial_correlation_post(
    processed_datasets$Depleted_discovery_cohort$corrected,
    post_clinical_data,
    clin_var,
    "Gender"
)

Olink_partial_cor_frame_post <- perform_partial_correlation_post(
    processed_datasets$Olink_discovery_cohort[, -1],
    post_clinical_data,
    clin_var,
    "Gender"
)


datasets_partial_clinical_correlation <- list("MagNet_baseline" = MagNet_partial_cor_frame, "MagNet_insulin" = MagNet_partial_cor_frame_post, "Neat_baseline" = Neat_partial_cor_frame, "Neat_insulin" = Neat_partial_cor_frame_post,
                                      "PCA_baseline" = PCA_partial_cor_frame, "PCA_insulin" = PCA_partial_cor_frame_post,"Depleted_baseline" = Depleted_partial_cor_frame, "Depleted_insulin" = Depleted_partial_cor_frame_post,
                                      "Olink_baseline" = Olink_partial_cor_frame, "Olink_insulin" = Olink_partial_cor_frame_post)



saveRDS(datasets_partial_clinical_correlation, file = "data/datasets_partial_clinical_correlation.rds")
