# ==============================================================================
# Tissue and cell-type enrichment of covariate-adjusted M-value-associated
# proteins
#   Extended figure 3A   Tissue enrichment volcano
#   Extended figure 3B   Cell-type enrichment volcano
# Inputs   data-raw/Mvalue_covariate_adjusted_tissue.tsv
#          data-raw/hpa_24.tsv
# ==============================================================================

# Load necessary libraries
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(stringr)

# Load datasets
test_hpa <- readr::read_tsv("data-raw/Mvalue_covariate_adjusted_tissue.tsv")

hpa_data <- readr::read_tsv("data-raw/hpa_24.tsv")

# Separate multiple Uniprot IDs into individual rows
mult_uniprot <- hpa_data %>%
    filter(stringr::str_detect(Uniprot, ", ")) %>%
    tidyr::separate_rows(Uniprot, sep = ", ")

hpa_data <- dplyr::bind_rows(hpa_data %>% filter(!stringr::str_detect(Uniprot, ", ")), mult_uniprot)

# Rename and merge datasets
test_hpa <- test_hpa %>%
    dplyr::rename("Uniprot" = "Gene_ID")

colnames(hpa_data) <- gsub(" ", ".", colnames(hpa_data))

merged_data <- test_hpa %>%
    dplyr::left_join(hpa_data, by = "Uniprot")


# Enrichment analysis for cell types and tissues, using UP vs DOWN in Regulation
enrichment_analysis_with_proteins <- function(data, specificity_column, specificity_label, regulation_type) {
    data <- data %>%
        filter(!is.na(.data[[specificity_column]])) %>%  # Exclude NA values
        mutate(is_associated = ifelse(Regulation == regulation_type, 1, 0))

    # Parse and expand specificity
    specificity_data <- data %>%
        mutate(specificity_split = strsplit(.data[[specificity_column]], ";")) %>%
        unnest(cols = specificity_split) %>%
        separate(specificity_split, into = c(specificity_label, "nTPM"), sep = ": ") %>%
        mutate(nTPM = as.numeric(nTPM)) %>%
        filter(!is.na(.data[[specificity_label]]))  # Exclude NA after splitting

    # Calculate totals for Fisher's test
    total_associated <- sum(data$is_associated == 1, na.rm = TRUE)
    total_not_associated <- sum(data$is_associated == 0, na.rm = TRUE)

    # Perform enrichment analysis
    enrichment_results <- specificity_data %>%
        group_by(.data[[specificity_label]]) %>%
        summarise(
            associated_count = sum(is_associated == 1, na.rm = TRUE),
            not_associated_count = sum(is_associated == 0, na.rm = TRUE),
            associated_not_in_specificity = total_associated - associated_count,
            not_associated_not_in_specificity = total_not_associated - not_associated_count,
            proteins = list(Uniprot[is_associated == 1])
        ) %>%
        rowwise() %>%
        mutate(
            fisher_p = fisher.test(matrix(c(associated_count, not_associated_count,
                                            associated_not_in_specificity, not_associated_not_in_specificity),
                                          ncol = 2))$p.value,
            fisher_OR = fisher.test(matrix(c(associated_count, not_associated_count,
                                             associated_not_in_specificity, not_associated_not_in_specificity),
                                           ncol = 2))$estimate
        ) %>%
        ungroup() %>%
        mutate(
            padj = p.adjust(fisher_p, method = "BH"),
            specificity_label = specificity_label,
            regulation = regulation_type
        )

    return(enrichment_results)
}

# Perform enrichment analysis separately for "UP" and "DOWN" regulation
celltype_enrichment_up <- enrichment_analysis_with_proteins(merged_data,
                                                            "RNA.single.cell.type.specific.nTPM",
                                                            "celltype",
                                                            "UP")

celltype_enrichment_down <- enrichment_analysis_with_proteins(merged_data,
                                                              "RNA.single.cell.type.specific.nTPM",
                                                              "celltype",
                                                              "DOWN")

tissue_enrichment_up <- enrichment_analysis_with_proteins(merged_data,
                                                          "RNA.tissue.specific.nTPM",
                                                          "tissue",
                                                          "UP")

tissue_enrichment_down <- enrichment_analysis_with_proteins(merged_data,
                                                            "RNA.tissue.specific.nTPM",
                                                            "tissue",
                                                            "DOWN")

# Combine enrichment results
celltype_enrichment_combined <- bind_rows(celltype_enrichment_up, celltype_enrichment_down) %>%
    filter(!is.na(celltype))  # Exclude NA values in celltype

tissue_enrichment_combined <- bind_rows(tissue_enrichment_up, tissue_enrichment_down) %>%
    filter(!is.na(tissue))  # Exclude NA values in tissue


# ==============================================================================
# Extended figure 3A - tissue enrichment volcano
# ==============================================================================

# Define tissues to exclude
tissues_to_exclude <- c("testis", "placenta", "vagina", "prostate", "fallopian tube",
                        "seminal vesicle", "ovary", "endometrium 1", "epididymis", "cervix")

# Filter out unwanted tissues
filtered_data <- tissue_enrichment_combined %>%
    filter(!tissue %in% tissues_to_exclude)

# Prepare data for plotting
plot_data <- filtered_data %>%
    mutate(
        log2_OR = log2(fisher_OR),   # Log-transform the odds ratio
        neg_log10_p = -log10(padj),  # Compute -log10(p-adj) for significance
        shape = ifelse(regulation == "UP", 24, 21)  # Triangle for UP, Circle for DOWN
    )

# Plot using shape for regulation and color for significance
ggplot(plot_data[!plot_data$log2_OR == "-Inf",], aes(x = log2_OR, y = neg_log10_p)) +
    geom_point(aes(shape = regulation, fill = neg_log10_p), size = 5, alpha = 0.8, color = "black") +  # Shape and fill
    geom_text(aes(label = tissue), vjust = -0.5, hjust = 0.5, size = 4, check_overlap = TRUE) +  # Add tissue labels
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") + # Reference line at OR = 1 (log2(OR) = 0)
    scale_shape_manual(values = c("UP" = 24, "DOWN" = 21)) +  # Triangle for UP, Circle for DOWN
    scale_fill_gradient(low = "lightblue", high = "darkred") +  # Color by significance (-log10 p-adj)
    labs(
        x = "Log2(Odds Ratio) (Enrichment vs. Depletion)",
        y = "-log10(p-adj value)",
        title = "Tissue-Specific Protein Enrichment & Depletion",
        shape = "Regulation",
        fill = "-log10(p-adj)"
    ) +
    theme_minimal() +
    theme(
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.title.x = element_text(color = "black"),
        axis.title.y = element_text(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")
    )


# ==============================================================================
# Extended figure 3B - cell-type enrichment volcano
# ==============================================================================

plot_data <- celltype_enrichment_combined %>%
    mutate(
        log2_OR = log2(fisher_OR),   # Log-transform the odds ratio
        neg_log10_p = -log10(padj),  # Compute -log10(p-adj) for significance
        shape = ifelse(regulation == "UP", 24, 21)  # Triangle for UP, Circle for DOWN
    )

# Plot using shape for regulation and color for significance
ggplot(plot_data[!plot_data$log2_OR == "-Inf",], aes(x = log2_OR, y = neg_log10_p)) +
    geom_point(aes(shape = regulation, fill = neg_log10_p), size = 5, alpha = 0.8, color = "black") +  # Shape and fill
    geom_text(aes(label = celltype), vjust = -0.5, hjust = 0.5, size = 4, check_overlap = TRUE) +  # Add tissue labels
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") + # Reference line at OR = 1 (log2(OR) = 0)
    scale_shape_manual(values = c("UP" = 24, "DOWN" = 21)) +  # Triangle for UP, Circle for DOWN
    scale_fill_gradient(low = "lightblue", high = "darkred") +  # Color by significance (-log10 p-adj)
    labs(
        x = "Log2(Odds Ratio) (Enrichment vs. Depletion)",
        y = "-log10(p-adj value)",
        title = "Cell-Specific Protein Enrichment & Depletion",
        shape = "Regulation",
        fill = "-log10(p-adj)"
    ) +
    theme_minimal() +
    theme(
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.title.x = element_text(color = "black"),
        axis.title.y = element_text(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")
    )
