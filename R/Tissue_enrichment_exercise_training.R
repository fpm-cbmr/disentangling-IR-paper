suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
    library(stringr)
    library(ggplot2)
    library(forcats)
})

hpa_file        <- "data-raw/hpa_24.tsv"
enrichment_file <- "data-raw/HIIT_training_enrichment.txt"

## ---------- Input: training-responsive proteins with direction -------------
if (file.exists(enrichment_file)) {
    test_hpa <- read.delim(enrichment_file)
} else if (exists("protein_df")) {
    test_hpa <- protein_df
} else {
    stop("Run Validation_training_script_part2_Fig4.R first, or provide 'protein_df'.")
}

# Standardize the join key
test_hpa <- test_hpa %>%
    dplyr::rename(Uniprot = UniProt_ID) %>%
    mutate(Uniprot = toupper(trimws(Uniprot)))

## ---------- Load & prep HPA -----------------------------------------------
hpa_data <- read_tsv(hpa_file, show_col_types = FALSE)
colnames(hpa_data) <- gsub(" ", ".", colnames(hpa_data))
if (!("Uniprot" %in% names(hpa_data))) stop("HPA file must have column 'Uniprot'.")

# Split rows carrying several UniProt IDs
mult_uniprot <- hpa_data %>%
    filter(str_detect(Uniprot, ", ")) %>%
    separate_rows(Uniprot, sep = ", ")

hpa_data <- bind_rows(hpa_data %>% filter(!str_detect(Uniprot, ", ")), mult_uniprot) %>%
    mutate(Uniprot = toupper(trimws(Uniprot)))

## ---------- Merge ----------------------------------------------------------
merged_data <- test_hpa %>% left_join(hpa_data, by = "Uniprot")

## ---------- Fisher's exact test, one direction at a time -------------------
enrichment_one_dir <- function(data, specificity_column, specificity_label, desired_direction) {
    base <- data %>%
        filter(!is.na(.data[[specificity_column]])) %>%
        mutate(is_associated = if_else(Is_Significant == "Yes" & Direction == desired_direction, 1L, 0L))

    if (nrow(base) == 0) {
        return(tibble(
            !!specificity_label := character(),
            associated_count = integer(), not_associated_count = integer(),
            associated_not_in_specificity = integer(), not_associated_not_in_specificity = integer(),
            fisher_p = numeric(), fisher_OR = numeric(), padj = numeric(),
            regulation = character(), proteins = list()
        ))
    }

    spec_df <- base %>%
        mutate(specificity_split = strsplit(.data[[specificity_column]], ";")) %>%
        unnest(cols = specificity_split) %>%
        separate(specificity_split, into = c({{ specificity_label }}, "nTPM"), sep = ": ") %>%
        mutate(nTPM = suppressWarnings(as.numeric(nTPM))) %>%
        filter(!is.na(.data[[specificity_label]])) %>%
        distinct(Uniprot, .data[[specificity_label]], .keep_all = TRUE)

    if (nrow(spec_df) == 0) {
        return(tibble(
            !!specificity_label := character(),
            associated_count = integer(), not_associated_count = integer(),
            associated_not_in_specificity = integer(), not_associated_not_in_specificity = integer(),
            fisher_p = numeric(), fisher_OR = numeric(), padj = numeric(),
            regulation = character(), proteins = list()
        ))
    }

    total_assoc     <- n_distinct(spec_df$Uniprot[spec_df$is_associated == 1])
    total_not_assoc <- n_distinct(spec_df$Uniprot[spec_df$is_associated == 0])

    spec_df %>%
        group_by(.data[[specificity_label]]) %>%
        summarise(
            associated_count     = n_distinct(Uniprot[is_associated == 1]),
            not_associated_count = n_distinct(Uniprot[is_associated == 0]),
            proteins             = list(unique(Uniprot[is_associated == 1])),
            .groups = "drop"
        ) %>%
        mutate(
            associated_not_in_specificity     = total_assoc     - associated_count,
            not_associated_not_in_specificity = total_not_assoc - not_associated_count
        ) %>%
        rowwise() %>%
        mutate(
            fisher_p  = fisher.test(matrix(
                c(associated_count, not_associated_count,
                  associated_not_in_specificity, not_associated_not_in_specificity),
                ncol = 2))$p.value,
            fisher_OR = as.numeric(fisher.test(matrix(
                c(associated_count, not_associated_count,
                  associated_not_in_specificity, not_associated_not_in_specificity),
                ncol = 2))$estimate)
        ) %>%
        ungroup() %>%
        mutate(
            padj       = p.adjust(fisher_p, method = "BH"),
            regulation = desired_direction
        )
}

## ---------- Run both directions and combine -------------------------------
tissue_up   <- enrichment_one_dir(merged_data, "RNA.tissue.specific.nTPM", "tissue", "Up")
tissue_down <- enrichment_one_dir(merged_data, "RNA.tissue.specific.nTPM", "tissue", "Down")
tissue_enrichment_combined <- bind_rows(tissue_up, tissue_down) %>%
    filter(!is.na(tissue))

celltype_up   <- enrichment_one_dir(merged_data, "RNA.single.cell.type.specific.nTPM", "celltype", "Up")
celltype_down <- enrichment_one_dir(merged_data, "RNA.single.cell.type.specific.nTPM", "celltype", "Down")
celltype_enrichment_combined <- bind_rows(celltype_up, celltype_down) %>%
    filter(!is.na(celltype))

## ---------- Volcano: shape encodes direction ------------------------------
volcano_both <- function(enrichment_data, label_col, title_prefix = "", exclude_labels = NULL) {
    df <- enrichment_data %>%
        { if (!is.null(exclude_labels)) filter(., !tolower(.data[[label_col]]) %in% tolower(exclude_labels)) else . } %>%
        mutate(
            log2_OR     = log2(fisher_OR),
            neg_log10_p = -log10(padj)
        ) %>%
        filter(is.finite(log2_OR))

    ggplot(df, aes(x = log2_OR, y = neg_log10_p)) +
        geom_point(aes(shape = regulation, fill = neg_log10_p), size = 5, alpha = 0.85, color = "black") +
        geom_text(aes(label = .data[[label_col]]), vjust = -0.5, size = 3.4, check_overlap = TRUE) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
        scale_shape_manual(values = c(Up = 24, Down = 21)) +   # triangle = Up, circle = Down
        scale_fill_gradient(low = "lightblue", high = "darkred") +
        labs(
            x     = "Log2(Odds Ratio) (Enrichment vs. Depletion)",
            y     = "-log10(FDR)",
            title = paste0(title_prefix, " (UP vs DOWN)"),
            shape = "Direction",
            fill  = "-log10(FDR)"
        ) +
        theme_minimal() +
        theme(
            axis.text        = element_text(color = "black"),
            axis.title       = element_text(color = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line        = element_line(color = "black")
        )
}

# Reproductive and sex-specific tissues are excluded from the tissue volcano
tissues_to_exclude <- c("testis", "placenta", "vagina", "prostate", "fallopian tube",
                        "seminal vesicle", "ovary", "endometrium 1", "epididymis", "cervix")

## MAIN FIGURE, PANEL d -----------------------------------------------------
tissue_volcano <- volcano_both(tissue_enrichment_combined, "tissue",
                               "Tissue Enrichment", exclude_labels = tissues_to_exclude)
tissue_volcano

## Cell-type-resolved equivalent (not a published panel) --------------------
celltype_volcano <- volcano_both(celltype_enrichment_combined, "celltype",
                                 "Cell-type Enrichment")
celltype_volcano
