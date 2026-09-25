# ==============================================================================
# Tissue and cell-type enrichment of M-value-associated plasma proteins
#   Figure a        Tissue origin map        -> m_value_tissue_heatmap.pdf
#   Extended fig c  Cell-type volcano        -> printed to the device
#
# Extended figure a and b are not produced by this script - see
# Tissue_celltype_enrichment_related_to_extended_Fig3.R
# Inputs   data-raw/M_value_enrichment.txt
#          data-raw/hpa_24.tsv
# ==============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(ggbreak)
library(ComplexHeatmap)
library(circlize)

# Load datasets
test_hpa <- readr::read_tsv("data-raw/M_value_enrichment.txt")
hpa_data <- readr::read_tsv("data-raw/hpa_24.tsv")

# Separate multiple Uniprot IDs into individual rows
mult_uniprot <- hpa_data %>%
    filter(stringr::str_detect(Uniprot, ", ")) %>%
    tidyr::separate_rows(Uniprot, sep = ", ")

hpa_data <- dplyr::bind_rows(hpa_data %>% filter(!stringr::str_detect(Uniprot, ", ")), mult_uniprot)

# rename uniprot
test_hpa <- test_hpa %>%
    dplyr::rename("Uniprot" = "UniProt_ID")

colnames(hpa_data) <- gsub(" ", ".", colnames(hpa_data))

# Merge datasets
merged_data <- test_hpa %>%
    dplyr::left_join(hpa_data, by = c("Uniprot" = "Uniprot"))


# Define protein groups
group_lists <- list(
    m_positive = merged_data %>%
        dplyr::filter((Is_Significant == "Yes" & Direction == "Up")) %>%
        dplyr::pull(GeneID),

    m_negative = merged_data %>%
        dplyr::filter((Is_Significant == "Yes" & Direction == "Down")) %>%
        dplyr::pull(GeneID),

    m_all = merged_data %>%
        dplyr::filter((Is_Significant == "Yes")) %>%
        dplyr::pull(GeneID)
)

m_positive <- group_lists$m_positive
m_negative <- group_lists$m_negative
m_all <- group_lists$m_all

association_colors <- c(
    "m_positive" = "#882255",
    "m_negative"   = "#88CCEE",
    "m_all"   = "#332288"
)

# Replace NA values with 0 for numeric columns
merged_data <- merged_data %>%
    mutate(across(where(is.numeric), ~ replace_na(., 0)))

new_lists <- list(
    m_positive = m_positive,
    m_negative = m_negative,
    m_all = m_all
)

# Add group annotations
merged_data <- merged_data %>%
    mutate(
        m_positive = ifelse(GeneID %in% new_lists$m_positive, 1, 0),
        m_negative = ifelse(GeneID %in% new_lists$m_negative, 1, 0),
        m_all = ifelse(GeneID %in% new_lists$m_all, 1, 0)
    )


### cell type enrichment analysis ###
results_list <- lapply(names(new_lists), function(group_name) {
    group_data <- merged_data %>%
        mutate(is_associated = get(group_name))

    celltype_data <- group_data %>%
        mutate(RNA_celltype_split = strsplit(RNA.single.cell.type.specific.nTPM, ";")) %>%
        unnest(cols = RNA_celltype_split) %>%
        separate(RNA_celltype_split, into = c("celltype", "nTPM"), sep = ": ") %>%
        mutate(nTPM = as.numeric(nTPM))

    total_associated <- sum(group_data$is_associated == 1, na.rm = TRUE)
    total_not_associated <- sum(group_data$is_associated == 0, na.rm = TRUE)

    celltype_enrichment <- celltype_data %>%
        group_by(celltype) %>%
        summarise(
            associated_count = sum(is_associated == 1, na.rm = TRUE),
            not_associated_count = sum(is_associated == 0, na.rm = TRUE),
            associated_not_in_celltype = total_associated - associated_count,
            not_associated_not_in_celltype = total_not_associated - not_associated_count
        ) %>%
        rowwise() %>%
        mutate(
            fisher_p = fisher.test(
                matrix(c(associated_count, not_associated_count,
                         associated_not_in_celltype, not_associated_not_in_celltype),
                       ncol = 2)
            )$p.value,
            fisher_OR = fisher.test(
                matrix(c(associated_count, not_associated_count,
                         associated_not_in_celltype, not_associated_not_in_celltype),
                       ncol = 2)
            )$estimate
        ) %>%
        ungroup() %>%
        mutate(
            padj = p.adjust(fisher_p, method = "BH"),
            group = group_name
        )

    write_csv(celltype_enrichment, paste0("enrichment_celltype_", group_name, ".csv"))

    return(celltype_enrichment)
})

combined_cell_results <- do.call(rbind, results_list)


### tissue enrichment
tissue_results_list <- lapply(names(new_lists), function(group_name) {
    group_data <- merged_data %>%
        mutate(is_associated = ifelse(GeneID %in% new_lists[[group_name]], 1, 0))

    tissue_data <- group_data %>%
        mutate(RNA_tissue_split = strsplit(RNA.tissue.specific.nTPM, ";")) %>%
        unnest(cols = RNA_tissue_split) %>%
        separate(RNA_tissue_split, into = c("tissue", "nTPM"), sep = ": ") %>%
        mutate(nTPM = as.numeric(nTPM))

    total_associated <- sum(group_data$is_associated == 1, na.rm = TRUE)
    total_not_associated <- sum(group_data$is_associated == 0, na.rm = TRUE)

    tissue_enrichment <- tissue_data %>%
        group_by(tissue) %>%
        summarise(
            associated_count = sum(is_associated == 1, na.rm = TRUE),
            not_associated_count = sum(is_associated == 0, na.rm = TRUE),
            associated_not_in_tissue = total_associated - associated_count,
            not_associated_not_in_tissue = total_not_associated - not_associated_count
        ) %>%
        rowwise() %>%
        mutate(
            fisher_p = fisher.test(
                matrix(c(associated_count, not_associated_count,
                         associated_not_in_tissue, not_associated_not_in_tissue),
                       ncol = 2)
            )$p.value,
            fisher_OR = fisher.test(
                matrix(c(associated_count, not_associated_count,
                         associated_not_in_tissue, not_associated_not_in_tissue),
                       ncol = 2)
            )$estimate
        ) %>%
        ungroup() %>%
        mutate(
            padj = p.adjust(fisher_p, method = "BH"),
            group = group_name
        )

    write_csv(tissue_enrichment, paste0("enrichment_tissue_", group_name, ".csv"))

    return(tissue_enrichment)
})

combined_tissue_results <- do.call(rbind, tissue_results_list)


### build the OR / padj matrices
prepare_heatmap_data_cell <- function(data) {
    logOR_matrix <- data %>%
        dplyr::select(celltype, group, fisher_OR) %>%
        pivot_wider(names_from = group, values_from = fisher_OR, values_fill = 0) %>%
        tibble::column_to_rownames("celltype") %>%
        as.matrix()

    pvalue_matrix <- data %>%
        dplyr::select(celltype, group, padj) %>%
        pivot_wider(names_from = group, values_from = padj, values_fill = Inf) %>%
        tibble::column_to_rownames("celltype") %>%
        as.matrix()

    list(logOR_matrix = logOR_matrix, pvalue_matrix = pvalue_matrix)
}

prepare_heatmap_data_tissue <- function(data) {
    logOR_matrix <- data %>%
        dplyr::select(tissue, group, fisher_OR) %>%
        pivot_wider(names_from = group, values_from = fisher_OR, values_fill = 0) %>%
        tibble::column_to_rownames("tissue") %>%
        as.matrix()

    pvalue_matrix <- data %>%
        dplyr::select(tissue, group, padj) %>%
        pivot_wider(names_from = group, values_from = padj, values_fill = Inf) %>%
        tibble::column_to_rownames("tissue") %>%
        as.matrix()

    list(logOR_matrix = logOR_matrix, pvalue_matrix = pvalue_matrix)
}

heatmap_cell   <- prepare_heatmap_data_cell(combined_cell_results %>% filter(!is.na(celltype)))
heatmap_tissue <- prepare_heatmap_data_tissue(combined_tissue_results %>% filter(!is.na(tissue)))


# ==============================================================================
# Figure a - tissue origin map
# ==============================================================================

plot_heatmap_tissue <- function(
        logOR_matrix,
        pvalue_matrix,
        title,
        output_file,
        association_colors
) {
    # Replace zeros with NA to handle log2 transformation
    logOR_matrix[logOR_matrix == 0] <- NA

    # Apply log2 transformation
    log2OR_matrix <- log2(logOR_matrix)

    # Create a temporary matrix with NA replaced by 0 for clustering
    log2OR_matrix_temp <- log2OR_matrix
    log2OR_matrix_temp[is.na(log2OR_matrix_temp)] <- 0

    # Perform clustering using the temporary matrix
    row_dendrogram <- hclust(dist(log2OR_matrix_temp))
    col_dendrogram <- hclust(dist(t(log2OR_matrix_temp)))

    # Restore original NA matrix for plotting
    log2OR_matrix_plot <- log2OR_matrix

    # Define column annotations
    col_annotation <- HeatmapAnnotation(
        Group = colnames(log2OR_matrix_plot),
        col = list(Group = association_colors)
    )

    # Define the heatmap
    heatmap <- Heatmap(
        log2OR_matrix_plot,
        name = "Log2(OR)",
        col = colorRamp2(c(-2, 0, 2), c("#006ae3", "white", "#d70000")),
        na_col = "grey",
        top_annotation = col_annotation,
        column_names_gp = gpar(fontsize = 0),
        show_column_names = FALSE,
        cluster_rows = as.dendrogram(row_dendrogram),
        cluster_columns = as.dendrogram(col_dendrogram),
        cell_fun = function(j, i, x, y, width, height, fill) {
            if (!is.na(pvalue_matrix[i, j]) && pvalue_matrix[i, j] < 0.05) {
                grid.text("*", x, y, gp = gpar(fontsize = 16, col = "black"))
            }
        },
        heatmap_legend_param = list(
            title_position = "topcenter",
            title_gp = gpar(fontsize = 10, fontface = "bold")
        )
    )

    # Save and display the heatmap
    tryCatch({
        pdf(output_file, width = 8, height = 12)
        draw(
            heatmap,
            heatmap_legend_side = "right",
            annotation_legend_side = "right",
            merge_legends = TRUE
        )
        dev.off()
        message("Heatmap saved as: ", output_file)
    }, error = function(e) {
        message("Error generating heatmap: ", e$message)
        if (dev.cur() != 1) dev.off()
    })

    grid.newpage()
    draw(
        heatmap,
        heatmap_legend_side = "right",
        annotation_legend_side = "right",
        merge_legends = TRUE
    )
}

# Exclude specified tissues
tissues_to_exclude <- c("testis", "placenta", "vagina","prostate","fallopian tube","seminal vesicle","ovary","endometrium 1","epididymis","cervix")

filtered_logOR_matrix <- heatmap_tissue$logOR_matrix[!rownames(heatmap_tissue$logOR_matrix) %in% tissues_to_exclude, , drop = FALSE]
filtered_pvalue_matrix <- heatmap_tissue$pvalue_matrix[!rownames(heatmap_tissue$pvalue_matrix) %in% tissues_to_exclude, , drop = FALSE]

plot_heatmap_tissue(
    logOR_matrix = filtered_logOR_matrix,
    pvalue_matrix = filtered_pvalue_matrix,
    title = "Tissue Enrichment Across Groups (Filtered)",
    output_file = "m_value_tissue_heatmap.pdf",
    association_colors = association_colors
)


# ==============================================================================
# Extended figure c - cell-resolved volcano
# ==============================================================================

# NOTE: the column indices below assume the matrix columns are in the order
# m_positive, m_negative, m_all. That holds because new_lists is defined in that
# order, but it is positional - check if you ever reorder the groups.
Comb_cell_direction <- as.data.frame(heatmap_cell$pvalue_matrix)
Comb_cell_direction$celltype <- rownames(heatmap_cell$pvalue_matrix)
colnames(Comb_cell_direction)[1:3] <- c("p.val_positive","p.val_negative","p.val_all")
Comb_cell_direction <- cbind(Comb_cell_direction,heatmap_cell$logOR_matrix)
colnames(Comb_cell_direction)[5:7] <- c("OR.positive","OR.negative","OR.all")

# Reshape the data to long format
Comb_cell_direction_long <- Comb_cell_direction %>%
    dplyr::select(celltype, p.val_positive, p.val_negative, OR.positive, OR.negative) %>%
    pivot_longer(cols = starts_with("p.val"),
                 names_to = "direction",
                 values_to = "p_value") %>%
    mutate(OR_value = ifelse(direction == "p.val_positive", OR.positive, OR.negative),
           direction = ifelse(direction == "p.val_positive", "positive", "negative")) %>%
    mutate(log2_OR = log2(OR_value),
           neg_log10_p = -log10(p_value))

significance_threshold <- -log10(0.05)

# Remove rows with -Inf Log2(OR)
Comb_cell_direction_long <- Comb_cell_direction_long %>%
    filter(log2_OR != -Inf)

# Identify top 5 most significant for each direction and Log2(OR) sign
top_positive_posOR <- Comb_cell_direction_long %>%
    filter(direction == "positive", log2_OR > 0) %>%
    arrange(p_value) %>%
    dplyr::slice(1:5)

top_positive_negOR <- Comb_cell_direction_long %>%
    filter(direction == "positive", log2_OR < 0) %>%
    arrange(p_value) %>%
    dplyr::slice(1:5)

top_negative_posOR <- Comb_cell_direction_long %>%
    filter(direction == "negative", log2_OR > 0) %>%
    arrange(p_value) %>%
    dplyr::slice(1:5)

top_negative_negOR <- Comb_cell_direction_long %>%
    filter(direction == "negative", log2_OR < 0) %>%
    arrange(p_value) %>%
    dplyr::slice(1:5)

Comb_cell_direction_long <- Comb_cell_direction_long %>%
    mutate(Label_ID = ifelse(celltype %in% c(
        top_positive_posOR$celltype, top_positive_negOR$celltype,
        top_negative_posOR$celltype, top_negative_negOR$celltype),
        "Yes", ""
    ))

# Create the volcano plot
ggplot(Comb_cell_direction_long, aes(x = log2_OR, y = neg_log10_p, shape = direction, color = direction)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_text(data = Comb_cell_direction_long %>% filter(Label_ID == "Yes"),
              aes(label = celltype), hjust = 0.5, vjust = -0.5, size = 4, color = "black") +
    scale_shape_manual(values = c(16, 17)) +
    scale_color_manual(values = c("steelblue", "darkred")) +
    geom_hline(yintercept = significance_threshold, linetype = "dashed", color = "black", size = 1) +
    labs(x = "Log2(Odds Ratio)", y = "-Log10(adj. P-value)", title = NULL) +
    scale_y_break(c(13, 30), space = 1) +
    expand_limits(y = 34) +
    xlim(-3, 3) +
    guides(color = "none", shape = "none") +
    theme_classic() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        axis.line.x.top = element_blank(),
        axis.line.x.bottom = element_line(color = "black"),
        axis.ticks.x.top = element_blank(),
        axis.text.x.top = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black")
    )


