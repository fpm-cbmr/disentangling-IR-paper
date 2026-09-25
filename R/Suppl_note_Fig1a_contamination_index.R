# ==============================================================================
# Erythrocyte and platelet contamination indices
# Contamination index = fraction of total (non-log) MS signal carried by a panel
# of erythrocyte- or platelet-specific proteins, z-scored across samples.
# Shown for the neat (undepleted) plasma workflow, discovery cohort (n = 150).

library(dplyr)
library(tidyr)
library(ggplot2)

WORKFLOW  <- "Neat_discovery_cohort"
Z_CUTOFF  <- 2


# ---- 1. Contamination marker panels ------------------------------------------

# Gene.names fields may hold several semicolon-separated synonyms
split_semicolon <- function(x) {
    x <- x[!is.na(x)]
    unique(trimws(unlist(strsplit(x, ";", fixed = TRUE))))
}

erythrocyte_list <- read.delim("data-raw/Eythrocyte_contamination_list.txt")
platelet_list    <- read.delim("data-raw/Platelet_contamination_list.txt")

# Excluded: not specific to the compartment in plasma
erythrocyte_list <- erythrocyte_list[!grepl("NIF3L1", erythrocyte_list$Gene.names), ]
platelet_list    <- platelet_list[!grepl("F13A1",  platelet_list$Gene.names), ]

ery_genes <- split_semicolon(erythrocyte_list$Gene.names)
plt_genes <- split_semicolon(platelet_list$Gene.names)


# ---- 2. Expression matrix ----------------------------------------------------

processed_datasets <- readRDS("data/processed_datasets.rds")

expr_mat <- as.matrix(processed_datasets[[WORKFLOW]]$raw)
expr_mat[is.nan(expr_mat)] <- NA

# rownames are "UniProtID_GENE"
gene_of_row <- sub("^[^_]+_", "", rownames(expr_mat))

ery_hits <- rownames(expr_mat)[gene_of_row %in% ery_genes]
plt_hits <- rownames(expr_mat)[gene_of_row %in% plt_genes]

message("Erythrocyte panel: ", length(ery_hits), " features matched of ",
        length(ery_genes), " genes")
message("Platelet panel:    ", length(plt_hits), " features matched of ",
        length(plt_genes), " genes")


# ---- 3. Contamination index --------------------------------------------------

contamination_index <- function(expr_mat, panel_features) {
    panel_features <- intersect(panel_features, rownames(expr_mat))
    if (length(panel_features) == 0) stop("No panel features found in expression matrix.")
    colSums(expr_mat[panel_features, , drop = FALSE], na.rm = TRUE) /
        colSums(expr_mat, na.rm = TRUE)
}

contamination_scores <- data.frame(
    Sample            = colnames(expr_mat),
    Erythrocyte_score = contamination_index(expr_mat, ery_hits),
    Platelet_score    = contamination_index(expr_mat, plt_hits),
    stringsAsFactors  = FALSE
) %>%
    mutate(
        Erythrocyte_z = as.numeric(scale(Erythrocyte_score)),
        Platelet_z    = as.numeric(scale(Platelet_score))
    )


# ---- 4. Figure ---------------------------------------------------------------

panel_levels <- c("Erythrocyte contamination", "Platelet contamination")

plot_df <- contamination_scores %>%
    dplyr::select(Sample, Erythrocyte_z, Platelet_z) %>%
    pivot_longer(
        cols      = c(Erythrocyte_z, Platelet_z),
        names_to  = "Panel",
        values_to = "Zscore"
    ) %>%
    mutate(
        Panel  = factor(c(Erythrocyte_z = panel_levels[1],
                          Platelet_z    = panel_levels[2])[Panel],
                        levels = panel_levels),
        Sample = factor(Sample, levels = colnames(expr_mat))
    )

p <- ggplot(plot_df, aes(x = Sample, y = Zscore)) +
    geom_point(
        shape  = 21,
        size   = 2.3,
        alpha  = 0.4,
        fill   = "grey70",
        colour = "black",
        stroke = 0.4
    ) +
    geom_hline(yintercept = Z_CUTOFF, linetype = "dashed", linewidth = 0.25) +
    facet_wrap(~ Panel, ncol = 1, scales = "free_y") +
    theme_classic() +
    theme(
        axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank(),
        axis.text.y      = element_text(colour = "black"),
        axis.title       = element_text(colour = "black"),
        strip.text       = element_text(colour = "black"),
        strip.background = element_blank(),
        panel.spacing    = grid::unit(1, "cm")
    ) +
    labs(x = "Samples", y = "z-score")

