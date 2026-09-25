# ==============================================================================
# Differential protein regulation - discovery cohort
# ==============================================================================

library(limma)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

WORKFLOWS <- c("MagNet", "PCA", "Neat", "Depleted", "Olink")

WORKFLOW_COLOURS <- c(
    Depleted = "darkred",
    MagNet   = "darkgreen",
    Neat     = "darkgrey",
    Olink    = "darkorange",
    PCA      = "darkblue"
)

RESULT_DIR <- "data/Differential_regulation"


# ==============================================================================
# 1. Preprocessing, model fitting, and saving results
# ==============================================================================

processed_datasets  <- readRDS("data/processed_datasets.rds")
Sample_ID_filt      <- readRDS("data/Sample_ID_filt.rds")
Clinical_data_filt   <- readRDS("data/Clinical_data_filt")

rownames(processed_datasets$Olink_discovery_cohort) <-
    processed_datasets$Olink_discovery_cohort$ID

# Expression matrix per workflow. Olink is a data frame whose first column is
# the ID, so that column is dropped; the MS workflows use the batch-corrected,
# median-scaled slot.
expr_of <- function(wf) {
    if (wf == "Olink") {
        processed_datasets$Olink_discovery_cohort[, -1]
    } else {
        processed_datasets[[paste0(wf, "_discovery_cohort")]]$corrected
    }
}

Disease <- factor(Sample_ID_filt$disease, levels = c("NGT", "T2D"))
Gender  <- factor(Sample_ID_filt$Gender,  levels = c("Male", "Female"))
Clamp   <- factor(Sample_ID_filt$Clamp,   levels = c("Pre", "Post"))
BMI     <- Clinical_data_filt$BMI..kg.m2.
Age     <- Clinical_data_filt$Age
Subject_ID <- Sample_ID_filt$Subject_ID

design <- model.matrix(~ Disease + BMI + Gender + Age + Clamp)
colnames(design) <- c("Intercept", "Disease", "BMI", "Gender", "Age", "Clamp")

COEF_DISEASE <- which(colnames(design) == "Disease")   # 2
COEF_INSULIN <- which(colnames(design) == "Clamp")     # 6

# Gene/protein ID fields differ between the MS workflows ("UniProt_GENE") and
# Olink ("Assay_Block_UniProt").
add_ids <- function(df, wf) {
    rn <- rownames(df)
    if (wf == "Olink") {
        df$Gene_ID    <- word(rn, 1, sep = "_")
        df$Protein_ID <- word(rn, 3, sep = "_")
    } else {
        df$Gene_ID    <- word(rn, 2, sep = "_")
        df$Protein_ID <- word(rn, 1, sep = "_")
    }
    df
}

# One fit per workflow, each with its OWN duplicateCorrelation consensus.
fits <- list()
for (wf in WORKFLOWS) {
    message("Fitting: ", wf)
    e <- expr_of(wf)
    corfit <- duplicateCorrelation(e, design, block = Subject_ID)
    fits[[wf]] <- eBayes(
        lmFit(e, design, block = Subject_ID, correlation = corfit$consensus)
    )
}

dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

extract_and_save <- function(coef, tag) {
    out <- list()
    for (wf in WORKFLOWS) {
        res <- topTable(fits[[wf]], coef = coef, number = Inf, sort.by = "none")
        res <- add_ids(res, wf)
        out[[wf]] <- res
        saveRDS(res, file.path(
            RESULT_DIR,
            paste0(wf, "_differential_result_", tag, "_BMI_Gender_Clamp_adjusted.rds")
        ))
    }
    out
}

results_GROUP   <- extract_and_save(COEF_DISEASE, "GROUP")
results_INSULIN <- extract_and_save(COEF_INSULIN, "INSULIN")


# ==============================================================================
# 2. Insulin (Pre vs Post clamp) figures
# ==============================================================================

# ---- 2a. Volcano of insulin-regulated proteins -------------------------------

all_results <- bind_rows(lapply(WORKFLOWS, function(wf) {
    transform(results_INSULIN[[wf]], Workflow = wf)
}))

significant_results_Insulin <- subset(all_results, adj.P.Val < 0.05)

message("\nSignificant insulin-regulated proteins (adj.P < 0.05):")
print(table(significant_results_Insulin$Workflow))

# Where a protein is significant in several workflows, keep the best p-value
best_results_insulin <- significant_results_Insulin %>%
    group_by(Protein_ID) %>%
    slice_min(P.Value, n = 1) %>%
    ungroup()

insulin_volcano <- ggplot(
    best_results_insulin,
    aes(x = logFC, y = -log10(P.Value), colour = Workflow, label = Gene_ID)
) +
    geom_point(alpha = 0.5, size = 3) +
    geom_text_repel(colour = "black", box.padding = 0.5, point.padding = 0.5,
                    segment.color = "grey50") +
    scale_colour_manual(values = WORKFLOW_COLOURS) +
    xlim(-1, 2) +
    labs(x = "LogFC (Post-Pre Clamp)", y = "-Log10(p-value)", colour = "Workflow") +
    theme_minimal() +
    theme(legend.position = "right",
          panel.grid  = element_blank(),
          axis.line   = element_line(colour = "black"))

print(insulin_volcano)


# ---- 2b. Insulin logFC against insulin-sensitivity associations --------------
#
# The M-value estimates come from the PARTIAL (covariate-adjusted) clinical
# correlation analysis, not the unadjusted one. Elsewhere this same file is
# loaded into a variable called `datasets_clinical_correlation`, which makes it
# look unadjusted - it is not. Say "partial correlation" in the figure legend.
# The unadjusted data/datasets_clinical_correlation.rds is NOT the right input:
# it uses the old SAX/Undepleted names and has no Depleted entry.

clin_cor <- readRDS("data/datasets_partial_clinical_correlation.rds")

cor_pre <- list(
    MagNet   = clin_cor$MagNet_baseline,
    PCA      = clin_cor$PCA_baseline,
    Neat     = clin_cor$Neat_baseline,
    Depleted = clin_cor$Depleted_baseline,
    Olink    = clin_cor$Olink_baseline
)
cor_post <- list(
    MagNet   = clin_cor$MagNet_insulin,
    PCA      = clin_cor$PCA_insulin,
    Neat     = clin_cor$Neat_insulin,
    Depleted = clin_cor$Depleted_insulin,
    Olink    = clin_cor$Olink_insulin
)

stopifnot(!vapply(c(cor_pre, cor_post), is.null, logical(1)))

prepare_cor <- function(pre_df, post_df, wf) {
    if (is.null(pre_df) || is.null(post_df)) return(NULL)
    id_field <- if (wf == "Olink") 3 else 1
    gene_field <- if (wf == "Olink") 1 else 2
    inner_join(pre_df, post_df, by = "ID", suffix = c(".Pre", ".Post")) %>%
        mutate(
            UniProt_ID = word(ID, id_field,   sep = "_"),
            Gene_ID    = word(ID, gene_field, sep = "_"),
            Workflow   = wf
        )
}

combined_all <- bind_rows(lapply(WORKFLOWS, function(wf)
    prepare_cor(cor_pre[[wf]], cor_post[[wf]], wf)))

# The original joined these two tables by ROW POSITION via cbind(), with the
# workflow blocks concatenated in different orders on each side. Join on the
# protein ID instead, so rows cannot be mis-paired.
insulin_sig <- significant_results_Insulin

LogFC_correlation_insulin <- insulin_sig %>%
    inner_join(
        combined_all %>%
            select(Workflow, UniProt_ID, Gene_ID,
                   Pre  = `M..µmol..kg.min.._estimate.Pre`,
                   Post = `M..µmol..kg.min.._estimate.Post`,
                   adj_Pre  = `M..µmol..kg.min.._adj.p.value.Pre`,
                   adj_Post = `M..µmol..kg.min.._adj.p.value.Post`),
        by = c("Workflow", "Protein_ID" = "UniProt_ID", "Gene_ID")
    ) %>%
    filter(adj_Pre < 0.05 | adj_Post < 0.05)

if (nrow(LogFC_correlation_insulin) == 0) {
    stop("No rows survived the ID join - check the Workflow / Protein_ID / Gene_ID ",
         "fields in both tables before interpreting this figure.")
}
message("Proteins in the overlay figure: ", nrow(LogFC_correlation_insulin))

long_plot_data <- LogFC_correlation_insulin %>%
    select(UniProt_ID = Protein_ID, Gene_ID, logFC, Workflow, Pre, Post) %>%
    pivot_longer(cols = c(Pre, Post),
                 names_to = "Timepoint", values_to = "M_value_estimate") %>%
    mutate(Timepoint = factor(Timepoint, levels = c("Pre", "Post")))

# Exclusions carried over from the original figure
plot_data <- long_plot_data %>%
    filter(UniProt_ID != "P14780", Workflow != "Depleted")

insulin_IS_overlay <- ggplot(
    plot_data,
    aes(x = M_value_estimate, y = logFC, group = UniProt_ID,
        colour = Timepoint, shape = Timepoint)
) +
    geom_point(size = 4, alpha = 0.5, stroke = NA) +
    geom_line(alpha = 0.6, colour = "grey") +
    geom_text(data = ~ subset(., Timepoint == "Post"),
              aes(label = Gene_ID), colour = "black", size = 3.5, hjust = -0.1) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_colour_manual(values = c("darkgrey", "darkred")) +
    scale_shape_manual(values = c(Pre = 16, Post = 17)) +
    labs(x = "M-value estimate", y = "log2 Fold Change (logFC)",
         colour = "Timepoint", shape = "Timepoint") +
    theme_minimal() +
    theme(legend.position = c(1, 1),
          legend.justification = c(1, 1),
          panel.grid  = element_blank(),
          axis.line   = element_line(colour = "black"),
          axis.text   = element_text(colour = "black", size = 14),
          axis.title  = element_text(colour = "black", size = 16),
          axis.ticks  = element_line(colour = "black"))

print(insulin_IS_overlay)


# ==============================================================================
# 3. Limma with M-value 
# ==============================================================================

# The original defined design_Mvalue twice, the second overwriting the first:
#   ~ Mvalue + BMI + Gender + Age + Clamp     (dead)
#   ~ Mvalue + BF  + Gender + Age + Clamp     (the one that ran)
# So the published M-value model is adjusted for BODY FAT %, not BMI. Set this
# to "BMI" only if you intend to change the model.
MVALUE_ADIPOSITY_COVARIATE <- "BF"      # "BF" or "BMI"

Mvalue <- Clinical_data_filt$M..µmol..kg.min..
BF     <- Clinical_data_filt$Body.fat....

na_mv <- which(is.na(Mvalue))
message("\nSamples dropped for missing M-value: ", length(na_mv))

adiposity <- if (MVALUE_ADIPOSITY_COVARIATE == "BF") BF else BMI

design_Mvalue <- model.matrix(~ Mvalue + adiposity + Gender + Age + Clamp)
colnames(design_Mvalue) <- c("Intercept", "Mvalue", MVALUE_ADIPOSITY_COVARIATE,
                             "Gender", "Age", "Clamp")

COEF_MVALUE <- which(colnames(design_Mvalue) == "Mvalue")   # 2

Subject_ID_Mvalue <- Sample_ID_filt[-na_mv, ]$Subject_ID

# Olink's first column is the ID, so sample j sits in column j+1. The original
# hard-coded -c(1, 27); this derives it from the NA positions instead.
expr_of_mv <- function(wf) {
    if (wf == "Olink") {
        processed_datasets$Olink_discovery_cohort[, -c(1, na_mv + 1)]
    } else {
        processed_datasets[[paste0(wf, "_discovery_cohort")]]$corrected[, -na_mv]
    }
}

results_Mvalue <- list()
for (wf in WORKFLOWS) {
    message("Fitting (M-value): ", wf)
    e <- expr_of_mv(wf)
    corfit <- duplicateCorrelation(e, design_Mvalue, block = Subject_ID_Mvalue)
    fit <- eBayes(lmFit(e, design_Mvalue, block = Subject_ID_Mvalue,
                        correlation = corfit$consensus))
    res <- topTable(fit, coef = COEF_MVALUE, number = Inf, sort.by = "none")
    res <- add_ids(res, wf)
    res$Workflow <- wf
    results_Mvalue[[wf]] <- res
}

combined_results_Mvalue <- bind_rows(results_Mvalue)

best_results_Mvalue <- combined_results_Mvalue %>%
    filter(adj.P.Val < 0.05) %>%
    group_by(Protein_ID) %>%
    slice_min(adj.P.Val, n = 1) %>%
    ungroup() %>%
    arrange(logFC) %>%
    mutate(Rank = row_number(),
           PercentageChange = (2^logFC - 1) * 100)

message("\nProteins significantly associated with M-value: ",
        nrow(best_results_Mvalue))

mvalue_theme <- theme_minimal() +
    theme(axis.line  = element_line(colour = "black"),
          panel.grid = element_blank(),
          axis.text  = element_text(colour = "black"),
          axis.title = element_text(colour = "black"))


mvalue_volcano <- ggplot(
    best_results_Mvalue,
    aes(x = logFC, y = -log10(P.Value), colour = Workflow, label = Gene_ID)
) +
    geom_point(size = 4, alpha = 0.5) +
    geom_text_repel(size = 3, box.padding = 0.2, point.padding = 0.2,
                    segment.color = "grey50") +
    scale_colour_manual(values = WORKFLOW_COLOURS) +
    labs(x = "Log2 FC per 1 umol/hr/kg M-value", y = "(-log10 P.value)",
         colour = "Workflow") +
    mvalue_theme

mvalue_ranked_pct_plot <- ggplot(
    best_results_Mvalue,
    aes(x = Rank, y = logFC, colour = -log10(P.Value),
        label = paste0(Gene_ID, " (", round(PercentageChange, 1), "%)"))
) +
    geom_point(size = 4, alpha = 0.5) +
    geom_text_repel(size = 3, box.padding = 0.2, point.padding = 0.2,
                    segment.color = "grey50") +
    scale_colour_gradient(low = "grey", high = "darkred", name = "-log10(P.Value)") +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    labs(x = "Ranked Proteins", y = "Log2 FC per 1 unit M-value") +
    mvalue_theme

print(mvalue_volcano)


