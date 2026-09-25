library(ggplot2)
library(dplyr)
library(stringr)
library(ppcor)

## Inputs
Protein_correlations <- read.delim("data-raw/UKBB_protein_Mvalue_correlations.csv", sep = ",")

ml_input           <- readRDS("data/ML_Mvalue_prediction_data/Mvalue_prediction_input.rds")
combined_clinical  <- ml_input$combined_clinical
proteome_corrected <- ml_input$proteome_corrected


partial_results_sex <- data.frame(Protein = rownames(proteome_corrected),
                                  rho = NA, pval = NA, adj_pval = NA)

for (i in seq_len(nrow(proteome_corrected))) {

  protein_values <- as.numeric(proteome_corrected[i, ])
  m_values       <- combined_clinical$M_value_z_by_cohort
  gender         <- factor(combined_clinical$Gender)

  complete_idx <- complete.cases(protein_values, m_values, gender)

  if (sum(complete_idx) >= 4) {   # require at least 4 data points
    test <- pcor.test(
      x = protein_values[complete_idx],
      y = m_values[complete_idx],
      z = data.frame(Gender = as.numeric(gender[complete_idx])),
      method = "spearman"
    )

    partial_results_sex$rho[i]  <- test$estimate
    partial_results_sex$pval[i] <- test$p.value
  }
}

partial_results_sex$adj_pval <- p.adjust(partial_results_sex$pval, method = "BH")
partial_results_sex$Gene     <- word(partial_results_sex$Protein, 1, sep = "_")



#### 2. Join to the UK Biobank correlations ####


# Match each clamp-cohort protein to its UK Biobank assay
Protein_UKBB_correlations <- Protein_correlations[
  match(partial_results_sex$Gene, Protein_correlations$assay_label),
]

# cbind keeps the clamp-cohort rho as `rho` and renames the UK Biobank rho to
# `rho.1`, because both source columns are called rho
UKBB_KI_overlap <- data.frame(partial_results_sex, Protein_UKBB_correlations)

# Keep only the proteins that matched a UK Biobank assay (one row per gene)
Protein_UKBB_correlations <- UKBB_KI_overlap[
  na.omit(match(partial_results_sex$Gene, UKBB_KI_overlap$assay_label)),
]

nrow(Protein_UKBB_correlations)          # number of plotted proteins
mean(Protein_UKBB_correlations$n)        # UK Biobank n behind the axis label

common_breaks <- seq(-1, 1, by = 0.5)   # same ticks on both axes

# Correlation between the two cohorts' rho values (the published annotation)
rho_test <- cor.test(Protein_UKBB_correlations$rho,
                     Protein_UKBB_correlations$rho.1, method = "spearman")
rho_test

ggplot(Protein_UKBB_correlations, aes(x = rho, y = rho.1)) +
    geom_point(alpha = 0.6, size = 3, color = "darkblue", fill = "#e8eefa", shape = 21) +
    geom_smooth(method = "lm", se = TRUE, color = "#4564a3", fill = "#4564a3",
                alpha = 0.1, linetype = "dashed") +
    # ADDED so the panel comes out of R complete - value from rho_test above
    annotate("text", x = -0.15, y = 0.6,
             label = paste0("ρ = ", sprintf("%.2f", unname(rho_test$estimate))),
             size = 4.5) +
    labs(
        title = "M-value plasma proteome association",
        x = "Spearman ρ (UK Biobank, n = 20389 individuals)",
        y = "Spearman ρ (Clamp cohorts, n = 161 individuals)"
    ) +
    scale_x_continuous(limits = c(-0.7, 0.7), breaks = common_breaks, expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.7, 0.7), breaks = common_breaks, expand = c(0, 0)) +
    coord_fixed(ratio = 1) +
    theme(
        text             = element_text(color = "black"),
        axis.title       = element_text(color = "black", size = 14),
        axis.text        = element_text(color = "black", size = 12),
        plot.title       = element_text(color = "black", size = 14, hjust = 0.5),
        axis.line        = element_line(color = "black"),
        axis.ticks       = element_line(color = "black"),
        panel.grid       = element_blank(),
        panel.background = element_blank()
    )
