###############################################
## Tissue importance analysis (LASSO, Pre only)
###############################################

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(glmnet)

## ------------------------------------------------
## 0) INPUTS (must already exist in your session)
## ------------------------------------------------
## - new_lists$m_positive / m_negative / m_all
## - tissue_map: GeneID, Uniprot, RNA.tissue.specific.nTPM (HPA parsed)
## - sig_tissue_proteins_long: group, tissue, GeneID, Uniprot, fisher_OR, padj
##   (if needed: sig_tissue_proteins_long <- read.csv("significant_tissue_proteins.csv"))
## - expr_long: SampleID, Uniprot, Gene, Workflow, Intensity, M_value (Pre only)

## If expr_long not yet cleaned:
expr_long_clean <- expr_long %>%
  mutate(
    Uniprot  = toupper(str_trim(Uniprot)),
    SampleID = as.character(SampleID)
  ) %>%
  filter(!is.na(Intensity), !is.na(M_value))

## Membership of insulin-sensitivity protein sets
membership <- bind_rows(
  tibble(GeneID = new_lists$m_positive, group = "m_positive"),
  tibble(GeneID = new_lists$m_negative, group = "m_negative"),
  tibble(GeneID = new_lists$m_all,      group = "m_all")
) %>%
  distinct()

## ------------------------------------------------
## 1) Protein–M-value correlation (for ranking features)
## ------------------------------------------------
prot_cor <- expr_long_clean %>%
  group_by(Uniprot) %>%
  summarise(
    cor_M = suppressWarnings(
      cor(Intensity, M_value, method = "spearman", use = "complete.obs")
    ),
    .groups = "drop"
  ) %>%
  filter(!is.na(cor_M))

## ------------------------------------------------
## 2) Base tissue–protein mapping from enriched tissues
##    (exclude placenta & similar)
## ------------------------------------------------

exclude_tissues <- c(
  "placenta", "testis", "vagina", "prostate",
  "seminal vesicle", "ovary", "endometrium 1",
  "epididymis", "cervix"
)

sig_map_base <- sig_tissue_proteins_long %>%
  mutate(
    tissue  = str_trim(tolower(tissue)),
    Uniprot = toupper(str_trim(Uniprot))
  ) %>%
  filter(!tissue %in% exclude_tissues) %>%
  distinct(tissue, Uniprot)

## ------------------------------------------------
## 3) Force-include skeletal muscle
##    - take proteins mapped to "skeletal muscle" in tissue_map
##    - keep only those in your insulin sensitivity sets (membership)
## ------------------------------------------------

skeletal_muscle_map <- tissue_map %>%
  mutate(
    tissue  = str_trim(tolower(tissue)),
    Uniprot = toupper(str_trim(Uniprot))
  ) %>%
  filter(tissue == "skeletal muscle") %>%
  inner_join(membership,
             by = "GeneID",
             relationship = "many-to-many") %>%  # ok: gene can be in multiple groups
  transmute(tissue, Uniprot) %>%
  distinct()

## ------------------------------------------------
## 4) Final tissue_protein_map & tissue list
## ------------------------------------------------

tissue_protein_map <- bind_rows(
  sig_map_base,
  skeletal_muscle_map
) %>%
  distinct(tissue, Uniprot)

all_tissues <- sort(unique(tissue_protein_map$tissue))

## ------------------------------------------------
## 5) Helper: build per-tissue wide matrix
##    - uses all mapped proteins for that tissue (capped by max_features)
##    - only Pre samples (already ensured via expr_long_clean)
## ------------------------------------------------

make_tissue_df <- function(tissue_name,
                           min_samples = 40,
                           max_features = 300,
                           protein_cov_thresh = 0.6) {

  # 1) Proteins mapped to this tissue (from tissue_protein_map),
  #    ranked by |cor_M|
  tp <- tissue_protein_map %>%
    filter(tissue == tissue_name) %>%
    inner_join(prot_cor, by = "Uniprot") %>%
    arrange(desc(abs(cor_M))) %>%
    slice_head(n = max_features)

  if (nrow(tp) < 2) return(NULL)
  selected_uniprot <- unique(tp$Uniprot)

  # 2) All samples with M_value (Pre only, already in expr_long_clean)
  samples <- expr_long_clean %>%
    distinct(SampleID, M_value)

  # 3) Wide matrix: one row per SampleID, cols = Uniprot
  wide_raw <- expr_long_clean %>%
    filter(Uniprot %in% selected_uniprot) %>%
    group_by(SampleID, Uniprot) %>%
    summarise(Intensity = mean(Intensity, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Uniprot, values_from = Intensity)

  # 4) Ensure all samples with M are present (even if some proteins missing)
  wide_df <- samples %>%
    left_join(wide_raw, by = "SampleID")

  # 5) Filter proteins with too low coverage across samples
  protein_cols <- setdiff(colnames(wide_df), c("SampleID", "M_value"))
  if (length(protein_cols) < 2) return(NULL)

  cov_ok <- sapply(wide_df[ , protein_cols, drop = FALSE], function(x) {
    mean(!is.na(x))
  })
  keep_proteins <- names(cov_ok)[cov_ok >= protein_cov_thresh]

  if (length(keep_proteins) < 2) return(NULL)

  wide_df <- wide_df %>%
    dplyr::select(SampleID, M_value, all_of(keep_proteins))

  # 6) Impute remaining NAs per protein (median imputation)
  for (p in keep_proteins) {
    vals <- wide_df[[p]]
    if (anyNA(vals)) {
      med <- median(vals, na.rm = TRUE)
      # if all NA, skip (shouldn't happen after cov filter)
      if (!is.na(med)) {
        vals[is.na(vals)] <- med
        wide_df[[p]] <- vals
      }
    }
  }

  # 7) Check sample size
  if (nrow(wide_df) < min_samples) return(NULL)

  list(
    tissue = tissue_name,
    df     = wide_df
  )
}


## ------------------------------------------------
## 6) Helper: CV R² from cv.glmnet (LASSO / elastic net)
##    - Uses cv.glmnet$cvm (MSE) → R² = 1 - MSE_min / var(y)
##    - No predict() dispatch issues
## ------------------------------------------------

cv_r2_lasso <- function(df,
                        alpha = 1,      # 1 = LASSO, 0.5 = elastic net
                        nfolds = 10) {

  pred_cols <- setdiff(colnames(df), c("SampleID", "M_value"))
  if (length(pred_cols) < 1) return(NA_real_)

  x <- as.matrix(df[, pred_cols, drop = FALSE])
  x <- scale(x)  # standardize predictors
  y <- df$M_value

  if (length(unique(y)) < 3) return(NA_real_)

  set.seed(123)
  fit <- try(
    glmnet::cv.glmnet(
      x, y,
      alpha  = alpha,
      nfolds = nfolds,
      family = "gaussian"
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) return(NA_real_)

  mse_min <- min(fit$cvm, na.rm = TRUE)
  var_y   <- stats::var(y, na.rm = TRUE)
  if (is.na(var_y) || var_y == 0) return(NA_real_)

  r2 <- 1 - mse_min / var_y
  as.numeric(r2)
}

## ------------------------------------------------
## 7) Run models for all tissues
## ------------------------------------------------

tissue_models <- purrr::map(all_tissues, make_tissue_df)

tissue_scores_lasso <- purrr::map_dfr(tissue_models, function(x) {
  if (is.null(x)) return(NULL)
  df <- x$df

  tibble(
    tissue     = x$tissue,
    n_proteins = ncol(df) - 2,  # minus SampleID + M_value
    n_samples  = nrow(df),
    cv_R2      = cv_r2_lasso(df)
  )
}) %>%
  filter(!is.na(cv_R2)) %>%
  mutate(
    cv_R2 = ifelse(cv_R2 < 0, 0, cv_R2)  # negative = no predictive value
  ) %>%
  arrange(desc(cv_R2))

## Inspect table
print(tissue_scores_lasso)

## ------------------------------------------------
## 8) Plot: LASSO-based tissue importance
## ------------------------------------------------

ggplot(tissue_scores_lasso,
       aes(y = reorder(tissue, cv_R2), x = cv_R2, fill = cv_R2)) +
  geom_col(alpha = 0.46, color = "black", width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", cv_R2)),
            hjust = -0.2,   # push labels a bit to the right of the bar
            size = 4,
            color = "black") +
  scale_fill_gradient(low = "grey70", high = "darkred") +
  labs(
    x = expression(paste("Cross-validated ", R^2)),
    y = NULL,
    title = "Tissue insulin sensitivity importance score"
  ) +
  coord_cartesian(xlim = c(0, max(tissue_scores_lasso$cv_R2) * 1.15)) +  # add space for labels
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),                 # remove gridlines
    axis.text.x = element_text(color = "black"),  # x-axis labels horizontal
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_text(color = "black", face = "bold"),
    axis.title.y = element_text(color = "black", face = "bold"),
    plot.title = element_text(color = "black", face = "bold", hjust = 0.5),
    axis.line = element_line(color = "black", linewidth = 0.4),
    legend.position = "none"
  )


