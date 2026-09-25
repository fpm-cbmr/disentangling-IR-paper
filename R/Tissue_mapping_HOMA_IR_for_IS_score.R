# ==============================================================================
# HOMA_IR - gene/tissue mapping and preprocessing for Tissue_IS_score_prediction.R
# ==============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(tibble)

IR_INDEX <- "HOMA_IR"

dir.create("results", showWarnings = FALSE)
OUT <- function(x) file.path("results", paste0(IR_INDEX, "_", x))


# ---- 1. Load and merge HPA annotation ---------------------------------------
# The subset comes from the partial IRindex correlation

test_hpa <- readRDS("data/IRindices_HOMA_IR_subset.rds")

hpa_data <- readr::read_tsv("data-raw/hpa_24.tsv", show_col_types = FALSE)

# Separate multiple Uniprot IDs into individual rows
mult_uniprot <- hpa_data %>%
  filter(stringr::str_detect(Uniprot, ", ")) %>%
  tidyr::separate_rows(Uniprot, sep = ", ")

hpa_data <- dplyr::bind_rows(hpa_data %>% filter(!stringr::str_detect(Uniprot, ", ")),
                             mult_uniprot)

test_hpa <- test_hpa %>%
  dplyr::rename("Uniprot" = "UniProt_ID")

colnames(hpa_data) <- gsub(" ", ".", colnames(hpa_data))

merged_data <- test_hpa %>%
  dplyr::left_join(hpa_data, by = c("Uniprot" = "Uniprot"))

stopifnot(all(c("Gene", "GeneID", "Significant", "Direction",
                "RNA.tissue.specific.nTPM") %in% names(merged_data)))


# ---- 2. Protein groups -------------------------------------------------------
# Groups are pulled from `Gene` (hpa_24.tsv) and membership is matched on
# `GeneID` (the partial-correlation table). Both are gene symbols but from
# different sources, so any that fail to match are reported below - a non-zero
# count silently shrinks the associated set and shifts every odds ratio.

group_lists <- list(
  m_positive = merged_data %>%
    dplyr::filter((Significant == "YES" & Direction == "Up")) %>%
    dplyr::pull(Gene),

  m_negative = merged_data %>%
    dplyr::filter((Significant == "YES" & Direction == "Down")) %>%
    dplyr::pull(Gene),

  m_all = merged_data %>%
    dplyr::filter((Significant == "YES")) %>%
    dplyr::pull(Gene)
)

new_lists <- group_lists

message("Group sizes (", IR_INDEX, "):")
print(vapply(new_lists, length, integer(1)))
for (g in names(new_lists)) {
  n_unmatched <- length(setdiff(new_lists[[g]], merged_data$GeneID))
  if (n_unmatched > 0) {
    message("  ", g, ": ", n_unmatched, " of ", length(new_lists[[g]]),
            " Gene symbols absent from GeneID")
  }
}


# ---- 3. Tissue enrichment (Fisher's exact) -----------------------------------
# Needed because sig_tissue_proteins_long is restricted to enriched tissues.

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

  write_csv(tissue_enrichment, OUT(paste0("enrichment_tissue_", group_name, ".csv")))

  return(tissue_enrichment)
})

combined_tissue_results <- do.call(rbind, tissue_results_list)


# ---- 4. Tissue mapping -------------------------------------------------------

membership <- bind_rows(
  tibble(GeneID = new_lists$m_positive, group = "m_positive"),
  tibble(GeneID = new_lists$m_negative, group = "m_negative"),
  tibble(GeneID = new_lists$m_all,      group = "m_all")
) %>% distinct()

# One row per GeneID-tissue
tissue_map <- merged_data %>%
  dplyr::select(GeneID, Uniprot, RNA.tissue.specific.nTPM) %>%
  filter(!is.na(RNA.tissue.specific.nTPM)) %>%
  separate_rows(RNA.tissue.specific.nTPM, sep = ";") %>%
  separate(RNA.tissue.specific.nTPM, into = c("tissue", "nTPM"), sep = ": ") %>%
  mutate(
    tissue = str_trim(tolower(tissue)),
    nTPM   = suppressWarnings(as.numeric(nTPM))
  ) %>%
  distinct(GeneID, tissue, .keep_all = TRUE)

padj_cut <- 0.05
min_or   <- 1

sig_tissues <- combined_tissue_results %>%
  filter(!is.na(padj), padj < padj_cut, fisher_OR > min_or) %>%
  mutate(tissue = str_trim(tolower(tissue))) %>%
  dplyr::select(group, tissue, fisher_OR, padj) %>%
  distinct()

sig_tissue_proteins_long <- sig_tissues %>%
  inner_join(tissue_map,  by = "tissue") %>%
  inner_join(membership,  by = c("group", "GeneID")) %>%
  dplyr::select(group, tissue, GeneID, Uniprot, fisher_OR, padj) %>%
  distinct() %>%
  arrange(group, tissue, GeneID)

sig_tissue_counts <- sig_tissue_proteins_long %>%
  count(group, tissue, name = "n_proteins") %>%
  arrange(group, desc(n_proteins))

write_csv(sig_tissue_proteins_long, OUT("significant_tissue_proteins.csv"))
write_csv(sig_tissue_counts,        OUT("significant_tissue_protein_counts.csv"))

message("Significantly enriched tissues: ", length(unique(sig_tissues$tissue)))


# ---- 5. expr_long: pre-clamp expression against the IR index ----------------

Discovery_proteome_data <- readRDS("data/processed_datasets.rds")
Discovery_clinical_data <- readRDS("data/Clinical_data_filt")
Discovery_sample_data   <- readRDS("data/Sample_ID_filt.rds")

## HOMA1-IR is already present in the clinical table

ms_mats <- list(
  Neat     = Discovery_proteome_data$Neat_discovery_cohort$corrected,
  PCA      = Discovery_proteome_data$PCA_discovery_cohort$corrected,
  MagNet   = Discovery_proteome_data$MagNet_discovery_cohort$corrected,
  Depleted = Discovery_proteome_data$Depleted_discovery_cohort$corrected
)
olink_df <- Discovery_proteome_data$Olink_discovery_cohort

ms_cols     <- unique(unlist(lapply(ms_mats, colnames)))
olink_cols  <- setdiff(colnames(olink_df), "ID")
all_samples <- unique(c(ms_cols, olink_cols))
all_samples <- all_samples[!is.na(all_samples)]

# Which metadata column holds the IDs used as matrix colnames
candidate_cols <- names(Discovery_sample_data)[sapply(Discovery_sample_data, function(x) is.character(x) || is.factor(x))]
overlap_sizes  <- sapply(candidate_cols, function(colname) {
  sum(as.character(Discovery_sample_data[[colname]]) %in% all_samples, na.rm = TRUE)
})
stopifnot(any(overlap_sizes > 0))
id_col <- names(overlap_sizes)[which.max(overlap_sizes)]
message("Using metadata ID column: ", id_col)

meta_df <- Discovery_sample_data %>%
  mutate(SampleID = as.character(.data[[id_col]]))

stopifnot(nrow(Discovery_clinical_data) == nrow(Discovery_sample_data))

meta_with_M <- meta_df %>%
  bind_cols(Discovery_clinical_data) %>%
  filter(Clamp == "Pre") %>%
  transmute(
    SampleID = SampleID,
    M_value  = as.numeric(HOMA1.IR)
  ) %>%
  filter(!is.na(SampleID), !is.na(M_value))

keep_cols <- intersect(all_samples, meta_with_M$SampleID)
stopifnot(length(keep_cols) > 0)

subset_ms    <- function(mat) mat[, intersect(colnames(mat), keep_cols), drop = FALSE]
subset_olink <- function(df)  df[, c("ID", intersect(colnames(df), keep_cols)), drop = FALSE]

Neat_mat_pre     <- subset_ms(ms_mats$Neat)
PCA_mat_pre      <- subset_ms(ms_mats$PCA)
MagNet_mat_pre   <- subset_ms(ms_mats$MagNet)
Depleted_mat_pre <- subset_ms(ms_mats$Depleted)
Olink_df_pre     <- subset_olink(olink_df)

M_values_df <- meta_with_M

parse_ms_id <- function(x) {
  tibble(
    id      = x,
    Uniprot = sub("_.*$", "", x),
    Gene    = sub("^[^_]+_", "", x)
  )
}
parse_olink_id <- function(x) {
  parts <- stringr::str_split_fixed(x, "_", 3)
  tibble(
    ID      = x,
    Gene    = parts[, 1],
    Uniprot = parts[, 3]
  )
}
long_from_ms <- function(mat, wf_name) {
  ids <- parse_ms_id(rownames(mat))
  as_tibble(mat, rownames = "id") %>%
    inner_join(ids, by = "id") %>%
    pivot_longer(-c(id, Uniprot, Gene), names_to = "SampleID", values_to = "Intensity") %>%
    mutate(Workflow = wf_name)
}
long_from_olink <- function(df) {
  idmap <- parse_olink_id(df$ID)
  df %>%
    inner_join(idmap, by = dplyr::join_by(ID)) %>%
    pivot_longer(-c(ID, Uniprot, Gene), names_to = "SampleID", values_to = "Intensity") %>%
    mutate(Workflow = "Olink_discovery_cohort")
}

expr_long <- bind_rows(
  long_from_ms(Neat_mat_pre,     "Neat_discovery_cohort"),
  long_from_ms(PCA_mat_pre,      "PCA_discovery_cohort"),
  long_from_ms(MagNet_mat_pre,   "MagNet_discovery_cohort"),
  long_from_ms(Depleted_mat_pre, "Depleted_discovery_cohort"),
  long_from_olink(Olink_df_pre)
) %>%
  mutate(
    Uniprot   = str_to_upper(str_trim(Uniprot)),
    SampleID  = as.character(SampleID),
    Intensity = suppressWarnings(as.numeric(Intensity))
  ) %>%
  inner_join(M_values_df, by = "SampleID") %>%
  filter(!is.na(Intensity), !is.na(M_value))

expr_long$M_value <- log2(expr_long$M_value)

attr(expr_long, "ir_index") <- IR_INDEX

message("\nexpr_long ready for ", IR_INDEX, ": ", nrow(expr_long), " rows, ",
        dplyr::n_distinct(expr_long$SampleID), " Pre samples.")
message("M_value column holds log2(", IR_INDEX, ") - not the clamp M-value.")
message("Now run Tissue_IS_score_prediction.R.")
