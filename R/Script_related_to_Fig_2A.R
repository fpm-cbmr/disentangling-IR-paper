#### Clinical correlation downstream analysis ####
# Libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

#datasets_clinical_correlation <- readRDS("data/datasets_clinical_correlation.rds")
datasets_clinical_correlation <- readRDS("data/datasets_partial_clinical_correlation.rds")

processed_datasets <- readRDS("data/processed_datasets.rds")
Clinical_data_filt <- readRDS("data/Clinical_data_filt")

MagNet_complete_cor_frame <- datasets_clinical_correlation$MagNet_baseline
PCA_complete_cor_frame <- datasets_clinical_correlation$PCA_baseline
Neat_complete_cor_frame <- datasets_clinical_correlation$Neat_baseline
Depleted_complete_cor_frame <- datasets_clinical_correlation$Depleted_baseline
Olink_complete_cor_frame <- datasets_clinical_correlation$Olink_baseline

MagNet_complete_cor_frame_post <- datasets_clinical_correlation$MagNet_insulin
PCA_complete_cor_frame_post <- datasets_clinical_correlation$PCA_insulin
Neat_complete_cor_frame_post <- datasets_clinical_correlation$Neat_insulin
Depleted_complete_cor_frame_post <- datasets_clinical_correlation$Depleted_insulin
Olink_complete_cor_frame_post <- datasets_clinical_correlation$Olink_insulin

clin_var <- c("M..µmol..kg.min..", "HbA1c..mmol.mol.", "BMI..kg.m2.", "Age",
              "FS.Insulin..mIE.L.", "HOMA1.IR", "Creatinine..µmol.L.", "P.TG..mmol.L.",
              "P.Chol..mmol.L.", "FFA..mmol.L.", "FP.Glucose..mmol.L.","X120.min.Glucose.OGTT..mmol.L.",
              "W.H.ratio", "Body.fat....")


# Total protein counts for each workflow
total_proteins <- c(
  MagNet = nrow(MagNet_complete_cor_frame),
  Neat = nrow(Neat_complete_cor_frame),
  PCA = nrow(PCA_complete_cor_frame),
  Depleted = nrow(Depleted_complete_cor_frame),
  Olink = nrow(Olink_complete_cor_frame)
)

# Define a function to calculate significant protein counts for each clinical variable
count_significant_proteins <- function(data, clin_var, threshold) {
  sapply(grep("_adj.p.value$", colnames(data), value = TRUE), function(col) {
    sum(data[[col]] < threshold, na.rm = TRUE)
  })
}

# Calculate counts for each workflow
p_value_threshold <- 0.05  # Define the significance threshold
MagNet_counts <- count_significant_proteins(MagNet_complete_cor_frame, clin_var, p_value_threshold)
Neat_counts <- count_significant_proteins(Neat_complete_cor_frame, clin_var, p_value_threshold)
PCA_counts <- count_significant_proteins(PCA_complete_cor_frame, clin_var, p_value_threshold)
Depleted_counts <- count_significant_proteins(Depleted_complete_cor_frame, clin_var, p_value_threshold)
Olink_counts <- count_significant_proteins(Olink_complete_cor_frame, clin_var, p_value_threshold)

# Combine counts into a single dataframe
significant_counts <- bind_rows(
  data.frame(Workflow = "MagNet", Clinical_Variable = names(MagNet_counts), Count = MagNet_counts, stringsAsFactors = FALSE),
  data.frame(Workflow = "Neat", Clinical_Variable = names(Neat_counts), Count = Neat_counts, stringsAsFactors = FALSE),
  data.frame(Workflow = "PCA", Clinical_Variable = names(PCA_counts), Count = PCA_counts, stringsAsFactors = FALSE),
  data.frame(Workflow = "Depleted", Clinical_Variable = names(Depleted_counts), Count = Depleted_counts, stringsAsFactors = FALSE),
  data.frame(Workflow = "Olink", Clinical_Variable = names(Olink_counts), Count = Olink_counts, stringsAsFactors = FALSE)
)


# Add the total proteins for each workflow
significant_counts <- significant_counts %>%
  mutate(Total_Proteins = case_when(
    Workflow == "MagNet" ~ total_proteins["MagNet"],
    Workflow == "Neat" ~ total_proteins["Neat"],
    Workflow == "PCA" ~ total_proteins["PCA"],
    Workflow == "Depleted" ~ total_proteins["Depleted"],
    Workflow == "Olink" ~ total_proteins["Olink"]
  ))


# Calculate counts per 1000 proteins
significant_counts$Count_per_1000 <- (significant_counts$Count / significant_counts$Total_Proteins) * 1000


# Define the mapping from old names to new names
rename_mapping <- c(
  "M..µmol..kg.min.." = "M-value",
  "HbA1c..mmol.mol." = "HbA1c",
  "BMI..kg.m2." = "BMI",
  "Age" = "Age",
  "FS.Insulin..mIE.L." = "Fasting Insulin",
  "HOMA1.IR" = "HOMA1-IR",
  "Creatinine..µmol.L." = "Creatinine",
  "P.TG..mmol.L." = "TG",
  "P.Chol..mmol.L." = "Cholesterole",
  "FFA..mmol.L." = "FFA",
  "FP.Glucose..mmol.L." = "Fasting glucose",
  "X120.min.Glucose.OGTT..mmol.L." = "OGTT_120min",
  "W.H.ratio" = "W.H.ratio",
  "Body.fat...." = "Body fat"
)


significant_counts <- significant_counts %>%
  mutate(Clean_Clinical_Variable = gsub("_adj.p.value.*", "", Clinical_Variable))

significant_counts <- significant_counts %>%
  mutate(Clean_Clinical_Variable = recode(Clean_Clinical_Variable, !!!rename_mapping))

significant_counts <- significant_counts[!significant_counts$Clean_Clinical_Variable == "OGTT_120min",]

# Plot circle map with Count_per_1000
#export in 5x7 landscape

ggplot(significant_counts[
    significant_counts$Clean_Clinical_Variable != "M-value" & significant_counts$Count > 0,
], aes(x = Clean_Clinical_Variable, y = Workflow, size = Count, color = Workflow)) +
    geom_point(alpha = 1, stroke = NA) +
    scale_size(range = c(1, 15), breaks = c(50, 100, 200, 300)) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    labs(
        title = "Plasma proteome x phenotype associations",
        x = NULL, y = NULL, size = "Count"
    ) +
    scale_color_manual(values = c("darkred", "darkgreen", "gray", "orange", "darkblue")) +
    guides(color = "none")


sum(significant_counts[significant_counts$Workflow == "Olink" & !significant_counts$Clean_Clinical_Variable == "M-value",]$Count)
sum(significant_counts[significant_counts$Workflow == "Depleted" & !significant_counts$Clean_Clinical_Variable == "M-value",]$Count)
sum(significant_counts[significant_counts$Workflow == "Neat" & !significant_counts$Clean_Clinical_Variable == "M-value",]$Count)
sum(significant_counts[significant_counts$Workflow == "PCA" & !significant_counts$Clean_Clinical_Variable == "M-value",]$Count)
sum(significant_counts[significant_counts$Workflow == "MagNet" & !significant_counts$Clean_Clinical_Variable == "M-value",]$Count)


#### Count unique baseline associations ####
library(dplyr)
library(purrr)
library(stringr)

# Helper function: get (Protein ID, Trait) pairs from each dataset
get_significant_pairs <- function(data, threshold = 0.05) {
    adj_p_cols <- grep("_adj.p.value$", colnames(data), value = TRUE)
    sig_pairs <- map_dfr(adj_p_cols, function(col) {
        base_trait <- gsub("_adj.p.value$", "", col)
        data %>%
            filter(.data[[col]] < threshold) %>%
            transmute(ID, Trait = base_trait)
    })
    return(sig_pairs)
}

# Apply to all workflows
all_significant_pairs <- bind_rows(
    get_significant_pairs(MagNet_complete_cor_frame),
    get_significant_pairs(PCA_complete_cor_frame),
    get_significant_pairs(Neat_complete_cor_frame),
    get_significant_pairs(Depleted_complete_cor_frame),
    get_significant_pairs(Olink_complete_cor_frame)
)

# Remove duplicates of same (protein, trait) pairs across workflows
unique_protein_trait_pairs <- distinct(all_significant_pairs)

# Count total number of unique (protein, trait) associations
total_unique_associations <- nrow(unique_protein_trait_pairs)

# Count number of unique proteins involved in any trait association
unique_proteins <- n_distinct(unique_protein_trait_pairs$ID)

# Count number of unique traits with at least one protein
unique_traits <- n_distinct(unique_protein_trait_pairs$Trait)

# Output summary
cat("✅ Total unique (protein–trait) associations:", total_unique_associations, "\n")
cat("🔢 Unique proteins involved:", unique_proteins, "\n")
cat("📊 Traits with significant associations:", unique_traits, "\n")



#### barplot summary of M-value associations ####
library(ggplot2)
library(patchwork)  # For combining plots


MagNet_complete_cor_frame <- datasets_clinical_correlation$MagNet_baseline
PCA_complete_cor_frame <- datasets_clinical_correlation$PCA_baseline
Neat_complete_cor_frame <- datasets_clinical_correlation$Neat_baseline
Depleted_complete_cor_frame <- datasets_clinical_correlation$Depleted_baseline
Olink_complete_cor_frame <- datasets_clinical_correlation$Olink_baseline

nrow(MagNet_complete_cor_frame[MagNet_complete_cor_frame$M..µmol..kg.min.._adj.p.value < 0.05 & MagNet_complete_cor_frame$M..µmol..kg.min.._estimate > 0,])
nrow(MagNet_complete_cor_frame[MagNet_complete_cor_frame$M..µmol..kg.min.._adj.p.value < 0.05 & MagNet_complete_cor_frame$M..µmol..kg.min.._estimate < 0,])

nrow(Neat_complete_cor_frame[Neat_complete_cor_frame$M..µmol..kg.min.._adj.p.value < 0.05 & Neat_complete_cor_frame$M..µmol..kg.min.._estimate > 0,])
nrow(Neat_complete_cor_frame[Neat_complete_cor_frame$M..µmol..kg.min.._adj.p.value < 0.05 & Neat_complete_cor_frame$M..µmol..kg.min.._estimate < 0,])

nrow(Depleted_complete_cor_frame[Depleted_complete_cor_frame$M..µmol..kg.min.._adj.p.value < 0.05 & Depleted_complete_cor_frame$M..µmol..kg.min.._estimate > 0,])
nrow(Depleted_complete_cor_frame[Depleted_complete_cor_frame$M..µmol..kg.min.._adj.p.value < 0.05 & Depleted_complete_cor_frame$M..µmol..kg.min.._estimate < 0,])

nrow(PCA_complete_cor_frame[PCA_complete_cor_frame$M..µmol..kg.min.._adj.p.value < 0.05 & PCA_complete_cor_frame$M..µmol..kg.min.._estimate > 0,])
nrow(PCA_complete_cor_frame[PCA_complete_cor_frame$M..µmol..kg.min.._adj.p.value < 0.05 & PCA_complete_cor_frame$M..µmol..kg.min.._estimate < 0,])

nrow(Olink_complete_cor_frame[Olink_complete_cor_frame$M..µmol..kg.min.._adj.p.value < 0.05 & Olink_complete_cor_frame$M..µmol..kg.min.._estimate > 0,])
nrow(Olink_complete_cor_frame[Olink_complete_cor_frame$M..µmol..kg.min.._adj.p.value < 0.05 & Olink_complete_cor_frame$M..µmol..kg.min.._estimate < 0,])


m_value_counts <- significant_counts[significant_counts$Clinical_Variable == "M..µmol..kg.min.._adj.p.value",]


# Order Workflows by Count (from lowest to highest)
m_value_counts$Workflow <- factor(m_value_counts$Workflow,
                                  levels = m_value_counts$Workflow[order(m_value_counts$Count)])

# Define a custom theme
custom_theme <- theme_minimal() +
    theme(legend.position = "none",
          axis.text = element_text(size = 12, color = "black"),  # Black axis text
          axis.title = element_text(size = 14, color = "black"),  # Black axis titles
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "black"),  # Black bold title
          axis.line = element_line(color = "black"),  # Add x and y axis lines
          panel.grid.major = element_blank(),  # Remove major grid lines
          panel.grid.minor = element_blank()   # Remove minor grid lines
    )

# First bar plot for Count (ordered)
p1 <- ggplot(m_value_counts, aes(x = Workflow, y = Count, fill = Workflow)) +
    geom_bar(stat = "identity", alpha = 0.7) +
    geom_text(aes(label = Count), vjust = -0.5, color = "black", size = 5) +  # Add labels above bars
    labs(title = "Insulin sensitivity associations",
         x = "",
         y = "Number of associations") +
    custom_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c("MagNet" = "darkgreen",
                                 "Neat" = "darkgrey",
                                 "PCA" = "darkblue",
                                 "Depleted" = "darkred",
                                 "Olink" = "darkorange")) +
    guides(fill = "none")  # Remove legend for cleaner look

# Second bar plot for Count per 1000 (keep the same order)
p2 <- ggplot(m_value_counts, aes(x = Workflow, y = Count_per_1000, fill = Workflow)) +
    geom_bar(stat = "identity", alpha = 0.7) +
    geom_text(aes(label = round(Count_per_1000, 0)), vjust = -0.5, color = "black", size = 5) +  # Add labels above bars
    labs(title = "Insulin sensitivity associations",
         x = "",
         y = "Number of associations per 1000 proteins") +
    custom_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c("MagNet" = "darkgreen",
                                 "Neat" = "darkgrey",
                                 "PCA" = "darkblue",
                                 "Depleted" = "darkred",
                                 "Olink" = "darkorange"))

# Combine the two plots side by side
p1 + p2

library(ggpattern)

#### U-plot of all M-value associations ####
MagNet_Mvalue_sign <- MagNet_complete_cor_frame[MagNet_complete_cor_frame$M..µmol..kg.min.._adj.p.value < 0.05,]
MagNet_Mvalue_sign$ID <- word(MagNet_Mvalue_sign$ID,2,sep="_")
Neat_Mvalue_sign <- Neat_complete_cor_frame[Neat_complete_cor_frame$M..µmol..kg.min.._adj.p.value < 0.05,]
Neat_Mvalue_sign$ID <- word(Neat_Mvalue_sign$ID,2,sep="_")
PCA_Mvalue_sign <- PCA_complete_cor_frame[PCA_complete_cor_frame$M..µmol..kg.min.._adj.p.value < 0.05,]
PCA_Mvalue_sign$ID <- word(PCA_Mvalue_sign$ID,2,sep="_")
Depleted_Mvalue_sign <- Depleted_complete_cor_frame[Depleted_complete_cor_frame$M..µmol..kg.min.._adj.p.value < 0.05,]
Depleted_Mvalue_sign$ID <- word(Depleted_Mvalue_sign$ID,2,sep="_")

Olink_Mvalue_sign <- Olink_complete_cor_frame[Olink_complete_cor_frame$M..µmol..kg.min.._adj.p.value < 0.05,]
Olink_Mvalue_sign$ID <- word(Olink_Mvalue_sign$ID,1,sep="_")


Sign_Plasma_Mvalue_background <- unique(c(word(MagNet_Mvalue_sign$ID,2,sep="_"),word(Neat_Mvalue_sign$ID,2,sep="_"), word(PCA_Mvalue_sign$ID,2,sep="_"),word(Depleted_Mvalue_sign$ID,2,sep="_"), word(Olink_Mvalue_sign$ID,1,sep="_")))

Sign_Plasma_Mvalue <- Reduce(intersect,list(c(MagNet_Mvalue_sign$ID,Neat_Mvalue_sign$ID, PCA_Mvalue_sign$ID, Depleted_Mvalue_sign$ID, Olink_Mvalue_sign$ID)))

common_genes <- Reduce(intersect, list(Olink_proteins, MagNet_proteins, Neat_proteins, PCA_proteins, Depleted_proteins))

# Combine all significant proteins into one data frame, adding a column for workflow
combined_significant <- bind_rows(
  MagNet_Mvalue_sign %>% mutate(Workflow = "MagNet"),
  Neat_Mvalue_sign %>% mutate(Workflow = "Neat"),
  PCA_Mvalue_sign %>% mutate(Workflow = "PCA"),
  Depleted_Mvalue_sign %>% mutate(Workflow = "Depleted"),
  Olink_Mvalue_sign %>% mutate(Workflow = "Olink")
)

# Filter to include only unique proteins with the most significant p-value
unique_significant <- combined_significant %>%
  group_by(ID) %>%
  slice_min(order_by = M..µmol..kg.min.._p.value, n = 1) %>%
  ungroup()


#### Export unique significance table ####
# Add UniProt_ID column to each workflow
MagNet_Mvalue_sign <- inner_join(MagNet_complete_cor_frame, MagNet_complete_cor_frame_post, by = "ID") %>%
    filter(M..µmol..kg.min.._adj.p.value.x < 0.05 | M..µmol..kg.min.._adj.p.value.y < 0.05)
MagNet_Mvalue_sign$UniProt_ID <- word(MagNet_Mvalue_sign$ID, 1, sep = "_")
MagNet_Mvalue_sign$ID <- word(MagNet_Mvalue_sign$ID, 2, sep = "_")

Neat_Mvalue_sign <- inner_join(Neat_complete_cor_frame, Neat_complete_cor_frame_post, by = "ID") %>%
    filter(M..µmol..kg.min.._adj.p.value.x < 0.05 | M..µmol..kg.min.._adj.p.value.y < 0.05)
Neat_Mvalue_sign$UniProt_ID <- word(Neat_Mvalue_sign$ID, 1, sep = "_")
Neat_Mvalue_sign$ID <- word(Neat_Mvalue_sign$ID, 2, sep = "_")

PCA_Mvalue_sign <- inner_join(PCA_complete_cor_frame, PCA_complete_cor_frame_post, by = "ID") %>%
    filter(M..µmol..kg.min.._adj.p.value.x < 0.05 | M..µmol..kg.min.._adj.p.value.y < 0.05)
PCA_Mvalue_sign$UniProt_ID <- word(PCA_Mvalue_sign$ID, 1, sep = "_")
PCA_Mvalue_sign$ID <- word(PCA_Mvalue_sign$ID, 2, sep = "_")

Depleted_Mvalue_sign <- inner_join(Depleted_complete_cor_frame, Depleted_complete_cor_frame_post, by = "ID") %>%
    filter(M..µmol..kg.min.._adj.p.value.x < 0.05 | M..µmol..kg.min.._adj.p.value.y < 0.05)
Depleted_Mvalue_sign$UniProt_ID <- word(Depleted_Mvalue_sign$ID, 1, sep = "_")
Depleted_Mvalue_sign$ID <- word(Depleted_Mvalue_sign$ID, 2, sep = "_")

Olink_Mvalue_sign <- inner_join(Olink_complete_cor_frame, Olink_complete_cor_frame_post, by = "ID") %>%
    filter(M..µmol..kg.min.._adj.p.value.x < 0.05 | M..µmol..kg.min.._adj.p.value.y < 0.05)
Olink_Mvalue_sign$UniProt_ID <- word(Olink_Mvalue_sign$ID, 3, sep = "_")
Olink_Mvalue_sign$ID <- word(Olink_Mvalue_sign$ID, 1, sep = "_")

# Combine all significant proteins into one data frame, adding Workflow column
combined_significant <- bind_rows(
    MagNet_Mvalue_sign %>% mutate(Workflow = "MagNet"),
    Neat_Mvalue_sign %>% mutate(Workflow = "Neat"),
    PCA_Mvalue_sign %>% mutate(Workflow = "PCA"),
    Depleted_Mvalue_sign %>% mutate(Workflow = "Depleted"),
    Olink_Mvalue_sign %>% mutate(Workflow = "Olink")
)

# Reorganize columns to ensure ID and UniProt_ID are together
combined_significant <- combined_significant %>%
    dplyr::select(ID, UniProt_ID, everything())


# Filter unique proteins with the most significant p-value across all workflows
All_unique_significant <- combined_significant %>%
    group_by(ID) %>%
    slice_min(order_by = pmin(M..µmol..kg.min.._p.value.x, M..µmol..kg.min.._p.value.y), n = 1) %>%
    ungroup()

# Add a column to indicate significance
All_unique_significant <- All_unique_significant %>%
    mutate(
        Significant = case_when(
            M..µmol..kg.min.._adj.p.value.x < 0.05 & M..µmol..kg.min.._adj.p.value.y < 0.05 ~ "Both",
            M..µmol..kg.min.._adj.p.value.x < 0.05 ~ "Pre_only",
            M..µmol..kg.min.._adj.p.value.y < 0.05 ~ "Post_only"
        )
    )
# Get the original column names
col_names <- colnames(All_unique_significant)

# Define fixed first four columns
first_cols <- c("ID", "UniProt_ID", "Significant", "Workflow")

# Define M-value columns (placed immediately after the first four)
m_value_cols <- col_names[c(3:5, 45:47)]  # Select manually to maintain order

# Extract remaining clinical variable indices (excluding first four and M-value)
clinical_indices_x <- seq(6, 43, by = 3)  # Start at 6, step by 3 for .x values
clinical_indices_y <- seq(48, 85, by = 3) # Start at 48, step by 3 for .y values

# Reorder clinical variables in p.value.x → adj.p.value.x → estimate.x → p.value.y → adj.p.value.y → estimate.y order
ordered_clinical_vars <- as.vector(rbind(
    clinical_indices_x, clinical_indices_x + 1, clinical_indices_x + 2,  # .x triplicates
    clinical_indices_y, clinical_indices_y + 1, clinical_indices_y + 2   # .y triplicates
))

# Generate final ordered column names
ordered_cols <- c(first_cols, m_value_cols, col_names[ordered_clinical_vars])

# Reorder the dataframe
All_unique_significant <- All_unique_significant[, ordered_cols]

# Check first few columns
head(All_unique_significant)


# Rename columns for clarity
colnames(All_unique_significant) <- gsub("\\.x$", "_pre", colnames(All_unique_significant))
colnames(All_unique_significant) <- gsub("\\.y$", "_post", colnames(All_unique_significant))


# Define the rename mapping
rename_mapping <- c(
    "M..µmol..kg.min.." = "M-value",
    "HbA1c..mmol.mol." = "HbA1c",
    "BMI..kg.m2." = "BMI",
    "Age" = "Age",
    "FS.Insulin..mIE.L." = "Fasting Insulin",
    "HOMA1.IR" = "HOMA1-IR",
    "Creatinine..µmol.L." = "Creatinine",
    "P.TG..mmol.L." = "TG",
    "P.Chol..mmol.L." = "Cholesterol",
    "FFA..mmol.L." = "FFA",
    "FP.Glucose..mmol.L." = "Fasting glucose",
    "X120.min.Glucose.OGTT..mmol.L." = "OGTT_120min",
    "W.H.ratio" = "W.H.ratio",
    "Body.fat...." = "Body fat"
)

# Function to remap column names
remap_column_names <- function(column_names, rename_mapping) {
    sapply(column_names, function(col) {
        # Extract the base variable name (before suffixes like "_pre" or "_post")
        base_name <- sub("_.*", "", col)
        # Check if the base name is in the rename_mapping
        new_base_name <- ifelse(base_name %in% names(rename_mapping), rename_mapping[base_name], base_name)
        # Reassemble the column name with the suffix
        suffix <- sub("^[^_]+", "", col)
        paste0(new_base_name, suffix)
    })
}

# Apply the remapping to the merged_data column names
colnames(All_unique_significant) <- remap_column_names(colnames(All_unique_significant), rename_mapping)

output_file <- "data-raw/All_unique_significant_results.txt"
# Write the data frame to the text file
write.table(
    All_unique_significant,
    file = output_file,
    sep = "\t",            # Use tab as the delimiter
    row.names = FALSE,     # Do not write row names
    col.names = TRUE,      # Include column headers
    quote = FALSE          # Do not include quotes around values
)



length(which(All_unique_significant$Significant == "Post_only"))
length(which(All_unique_significant$Significant == "Pre_only"))
length(which(All_unique_significant$Significant == "Both"))

#### export for disease/cell type specific enrichment ####

# Extract GeneID and ProteinID/UniProtID for each workflow
Olink_geneID <- word(Olink_complete_cor_frame$ID, 1, sep = "_")
PCA_geneID <- word(PCA_complete_cor_frame$ID, 2, sep = "_")
Depleted_geneID <- word(Depleted_complete_cor_frame$ID, 2, sep = "_")
Neat_geneID <- word(Neat_complete_cor_frame$ID, 2, sep = "_")
MagNet_geneID <- word(MagNet_complete_cor_frame$ID, 2, sep = "_")

# Combine all GeneIDs and UniProt IDs
GeneID_list <- c(Olink_geneID, PCA_geneID, Depleted_geneID, Neat_geneID, MagNet_geneID)

UniProtID_list <- c(
    word(Olink_complete_cor_frame$ID, 3, sep = "_"),
    word(PCA_complete_cor_frame$ID, 1, sep = "_"),
    word(Depleted_complete_cor_frame$ID, 1, sep = "_"),
    word(Neat_complete_cor_frame$ID, 1, sep = "_"),
    word(MagNet_complete_cor_frame$ID, 1, sep = "_")
)
All_unique_significant$`M-value_estimate_pre`
# Determine direction of change (up or down) based on M-value estimates
direction_list <- sapply(GeneID_list, function(gene) {
    if (gene %in% All_unique_significant$ID) {
        significant_row <- All_unique_significant[All_unique_significant$ID == gene, ]

        # Use M-value from "Pre" if significant in Pre-only or Both
        if (significant_row$Significant %in% c("Pre_only", "Both")) {
            mvalue <- significant_row$`M-value_estimate_pre`
        } else {
            mvalue <- significant_row$`M-value_estimate_post`
        }

        if (mvalue > 0) {
            return("Up")
        } else if (mvalue < 0) {
            return("Down")
        } else {
            return("No Change")
        }
    }
    return(NA)  # If GeneID is not found in All_unique_significant
})

# Combine the data into a data frame
protein_status_df <- data.frame(
    GeneID = GeneID_list,
    UniProt_ID = UniProtID_list,
    Is_Significant = ifelse(GeneID_list %in% All_unique_significant$ID, "Yes", "No"),
    Direction = direction_list
)

# Filter the dataframe for unique GeneIDs only
protein_status_df_unique <- protein_status_df %>%
    distinct(GeneID, .keep_all = TRUE)  # Keep only unique GeneIDs

# Print the final dataframe
print(head(protein_status_df_unique))


enrichment_file <- "data-raw/M_value_enrichment.txt"
# Write the data frame to the text file
write.table(
    protein_status_df_unique,
    file = enrichment_file,
    sep = "\t",            # Use tab as the delimiter
    row.names = FALSE,     # Do not write row names
    col.names = TRUE,      # Include column headers
    quote = FALSE          # Do not include quotes around values
)



library(ggbreak)
#export in 7.2x6cm landscape
ggplot(All_unique_significant,
       aes(x = `M-value_estimate_pre`,
           y = `M-value_estimate_post`,
           color = Workflow,
           shape = Significant)) +   # <-- add this

    geom_point(size = 2.5, alpha = 0.4, stroke = 0.8) +  # stroke works only for filled shapes

    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.8) +

    scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-0.8, 0.8),
                       breaks = seq(-0.8, 0.8, by = 0.2)) +

    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-0.8, 0.8),
                       breaks = seq(-0.8, 0.8, by = 0.2)) +

    scale_shape_manual(values = c(
        "Both" = 16,       # filled circle
        "Pre_only" = 17,   # triangle
        "Post_only" = 15   # square
    )) +

    scale_color_manual(values = c("darkred","darkgreen","darkgrey","orange","darkblue")) +

    labs(title = "Insulin sensitivity plasma proteome association",
         x = "Fasted (rho estimate)",
         y = "Insulin-stimulated (rho estimate)") +

    theme_bw() +
    theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        axis.text.x.top = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 14)
    )


#### ORA of IS-associated proteins ####







#### Overrepresentation of unique proteins for each workflow ####
Plasma_protein_background <- unique(c(Olink_proteins,MagNet_proteins, Neat_proteins, PCA_proteins,Depleted_proteins))

Olink_unique

## ORA
library(clusterProfiler)
library(org.Hs.eg.db)
go_enrich <- simplify(enrichGO(gene = All_unique_significant[All_unique_significant$`M-value_estimate_pre` < 0,]$ID,
                               universe = Plasma_gene_background,
                               OrgDb = org.Hs.eg.db,
                               keyType = 'SYMBOL',
                               readable = T,
                               ont = "CC",
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.05))

go_MF_Negative <- go_enrich@result
go_BP_Negative <- go_enrich@result
go_CC_Negative <- go_enrich@result

go_enrich <- simplify(enrichGO(gene = All_unique_significant[All_unique_significant$`M-value_estimate_pre` > 0,]$ID,
                               universe = Plasma_gene_background,
                               OrgDb = org.Hs.eg.db,
                               keyType = 'SYMBOL',
                               readable = T,
                               ont = "CC",
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.05))

go_MF_Positive <- go_enrich@result
go_BP_Positive <- go_enrich@result
go_CC_Positive <- go_enrich@result










#### Reactome pathway overrepresentation ####
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(dplyr)
library(forcats)
library(tidyr)

# **1️⃣ Convert Gene Symbols to ENTREZ IDs**
# Reactome requires ENTREZ IDs, so we must convert SYMBOLs → ENTREZIDs
genes_reactome_positive <- bitr(
    All_unique_significant[All_unique_significant$`M-value_estimate_pre` > 0,]$ID,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
)

genes_reactome_negative <- bitr(
    All_unique_significant[All_unique_significant$`M-value_estimate_pre` < 0,]$ID,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
)

# Convert the background genes to ENTREZ IDs
background_reactome <- bitr(
    Plasma_gene_background,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
)

# **2️⃣ Perform Reactome Pathway Over-Representation Analysis (ORA)**
reactome_enrich_positive <- enrichPathway(
    gene = genes_reactome_positive$ENTREZID,
    universe = background_reactome$ENTREZID,
    organism = "human",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE  # Convert ENTREZ IDs to gene symbols in results
)

# **2️⃣ Perform Reactome Pathway Over-Representation Analysis (ORA)**
reactome_enrich_negative <- enrichPathway(
    gene = genes_reactome_negative$ENTREZID,
    universe = background_reactome$ENTREZID,
    organism = "human",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE  # Convert ENTREZ IDs to gene symbols in results
)

# **3️⃣ Visualize Reactome Pathway Enrichment**
# Barplot of enriched Reactome pathways

Reactome_negative <- reactome_enrich_negative@result
Reactome_positive <- reactome_enrich_positive@result


# **1️⃣ Convert GeneRatio to Numeric**
convert_ratio <- function(ratio) {
    sapply(strsplit(ratio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
}

# Process Negative Reactome Enrichment
Reactome_negative <- reactome_enrich_negative@result %>%
    mutate(GeneRatio = convert_ratio(GeneRatio),
           BgRatio = convert_ratio(BgRatio),
           OddsRatio = (GeneRatio / (1 - GeneRatio)) / (BgRatio / (1 - BgRatio)),
           Log2_OddsRatio = log2(OddsRatio)) %>%  # Compute log2(Odds Ratio)
    arrange(p.adjust) %>%
    slice_head(n = 5)  # Select top 5 enriched terms

# Process Positive Reactome Enrichment
Reactome_positive <- reactome_enrich_positive@result %>%
    mutate(GeneRatio = convert_ratio(GeneRatio),
           BgRatio = convert_ratio(BgRatio),
           OddsRatio = (GeneRatio / (1 - GeneRatio)) / (BgRatio / (1 - BgRatio)),
           Log2_OddsRatio = log2(OddsRatio)) %>%  # Compute log2(Odds Ratio)
    arrange(p.adjust) %>%
    slice_head(n = 5)  # Select top 5 enriched terms

# Flip Log2 Odds Ratio for Negative Enrichment (to appear on left side)
Reactome_negative$Log2_OddsRatio <- -Reactome_negative$Log2_OddsRatio
Reactome_negative$Direction <- "Negative"
Reactome_positive$Direction <- "Positive"

# **2️⃣ Combine Data**
Reactome_combined <- bind_rows(Reactome_negative, Reactome_positive)

ggplot(Reactome_combined, aes(x = fct_reorder(Description, Log2_OddsRatio),
                              y = Log2_OddsRatio, fill = p.adjust)) +
    geom_bar(stat = "identity") +
    coord_flip() +  # Flip axes for better readability
    theme_minimal() +
    geom_vline(xintercept = 0, color = "black", linetype = "solid") +  # Keep only the central vertical axis
    geom_hline(yintercept = 0, color = "black", linetype = "solid") +  # X-axis at the bottom
    scale_fill_gradient(low = "grey", high = "darkred", name = "Adjusted p-value") +  # Adjusted p-value color gradient
    labs(title = "Top 5 Reactome Pathways (Positive & Negative Enrichment)",
         x = "Pathway",
         y = "Log2(Odds Ratio)") +
    theme(
        axis.line.x = element_line(color = "black"),  # Keep x-axis at bottom
        axis.line.y = element_blank(),  # Remove left y-axis
        text = element_text(color = "black"),  # Make all text black
        axis.text = element_text(color = "black"),  # Make tick labels black
        axis.title = element_text(color = "black"),  # Make axis titles black
        legend.text = element_text(color = "black"),  # Make legend text black
        legend.title = element_text(color = "black")  # Make legend title black
    )
