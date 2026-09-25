# Load required libraries
library(pheatmap)
library(VennDiagram)
library(purrr)
library(tidyverse)
library(PhosR)
library(ggplot2)
library(reshape2)

processed_datasets <- readRDS("data/processed_datasets.rds")

Olink_proteins <- unique(word(processed_datasets$Olink_discovery_cohort$ID,1,sep="_"))
MagNet_proteins <- unique(word(word(rownames(processed_datasets$MagNet_discovery_cohort$raw),2,sep="_"),1,sep=";"))
PCA_proteins <- unique(word(word(rownames(processed_datasets$PCA_discovery_cohort$raw),2,sep="_"),1,sep=";"))
Neat_proteins <- unique(word(word(rownames(processed_datasets$Neat_discovery_cohort$raw),2,sep="_"),1,sep=";"))
Depleted_proteins <- unique(word(word(rownames(processed_datasets$Depleted_discovery_cohort$raw),2,sep="_"),1,sep=";"))


#### unique proteins per workflow ####
# Proteins unique to the Olink workflow
Olink_unique <- setdiff(Olink_proteins, c(MagNet_proteins, Neat_proteins, PCA_proteins,Depleted_proteins))

# Proteins unique to the MagNet workflow
MagNet_unique <- setdiff(MagNet_proteins, c(Olink_proteins, Neat_proteins, PCA_proteins,Depleted_proteins))

# Proteins unique to the Neat workflow
Neat_unique <- setdiff(Neat_proteins, c(Olink_proteins, MagNet_proteins, PCA_proteins,Depleted_proteins))

# Proteins unique to the PCA workflow
PCA_unique <- setdiff(PCA_proteins, c(Olink_proteins, MagNet_proteins, Neat_proteins,Depleted_proteins))

# Proteins unique to the PCA workflow
Depleted_unique <- setdiff(Depleted_proteins, c(Olink_proteins, MagNet_proteins, Neat_proteins, PCA_proteins))


Plasma_protein_background <- unique(c(Olink_proteins,MagNet_proteins, Neat_proteins, PCA_proteins, Depleted_proteins))
common_genes <- Reduce(intersect, list(Olink_proteins, MagNet_proteins, Neat_proteins, PCA_proteins, Depleted_proteins))



Common_genes_MagNet <- processed_datasets$MagNet_discovery_cohort$corrected[which(word(word(rownames(processed_datasets$MagNet_discovery_cohort$corrected),2,sep="_"),1,sep=";") %in%  common_genes),]
rownames(Common_genes_MagNet) <- word(word(rownames(Common_genes_MagNet),2,sep="_"),sep=";")
Common_genes_MagNet <- Common_genes_MagNet[order(rownames(Common_genes_MagNet)),]

Common_genes_PCA <- processed_datasets$PCA_discovery_cohort$log2[which(word(word(rownames(processed_datasets$PCA_discovery_cohort$log2),2,sep="_"),1,sep=";") %in%  common_genes),]
rownames(Common_genes_PCA) <- word(word(rownames(Common_genes_PCA),2,sep="_"),sep=";")
Common_genes_PCA <- Common_genes_PCA[order(rownames(Common_genes_PCA)),]


Common_genes_Neat <- processed_datasets$Neat_discovery_cohort$log2[which(word(word(rownames(processed_datasets$Neat_discovery_cohort$log2),2,sep="_"),1,sep=";") %in%  common_genes),]
rownames(Common_genes_Neat) <- word(word(rownames(Common_genes_Neat),2,sep="_"),1,sep=";")
Common_genes_Neat <- Common_genes_Neat[order(rownames(Common_genes_Neat)),]

Common_genes_Depleted <- processed_datasets$Depleted_discovery_cohort$log2[which(word(word(rownames(processed_datasets$Depleted_discovery_cohort$log2),2,sep="_"),1,sep=";") %in%  common_genes),]
rownames(Common_genes_Depleted) <- word(word(rownames(Common_genes_Depleted),2,sep="_"),1,sep=";")
Common_genes_Depleted <- Common_genes_Depleted[order(rownames(Common_genes_Depleted)),]

Common_genes_Olink <- processed_datasets$Olink_discovery_cohort[which(word(processed_datasets$Olink_discovery_cohort$ID,1,sep="_") %in%  common_genes),]
rownames(Common_genes_Olink) <- word(Common_genes_Olink$ID,1,sep="_")
Common_genes_Olink <- as.matrix(Common_genes_Olink[,-1])

#### 155 shared complete quantified proteins ####

# List of all workflows
all_workflows <- list(Common_genes_Olink, Common_genes_PCA, Common_genes_Depleted, Common_genes_Neat, Common_genes_MagNet)

# Find proteins with any missing values across workflows
common_na_mask <- Reduce(`|`, lapply(all_workflows, function(x) rowSums(is.na(x)) > 0))

# Remove these proteins from all workflows to maintain alignment
clean_workflows <- lapply(all_workflows, function(x) x[!common_na_mask, ])
names(clean_workflows) <- c("Olink", "PCA", "Depleted", "Neat", "MagNet")

# Verify all workflows have the same dimensions
sapply(clean_workflows, dim)  # Should return (same number of proteins × 150 samples) for all



# Load required libraries
library(ggplot2)
library(pheatmap)

# Define the max CV threshold
max_cv_threshold <- 100  # Adjust as needed

# Initialize a dataframe to store CV values
cv_matrix <- data.frame(Protein = rownames(clean_workflows[[1]]))

# Loop through workflows, transform back to linear scale, and compute CV
for (workflow in names(clean_workflows)) {
    dataset <- clean_workflows[[workflow]]

    # Transform data back to linear scale (undo log2 transformation)
    dataset_linear <- 2^dataset

    # Compute CV per protein
    cv_values <- apply(dataset_linear, 1, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100)

    # Add CV values to the dataframe
    cv_matrix[[workflow]] <- cv_values
}

# Set protein names as row names and remove Protein column
rownames(cv_matrix) <- cv_matrix$Protein
cv_matrix <- cv_matrix[, -1]

# Filter out proteins exceeding the max allowed CV threshold
cv_matrix_filtered <- cv_matrix[apply(cv_matrix, 1, function(x) all(x <= max_cv_threshold, na.rm = TRUE)), ]

# Perform hierarchical clustering **once** on the full dataset
row_dendrogram <- hclust(dist(cv_matrix_filtered))

# Extract the protein order from the full heatmap
protein_order <- rownames(cv_matrix_filtered)[row_dendrogram$order]

# Generate the original heatmap **with clustering**
pheatmap(cv_matrix_filtered,
         cluster_rows = row_dendrogram,
         cluster_cols = TRUE,
         cutree_rows = 9,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         na_col = "grey",
         main = paste("Heatmap of CV Across Workflows (Max CV =", max_cv_threshold, ")"),
         border_color = NA,
         fontsize_row = 8,
         fontsize_col = 10,
         scale = "none",
         show_rownames = TRUE,
         show_colnames = TRUE,
         legend = TRUE
)









# Perform hierarchical clustering on the filtered matrix
row_dendrogram <- hclust(dist(cv_matrix_filtered))

# Extract clustering order
protein_order <- rownames(cv_matrix_filtered)[row_dendrogram$order]

# Cut tree into clusters (as before)
protein_clusters <- cutree(row_dendrogram, k = 9)

# Create a dataframe mapping proteins to clusters
cluster_df <- data.frame(
    Protein = names(protein_clusters),
    Cluster = protein_clusters,
    stringsAsFactors = FALSE
)



# Compute gaps from the full dataset
full_gaps_row <- cumsum(table(protein_clusters))

# Extract proteins from clusters 8, 6, 5, 9, and 4
selected_clusters <- c(8, 6, 5, 9, 4)
selected_proteins <- cluster_df$Protein[cluster_df$Cluster %in% selected_clusters]

# Preserve original clustering order
zoomed_protein_order <- protein_order[protein_order %in% selected_proteins]

# Reorder cv_matrix_zoomed to match the full heatmap order
cv_matrix_zoomed <- cv_matrix_filtered[zoomed_protein_order, , drop = FALSE]

# Extract cluster assignments for the zoomed-in dataset
zoomed_clusters <- protein_clusters[zoomed_protein_order]

saveRDS(zoomed_clusters,"zoomed_clusters_CV.rds")

# Recalculate `gaps_row` based on the original cluster gaps
zoomed_gaps_row <- which(diff(as.numeric(zoomed_clusters)) != 0)

# Add last cluster spacing if needed
if (!length(zoomed_gaps_row) == 0 && tail(zoomed_gaps_row, n=1) < nrow(cv_matrix_zoomed)) {
    zoomed_gaps_row <- c(zoomed_gaps_row, nrow(cv_matrix_zoomed))
}

# Generate the zoomed-in heatmap **with identical cluster spacing**
# Define the desired column order
desired_col_order <- c("PCA", "Depleted", "Neat", "Olink", "MagNet")

# Reorder columns in the zoomed dataset
cv_matrix_zoomed <- cv_matrix_zoomed[, desired_col_order]

pheatmap(cv_matrix_zoomed,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_row = zoomed_gaps_row,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         na_col = "grey",
         main = paste("Zoomed-in Heatmap: Clusters", paste(selected_clusters, collapse = ", ")),
         border_color = NA,
         fontsize_row = 8,
         fontsize_col = 10,
         scale = "none",
         show_rownames = TRUE,
         show_colnames = TRUE,
         legend = TRUE)







# Load required libraries
library(ggplot2)
library(reshape2)  # For melting data to long format
library(ggrepel)   # For better text labeling

# Define the max CV threshold (Adjust as needed)
max_cv_threshold <- 100  # Set your desired threshold here

# Initialize a dataframe to store CV values
cv_matrix <- data.frame(Protein = rownames(clean_workflows[[1]]))

# Loop through workflows, transform back to linear scale, and compute CV
for (workflow in names(clean_workflows)) {

  dataset <- clean_workflows[[workflow]]

  # Transform data back to linear scale (undo log2 transformation)
  dataset_linear <- 2^dataset

  # Compute CV per protein
  cv_values <- apply(dataset_linear, 1, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100)

  # Add CV values to the dataframe
  cv_matrix[[workflow]] <- cv_values
}

# Set protein names as row names and remove Protein column
rownames(cv_matrix) <- cv_matrix$Protein
cv_matrix <- cv_matrix[, -1]

# Convert the data to long format for ggplot
cv_long <- melt(cv_matrix, variable.name = "Workflow", value.name = "CV")
cv_long$Protein <- rep(rownames(clean_workflows[[1]]), length(names(clean_workflows)))  # Add protein names

# Filter points for labeling (CV > 150)
cv_high <- subset(cv_long[cv_long$CV < 100,], CV > 75)

library(ggplot2)
library(ggrepel)  # For labels

# Generate boxplot
ggplot(cv_long[cv_long$CV < 100,], aes(x = Workflow, y = CV, fill = Workflow)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # Boxplot without outliers
  geom_jitter(aes(color = Workflow), shape = 16, size = 3, alpha = 0.2, width = 0.2) +  # Jitter for visualization
  geom_text_repel(data = cv_high, aes(label = Protein),
                  size = 4, color = "black", max.overlaps = 20,
                  box.padding = 0.5, point.padding = 0.2) +  # Labels for CV > 150
  labs(title = "Distribution of CV Across Workflows",
       x = "",
       y = "Coefficient of Variation (%)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),  # Ensure Y-axis text is also black
    axis.title = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "black"),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black")  # Add x and y axis lines
  ) +
  scale_fill_manual(values = c("MagNet" = "darkgreen",
                               "Neat" = "darkgrey",
                               "PCA" = "darkblue",
                               "Depleted" = "darkred",
                               "Olink" = "darkorange")) +
  scale_color_manual(values = c("MagNet" = "darkgreen",
                                "Neat" = "darkgrey",
                                "PCA" = "darkblue",
                                "Depleted" = "darkred",
                                "Olink" = "darkorange"))  # Match colors for points

