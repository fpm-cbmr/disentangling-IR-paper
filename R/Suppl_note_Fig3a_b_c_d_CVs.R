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

# Step 1: Get a common sorted list of gene names
common_gene_names <- sort(intersect(
    rownames(Common_genes_MagNet),
    intersect(rownames(Common_genes_PCA),
              intersect(rownames(Common_genes_Neat),
                        intersect(rownames(Common_genes_Depleted),
                                  rownames(Common_genes_Olink))))
))

# Step 2: Reorder each dataset to match this order
Common_genes_MagNet <- Common_genes_MagNet[common_gene_names, , drop = FALSE]
Common_genes_PCA <- Common_genes_PCA[common_gene_names, , drop = FALSE]
Common_genes_Neat <- Common_genes_Neat[common_gene_names, , drop = FALSE]
Common_genes_Depleted <- Common_genes_Depleted[common_gene_names, , drop = FALSE]
Common_genes_Olink <- Common_genes_Olink[common_gene_names, , drop = FALSE]

# Step 3: Store in a list for downstream analysis
all_workflows <- list(
    Common_genes_Olink,
    Common_genes_PCA,
    Common_genes_Depleted,
    Common_genes_Neat,
    Common_genes_MagNet
)
names(all_workflows) <- c("Olink", "PCA", "Depleted", "Neat", "MagNet")

workflow_names <- names(all_workflows)
num_workflows  <- length(workflow_names)



# Find proteins (rows) that have missing values in any workflow
common_na_mask <- Reduce(`|`, lapply(all_workflows, function(x) rowSums(is.na(x)) > 0))

# Find protein names that should be kept (only complete proteins)
complete_proteins <- rownames(all_workflows[[1]])[!common_na_mask]

# Verify the number of proteins left
length(complete_proteins)  # Should now be consistent across all workflows



# Subset all workflows to only include proteins without missing values
clean_workflows <- lapply(all_workflows, function(x) x[complete_proteins, , drop = FALSE])

# Verify that all workflows now have the same dimensions
sapply(clean_workflows, dim)  # Should now return the same number of proteins × 150 samples




# Function to scale each protein within a workflow (row-wise across samples)
scale_workflow <- function(matrix) {
    t(apply(matrix, 1, scale))  # Z-score normalization across samples
}

# Apply scaling to each workflow
scaled_workflows <- lapply(clean_workflows, scale_workflow)

# Convert to data frames to allow correct indexing
scaled_workflows <- lapply(scaled_workflows, as.data.frame)

# Initialize correlation matrix
workflow_correlation_matrix <- matrix(NA, ncol = num_workflows, nrow = num_workflows)
colnames(workflow_correlation_matrix) <- rownames(workflow_correlation_matrix) <- workflow_names
for (i in 1:(num_workflows - 1)) {
    for (j in (i + 1):num_workflows) {
        print(paste("Computing correlation for:", workflow_names[i], "vs", workflow_names[j]))

        # Get common proteins and ensure they are sorted in the same order
        common_proteins <- sort(intersect(rownames(scaled_workflows[[i]]), rownames(scaled_workflows[[j]])))

        # Check if there are common proteins
        if (length(common_proteins) == 0) {
            print("⚠ No common proteins! Skipping this comparison.")
            next
        }

        # Extract numeric matrices and ensure rownames are explicitly set in the same order
        x <- as.matrix(scaled_workflows[[i]][common_proteins, , drop = FALSE])
        y <- as.matrix(scaled_workflows[[j]][common_proteins, , drop = FALSE])

        # Ensure rownames match in both matrices
        rownames(y) <- rownames(x)  # Force y to have the same row order as x

        # Check if rownames are actually aligned (fail if not)
        if (!all(rownames(x) == rownames(y))) {
            stop("❌ ERROR: Row names of x and y are still misaligned!")
        }

        # Compute correlation row-wise using mapply()
        cor_values <- mapply(function(row_x, row_y) {
            if (sum(!is.na(row_x)) > 2 && sum(!is.na(row_y)) > 2) {
                return(cor.test(row_x, row_y, method = "spearman", use = "pairwise.complete.obs")$estimate)
            } else {
                return(NA)
            }
        }, split(x, rownames(x)), split(y, rownames(y)))

        # Store mean correlation
        workflow_correlation_matrix[i, j] <- mean(cor_values, na.rm = TRUE)
        workflow_correlation_matrix[j, i] <- workflow_correlation_matrix[i, j]  # Make symmetric

        # Debugging: Print first few computed correlations
        print(paste("First computed correlation:", head(cor_values, 1)))
    }
}


# Convert correlation matrix to long format
cor_df <- melt(workflow_correlation_matrix)

# Rename columns for clarity
colnames(cor_df) <- c("Workflow1", "Workflow2", "Correlation")

# Keep only the lower triangle and **remove self-correlations (diagonal)**
cor_df <- cor_df[as.numeric(cor_df$Workflow1) > as.numeric(cor_df$Workflow2), ]

# Plot the Lower Triangle Correlation Heatmap
ggplot(cor_df, aes(x = Workflow1, y = Workflow2, fill = Correlation)) +
  geom_tile(color = "white") +  # White grid lines
  geom_text(aes(label = round(Correlation, 2)), color = "black", size = 6) +  # Show correlation values
  scale_fill_gradientn(colors = c("white", "darkred"),
                       limits = c(0.2, 0.8),  # Easy control of color range
                       breaks = seq(0, 1, by = 0.1),  # Customize tick marks
                       name = "Correlation") +  # Legend title
  scale_y_discrete(position = "right") +  # Move workflow labels to the right
  theme_minimal() +  # Clean theme
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(0.05, 0.95),  # **Move legend to top-left corner**
        legend.justification = c(0, 1),  # Align legend inside the plot
        legend.background = element_rect(fill = "white", color = "black")) +  # Optional: add a background box
  labs(title = "Workflow Correlation Heatmap (Lower Triangle)")









library(ggplot2)
library(reshape2)

# Create an empty list to store correlation values for each workflow pair
correlation_data <- list()

# Loop over each workflow pair to extract individual protein-wise correlations
for (i in 1:(num_workflows - 1)) {
    for (j in (i + 1):num_workflows) {
        cor_values <- sapply(1:155, function(p) {
            x <- as.numeric(scaled_workflows[[i]][p, ])
            y <- as.numeric(scaled_workflows[[j]][p, ])

            if (length(x) == 150 && length(y) == 150 && sum(!is.na(x)) > 2 && sum(!is.na(y)) > 2) {
                return(cor(x, y, method = "spearman", use = "pairwise.complete.obs"))
            } else {
                return(NA)
            }
        })

        # Extract protein names from rownames of scaled_workflows
        protein_names <- rownames(scaled_workflows[[i]])[1:155]

        # Store in a data frame
        correlation_data[[paste(workflow_names[i], "vs", workflow_names[j])]] <- data.frame(
            Workflow_Pair = paste(workflow_names[i], "vs", workflow_names[j]),
            Protein = protein_names,
            Correlation = cor_values
        )
    }
}

# Combine all workflow pairs into one data frame
correlation_df <- do.call(rbind, correlation_data)


# Compute mean correlation for each workflow pair and reorder factor levels
correlation_df$Workflow_Pair <- factor(correlation_df$Workflow_Pair,
                                       levels = names(sort(tapply(correlation_df$Correlation, correlation_df$Workflow_Pair, mean), decreasing = TRUE)))

ggplot(correlation_df, aes(x = Workflow_Pair, y = Correlation)) +
  # Transparent bars showing the mean correlation
  stat_summary(fun = mean, geom = "bar", alpha = 0.4, fill = "darkgrey", color = NA) +

  # Scatter plot of individual correlations
  geom_jitter(aes(color = Correlation < 0), width = 0.2, alpha = 0.3, size = 1) +

  # Add a horizontal dashed line at y = 0
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8) +

  # Set colors: Negative values = blue, Positive values = darkred
  scale_color_manual(values = c("darkred", "blue")) +

  # Aesthetics
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, color = "black"),  # X-axis labels black
    axis.text.y = element_text(size = 12, color = "black"),  # Y-axis labels black
    axis.title.x = element_text(size = 14, color = "black"),  # X-axis title black
    axis.title.y = element_text(size = 14, color = "black"),  # Y-axis title black
    panel.grid.major.x = element_blank(), aspect.ratio=1.4/1,
    legend.position = "none",  # Hide legend for dot colors
    axis.line = element_line(color = "black")  # **Black X and Y axis lines**
  ) +

  # Labels
  labs(y = "Protein-wise Correlation Coefficients",
       title = "Distribution of Correlations Across Workflows (Ranked by Mean Correlation)")


#### Examples of correlations/anti-correlations ####

library(ggplot2)
library(dplyr)

# Filter data for GC and LPA in the selected workflow pairs
library(ggplot2)
library(gridExtra)

# Define workflow indices
olink_idx <- 1
neat_idx <- 4
depleted_idx <- 3

# Extract raw scaled values for GC
gc_neat <- data.frame(
    Workflow = "Olink vs Neat",
    Protein = "GC",
    Olink = as.numeric(scaled_workflows[[olink_idx]]["GC", ]),
    Other = as.numeric(scaled_workflows[[neat_idx]]["GC", ])
)

gc_depleted <- data.frame(
    Workflow = "neat vs Depleted",
    Protein = "GC",
    Olink = as.numeric(scaled_workflows[[neat_idx]]["GC", ]),
    Other = as.numeric(scaled_workflows[[depleted_idx]]["GC", ])
)

# Extract raw scaled values for LPA
lpa_neat <- data.frame(
    Workflow = "Olink vs Neat",
    Protein = "LPA",
    Olink = as.numeric(scaled_workflows[[olink_idx]]["LPA", ]),
    Other = as.numeric(scaled_workflows[[neat_idx]]["LPA", ])
)

lpa_depleted <- data.frame(
    Workflow = "neat vs Depleted",
    Protein = "LPA",
    Olink = as.numeric(scaled_workflows[[neat_idx]]["LPA", ]),
    Other = as.numeric(scaled_workflows[[depleted_idx]]["LPA", ])
)

# Combine into one dataset for plotting
correlation_plot_data <- rbind(gc_neat, gc_depleted, lpa_neat, lpa_depleted)

# Print to verify extraction
print(head(correlation_plot_data))


library(ggplot2)
library(gridExtra)

# Base theme modification
custom_theme <- theme_minimal() +
    theme(
        axis.line = element_line(color = "black"),  # Black axes
        axis.text = element_text(size = 12, color = "black"),  # Black tick labels
        axis.title = element_text(size = 14, color = "black"),  # Black axis titles
        plot.title = element_text(size = 14, face = "bold", color = "black"),  # Black title
        panel.grid = element_blank(),  # Remove grid lines
        panel.border = element_blank()  # Remove border around plot
    )

# Scatter plot for GC: Olink vs Neat
plot_gc_neat <- ggplot(gc_neat[gc_neat$Olink > -2.5,], aes(x = Olink, y = Other)) +
    geom_point(color = "blue", alpha = 0.6) +
    labs(title = "GC: Olink vs Neat", x = "Olink (Scaled)", y = "Neat (Scaled)") +
    custom_theme

# Scatter plot for GC: Olink vs Depleted
plot_gc_depleted <- ggplot(gc_depleted, aes(x = Olink, y = Other)) +
    geom_point(color = "darkred", alpha = 0.6) +
    labs(title = "GC: Neat vs Depleted", x = "Neat (Scaled)", y = "Depleted (Scaled)") +
    custom_theme

# Scatter plot for LPA: Olink vs Neat
plot_lpa_neat <- ggplot(lpa_neat, aes(x = Olink, y = Other)) +
    geom_point(color = "darkred", alpha = 0.6) +
    labs(title = "LPA: Olink vs Neat", x = "Olink (Scaled)", y = "Neat (Scaled)") +
    custom_theme

# Scatter plot for LPA: Olink vs Depleted
plot_lpa_depleted <- ggplot(lpa_depleted, aes(x = Olink, y = Other)) +
    geom_point(color = "darkred", alpha = 0.6) +
    labs(title = "LPA: Neat vs Depleted", x = "Neat (Scaled)", y = "Depleted (Scaled)") +
    custom_theme


# Arrange all 4 plots in a 2x2 grid
grid.arrange(plot_lpa_neat,plot_gc_neat, ncol = 1)
grid.arrange(plot_lpa_depleted,plot_gc_depleted, ncol = 1)


















# Load in core datasets
library(VennDiagram)
library(purrr)
library(tidyverse)
library(PhosR)
processed_datasets <- readRDS("data/processed_datasets.rds")

Olink_proteins <- unique(word(processed_datasets$Olink_discovery_cohort$ID,1,sep="_"))
MagNet_proteins <- unique(word(word(rownames(processed_datasets$MagNet_discovery_cohort$raw),2,sep="_"),1,sep=";"))
PCA_proteins <- unique(word(word(rownames(processed_datasets$PCA_discovery_cohort$raw),2,sep="_"),1,sep=";"))
Neat_proteins <- unique(word(word(rownames(processed_datasets$Neat_discovery_cohort$raw),2,sep="_"),1,sep=";"))
Depleted_proteins <- unique(word(word(rownames(processed_datasets$Depleted_discovery_cohort$raw),2,sep="_"),1,sep=";"))


#### Workflow stacked unique protein gain identifications ####

library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)
library(stringr)


# Put into list in desired order
workflow_proteins <- list(
    Neat = Neat_proteins,
    PCA = PCA_proteins,
    Depleted = Depleted_proteins,
    MagNet = MagNet_proteins,
    Olink = Olink_proteins
)

# Initialize empty data frame
gain_df <- data.frame()

# Track cumulative union of previous workflows
cumulative_set <- character(0)

# Calculate gain per workflow
for (workflow in names(workflow_proteins)) {
    current_set <- workflow_proteins[[workflow]]
    new_proteins <- setdiff(current_set, cumulative_set)
    gain_df <- rbind(
        gain_df,
        data.frame(Workflow = workflow, New_Proteins = length(new_proteins))
    )
    cumulative_set <- union(cumulative_set, current_set)
}

# Force factor order
gain_df$Workflow <- factor(gain_df$Workflow, levels = c("Olink","MagNet","Depleted","PCA","Neat"))

ggplot(gain_df, aes(x = "Proteome", y = New_Proteins, fill = Workflow)) +
    geom_bar(stat = "identity", width = 0.5, alpha = 0.4) +
    geom_text(
        aes(label = Workflow),
        position = position_stack(vjust = 0.5),
        color = "black",
        size = 6
    ) +
    labs(
        y = "Cumulative Gain in Unique Proteins",
        x = NULL,
        title = NULL
    ) +
    scale_fill_manual(values = c(
        "Neat" = "darkgrey",
        "PCA" = "darkblue",
        "Depleted" = "darkred",
        "MagNet" = "darkgreen",
        "Olink" = "darkorange"
    )) +
    theme_bw(base_size = 26) +
    theme(
        panel.grid = element_blank(),               # Remove all grid lines
        axis.text = element_text(color = "black"),  # Make axis text black
        axis.ticks = element_line(color = "black"), # Optional: black axis ticks
        axis.text.x = element_blank(),              # Remove x-axis text
        axis.ticks.x = element_blank(),             # Remove x-axis ticks
        legend.position = "none"                    # Remove legend
    )




#### upset plot ####
# Load required libraries
library(ggplot2)
library(ComplexUpset)
library(tidyr)
library(dplyr)

# Same five sets as above, in the order the plot should display them
protein_lists <- workflow_proteins[c("Olink", "MagNet", "PCA", "Neat", "Depleted")]

all_proteins <- unique(unlist(protein_lists, use.names = FALSE))

upset_data <- cbind(
    Protein = all_proteins,
    as.data.frame(lapply(protein_lists, function(s) all_proteins %in% s))
)


# Define unique colors for each workflow
workflow_colors <- c("Olink" = "darkorange",
                     "MagNet" = "darkgreen",
                     "PCA" = "darkblue",
                     "Neat" = "darkgrey",
                     "Depleted" = "darkred")

# Create the UpSet plot without set size bars
upset(
    upset_data,
    intersect=c("Olink", "MagNet", "PCA", "Neat", "Depleted"),  # Workflow columns
    name="Workflow",
    width_ratio=0.1,

    # Custom intersection matrix with colored dots
    matrix=(
        intersection_matrix(geom=geom_point(shape='circle filled', size=4))
        + scale_color_manual(values=workflow_colors,
                             guide=guide_legend(override.aes=list(shape='circle')))
    ),

    # Remove set size bars
    set_sizes = FALSE,

    # Highlighting workflow sets with specific colors
    queries=list(
        upset_query(set="Olink", fill="darkorange", color="black"),
        upset_query(set="MagNet", fill="darkgreen",color="black"),
        upset_query(set="PCA", fill="darkblue",color="black"),
        upset_query(set="Neat", fill="darkgrey",color="black"),
        upset_query(set="Depleted", fill="darkred",color="black")
    )
) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 12, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          legend.position = "none")

