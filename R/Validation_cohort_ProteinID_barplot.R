processed_datasets_validation <- readRDS("data/processed_datasets_validation.rds")

library(stringr)
library(dplyr)
library(ggplot2)
library(tidyr)

# Group MS workflows together
workflow_groups <- list(
  `MS-based` = c("Neat_validation_cohort", "MagNet_validation_cohort",
                 "PCA_validation_cohort", "Depleted_validation_cohort"),
  `Olink` = c("Olink_validation_cohort")
)

# Extract gene lists and assign to new group labels
gene_lists_grouped <- list()

for (group in names(workflow_groups)) {
  genes <- unlist(lapply(workflow_groups[[group]], function(wf) {
    if (wf != "Olink_validation_cohort") {
      word(word(rownames(processed_datasets_validation[[wf]]$raw), 2, sep = "_"), 1, sep = ";")
    } else {
      word(rownames(processed_datasets_validation[[wf]]), 1, sep = "_")
    }
  }))
  gene_lists_grouped[[group]] <- unique(genes)
}

# Long format: Gene - Group
gene_df <- enframe(gene_lists_grouped, name = "Workflow", value = "Gene") %>%
  unnest(Gene)

# Count in how many groups each gene appears (MS-based vs Olink)
gene_counts <- gene_df %>%
  distinct(Workflow, Gene) %>%
  count(Gene, name = "n_workflows")

# Classify genes as Shared or Unique (across MS vs Olink)
annotated <- gene_df %>%
  left_join(gene_counts, by = "Gene") %>%
  mutate(Category = ifelse(n_workflows > 1, "Shared", "Unique"))

# Summarise for plotting
bar_data <- annotated %>%
  group_by(Workflow, Category) %>%
  summarise(n = n_distinct(Gene), .groups = "drop")

# Plot
ggplot(bar_data, aes(x = Workflow, y = n, fill = Category)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("Shared" = "lightgrey", "Unique" = "#707070")) +
  labs(
    title = "Protein Overlap: MS-based vs Olink",
    y = "Number of Proteins", x = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),                    # Remove gridlines
    axis.line = element_line(color = "black"),       # Add black x and y axes
    axis.ticks = element_line(color = "black"),      # Optional: also make ticks black
    axis.text = element_text(color = "black"),       # Make axis text black
    axis.title = element_text(color = "black"),
    legend.position = "top",
    plot.title = element_text(color = "black", face = "bold", hjust = 0.5)
  )

