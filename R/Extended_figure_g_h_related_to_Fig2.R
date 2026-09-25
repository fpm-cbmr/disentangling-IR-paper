library(ggplot2)
library(patchwork)
library(dplyr)

processed_datasets <- readRDS("data/processed_datasets.rds")
Sample_ID_filt <- readRDS("data/Sample_ID_filt.rds")
Clinical_data_filt <- readRDS("data/Clinical_data_filt")

Olink_discovery_cohort <- processed_datasets$Olink_discovery_cohort
# Extract the OSM protein values (excluding the ID column)
osm_values <- as.numeric(Olink_discovery_cohort[grepl("OSM_", Olink_discovery_cohort$ID), -1])

# Combine with sample metadata and clinical M-values
df <- data.frame(
  OSM = osm_values,
  Clamp = Sample_ID_filt$Clamp,
  M_value = Clinical_data_filt$M..µmol..kg.min..
)


df$Clamp <- factor(df$Clamp, levels=c("Pre","Post"))
df$Subject_ID <- Sample_ID_filt$Subject_ID

# Scatter plot: OSM vs M-value, by Clamp
scatter_plot_OSM <- ggplot(df, aes(x = M_value, y = OSM, color = Clamp)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "OSM vs M-value",
       x = "M-value (µmol/kg/min)",
       y = NULL) +
  scale_color_manual(values = c("Pre" = "grey", "Post" = "darkred")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.title.x = element_text(color = "black", size = 16),
    axis.ticks.x = element_line(color = "black"),
    # remove y-axis entirely
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = NA)
  )

# Boxplot + paired lines by Subject_ID
boxplot_OSM <- ggplot(df, aes(x = Clamp, y = OSM, group = Clamp)) +
  geom_boxplot(aes(fill = Clamp), alpha = 0.25, width = 0.6, outlier.shape = NA) +
  geom_line(color = "gray80", size = 0.6, aes(group=Subject_ID)) +
  geom_point(aes(color = Clamp), size = 3, alpha = 0.6) +
  scale_fill_manual(values = c("Pre" = "grey", "Post" = "darkred")) +
  scale_color_manual(values = c("Pre" = "grey", "Post" = "darkred")) +
  labs(y = "log2(NPX)") +
  theme_void() +
  theme(
    legend.position = "none",
    axis.title.y = element_text(color = "black", size = 16, angle = 90),
    axis.text.y = element_text(color = "black", size = 14),
    axis.line.y = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black")
  )



Olink_discovery_cohort <- processed_datasets$Olink_discovery_cohort
# Extract the OSM protein values (excluding the ID column)
GCG_values <- as.numeric(Olink_discovery_cohort[grepl("GCG_3", Olink_discovery_cohort$ID), -1])

# Combine with sample metadata and clinical M-values
df <- data.frame(
  GCG = GCG_values,
  Clamp = Sample_ID_filt$Clamp,
  M_value = Clinical_data_filt$M..µmol..kg.min..
)


df$Clamp <- factor(df$Clamp, levels=c("Pre","Post"))
df$Subject_ID <- Sample_ID_filt$Subject_ID

# Scatter plot: GCG vs M-value, by Clamp
scatter_plot_GCG <- ggplot(df, aes(x = M_value, y = GCG, color = Clamp)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "GCG vs M-value",
       x = "M-value (µmol/kg/min)",
       y = NULL) +
  scale_color_manual(values = c("Pre" = "grey", "Post" = "darkred")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.text.x = element_text(color = "black", size = 14),
    axis.title.x = element_text(color = "black", size = 16),
    axis.ticks.x = element_line(color = "black"),
    # remove y-axis entirely
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = NA)
  )

# Boxplot + paired lines by Subject_ID
boxplot_GCG <- ggplot(df, aes(x = Clamp, y = GCG, group = Clamp)) +
  geom_boxplot(aes(fill = Clamp), alpha = 0.25, width = 0.6, outlier.shape = NA) +
  geom_line(color = "gray80", size = 0.6, aes(group=Subject_ID)) +
  geom_point(aes(color = Clamp), size = 3, alpha = 0.6) +
  scale_fill_manual(values = c("Pre" = "grey", "Post" = "darkred")) +
  scale_color_manual(values = c("Pre" = "grey", "Post" = "darkred")) +
  labs(y = "log2(NPX)") +
  theme_void() +
  theme(
    legend.position = "none",
    axis.title.y = element_text(color = "black", size = 16, angle = 90),
    axis.text.y = element_text(color = "black", size = 14),
    axis.line.y = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black")
  )

# Combine plots side by side# Combine plots sideSubject_ID by side
boxplot_OSM + scatter_plot_OSM + boxplot_GCG + scatter_plot_GCG + plot_layout(widths = c(1, 2.5, 1, 2.5))



