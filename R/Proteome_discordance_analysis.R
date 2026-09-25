#### Model the AdipoIR - Mvalue and HOMA-IR - Mvalue discordance ####
library(dplyr)
library(dplyr)
library(tidyr)
library(readr)
library(usethis)
library(ggplot2)
library(ggpubr)
library(tools)
library(stringr)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(ComplexUpset)
library(clusterProfiler)
library(UpSetR)
library(org.Hs.eg.db)


# Full clinical data (ALL samples)
df_ir_all <- Discovery_clinical_data %>%
  mutate(
    AdipoIR_log2 = log2(AdipoIR)
  ) %>%
  dplyr::select(
    -any_of(c(
      "M_discordance",
      "AdipoIR_discordance",
      "M_fitted",
      "AdipoIR_fitted"
    ))
  )

# Baseline subset for modeling
baseline_idx <- Discovery_sample_data$Clamp == "Pre"
df_ir_baseline <- df_ir_all[baseline_idx, ]

# Fit models on baseline only
fit_M <- lm(
  M..µmol..kg.min.. ~ AdipoIR_log2 + Gender + Age + BMI..kg.m2.,
  data = df_ir_baseline,
  na.action = na.exclude
)

fit_Adipo <- lm(
  AdipoIR_log2 ~ M..µmol..kg.min.. + Gender + Age + BMI..kg.m2.,
  data = df_ir_baseline,
  na.action = na.exclude
)

# Add discordance to baseline
df_ir_baseline <- df_ir_baseline %>%
  mutate(
    M_discordance = resid(fit_M),
    AdipoIR_discordance = resid(fit_Adipo),
    M_fitted = fitted(fit_M),
    AdipoIR_fitted = fitted(fit_Adipo)
  )

# Propagate to ALL samples
df_ir <- df_ir_all %>%
  left_join(
    df_ir_baseline %>%
      dplyr::select(
        Subject_ID,
        M_discordance,
        AdipoIR_discordance,
        M_fitted,
        AdipoIR_fitted
      ),
    by = "Subject_ID"
  )





#### Discordance x proteome association ####
library(dplyr)
library(limma)

workflow_names <- c(
  "Neat_discovery_cohort",
  "PCA_discovery_cohort",
  "MagNet_discovery_cohort",
  "Depleted_discovery_cohort",
  "Olink_discovery_cohort"
)

run_limma_assoc <- function(prot, discordance, block, clamp, workflow_name) {
  prot <- as.matrix(prot)
  mode(prot) <- "numeric"

  if (ncol(prot) != length(discordance)) {
    stop(paste0(
      "Sample number mismatch for ", workflow_name,
      ": ncol(prot) = ", ncol(prot),
      ", length(discordance) = ", length(discordance)
    ))
  }

  valid_samples <- is.finite(discordance) & !is.na(block) & !is.na(clamp)

  prot <- prot[, valid_samples, drop = FALSE]
  discordance <- discordance[valid_samples]
  block <- block[valid_samples]
  clamp <- clamp[valid_samples]

  design <- model.matrix(~ discordance + clamp)

  corfit <- duplicateCorrelation(
    prot,
    design = design,
    block = block
  )

  fit <- lmFit(
    prot,
    design = design,
    block = block,
    correlation = corfit$consensus
  )

  fit <- eBayes(fit)

  tt <- topTable(
    fit,
    coef = "discordance",
    number = Inf,
    sort.by = "none"
  )

  tt$Protein <- rownames(tt)
  tt$Workflow <- workflow_name
  tt$N <- rowSums(is.finite(prot))
  tt$ConsensusCorrelation <- corfit$consensus

  tt %>%
    dplyr::select(
      Workflow,
      Protein,
      logFC,
      AveExpr,
      t,
      P.Value,
      adj.P.Val,
      B,
      N,
      ConsensusCorrelation
    )
}

# M discordance x proteome
y_M <- df_ir$M_discordance
block <- df_ir$Subject_ID
df_ir$Clamp <- Discovery_sample_data$Clamp
clamp <- factor(df_ir$Clamp)

results_all_workflows_M_limma <- lapply(workflow_names, function(wf) {
  message("Running M discordance limma: ", wf)

  if (wf == "Olink_discovery_cohort") {
    dat <- processed_datasets[[wf]]
    dat <- as.data.frame(dat)
    rownames(dat) <- dat[[1]]
    dat <- dat[, -1, drop = FALSE]
  } else {
    dat <- processed_datasets[[wf]]$corrected
  }

  run_limma_assoc(
    prot = dat,
    discordance = y_M,
    block = block,
    clamp = clamp,
    workflow_name = wf
  )
})

names(results_all_workflows_M_limma) <- workflow_names
res_M_discordance_all_limma_blocked <- bind_rows(results_all_workflows_M_limma)

# AdipoIR discordance x proteome
y_Adipo <- df_ir$AdipoIR_discordance
block <- df_ir$Subject_ID
clamp <- factor(df_ir$Clamp)

results_all_workflows_Adipo_limma <- lapply(workflow_names, function(wf) {
  message("Running AdipoIR discordance limma: ", wf)

  if (wf == "Olink_discovery_cohort") {
    dat <- processed_datasets[[wf]]
    dat <- as.data.frame(dat)
    rownames(dat) <- dat[[1]]
    dat <- dat[, -1, drop = FALSE]
  } else {
    dat <- processed_datasets[[wf]]$corrected
  }

  run_limma_assoc(
    prot = dat,
    discordance = y_Adipo,
    block = block,
    clamp = clamp,
    workflow_name = wf
  )
})

names(results_all_workflows_Adipo_limma) <- workflow_names
res_AdipoIR_discordance_all_limma_blocked <- bind_rows(results_all_workflows_Adipo_limma)


#### HOMA-IR ####
library(dplyr)
library(limma)

# 1. Prepare full clinical dataset
df_ir_all <- Discovery_clinical_data %>%
  mutate(
    HOMAIR_log2 = log2(HOMA1.IR)
  ) %>%
  dplyr::select(
    -any_of(c(
      "HOMA_discordance",
      "M_discordance_HOMA",
      "HOMA_fitted",
      "M_fitted_HOMA"
    ))
  )

# 2. Baseline subset (for discordance model)
baseline_idx <- Discovery_sample_data$Clamp == "Pre"
df_ir_baseline <- df_ir_all[baseline_idx, ]

# 3. Fit baseline models
# M ~ HOMA
fit_M_HOMA <- lm(
  M..µmol..kg.min.. ~ HOMAIR_log2 + Gender + Age + BMI..kg.m2.,
  data = df_ir_baseline,
  na.action = na.exclude
)

# HOMA ~ M
fit_HOMA <- lm(
  HOMAIR_log2 ~ M..µmol..kg.min.. + Gender + Age + BMI..kg.m2.,
  data = df_ir_baseline,
  na.action = na.exclude
)

# 4. Compute discordance (baseline)
df_ir_baseline <- df_ir_baseline %>%
  mutate(
    M_discordance_HOMA = resid(fit_M_HOMA),
    HOMA_discordance = resid(fit_HOMA),
    M_fitted_HOMA = fitted(fit_M_HOMA),
    HOMA_fitted = fitted(fit_HOMA)
  )

# 5. Propagate to ALL samples
df_ir <- df_ir_all %>%
  left_join(
    df_ir_baseline %>%
      dplyr::select(
        Subject_ID,
        M_discordance_HOMA,
        HOMA_discordance,
        M_fitted_HOMA,
        HOMA_fitted
      ),
    by = "Subject_ID"
  )

# 6. LIMMA FUNCTION
run_limma_assoc <- function(prot, discordance, block, clamp, workflow_name) {
  prot <- as.matrix(prot)
  mode(prot) <- "numeric"

  valid_samples <- is.finite(discordance) & !is.na(block) & !is.na(clamp)

  prot <- prot[, valid_samples, drop = FALSE]
  discordance <- discordance[valid_samples]
  block <- block[valid_samples]
  clamp <- clamp[valid_samples]

  design <- model.matrix(~ discordance + clamp)

  corfit <- duplicateCorrelation(
    prot,
    design = design,
    block = block
  )

  fit <- lmFit(
    prot,
    design = design,
    block = block,
    correlation = corfit$consensus
  )

  fit <- eBayes(fit)

  tt <- topTable(
    fit,
    coef = "discordance",
    number = Inf,
    sort.by = "none"
  )

  tt$Protein <- rownames(tt)
  tt$Workflow <- workflow_name
  tt$N <- rowSums(is.finite(prot))
  tt$ConsensusCorrelation <- corfit$consensus

  tt %>%
    dplyr::select(
      Workflow,
      Protein,
      logFC,
      AveExpr,
      t,
      P.Value,
      adj.P.Val,
      B,
      N,
      ConsensusCorrelation
    )
}

workflow_names <- c(
  "Neat_discovery_cohort",
  "PCA_discovery_cohort",
  "MagNet_discovery_cohort",
  "Depleted_discovery_cohort",
  "Olink_discovery_cohort"
)

block <- df_ir$Subject_ID
df_ir$Clamp <- Discovery_sample_data$Clamp
clamp <- factor(df_ir$Clamp)

# 7. LIMMA: HOMA discordance
y_HOMA <- df_ir$HOMA_discordance

results_HOMA <- lapply(workflow_names, function(wf) {
  message("Running HOMA discordance: ", wf)

  if (wf == "Olink_discovery_cohort") {
    dat <- processed_datasets[[wf]]
    dat <- as.data.frame(dat)
    rownames(dat) <- dat[[1]]
    dat <- dat[, -1, drop = FALSE]
  } else {
    dat <- processed_datasets[[wf]]$corrected
  }

  run_limma_assoc(dat, y_HOMA, block, clamp, wf)
})

names(results_HOMA) <- workflow_names
res_HOMA_discordance_all_limma <- bind_rows(results_HOMA)

# 8. LIMMA: M discordance (HOMA-based)
y_M_HOMA <- df_ir$M_discordance_HOMA

results_M_HOMA <- lapply(workflow_names, function(wf) {
  message("Running M discordance (HOMA-based): ", wf)

  if (wf == "Olink_discovery_cohort") {
    dat <- processed_datasets[[wf]]
    dat <- as.data.frame(dat)
    rownames(dat) <- dat[[1]]
    dat <- dat[, -1, drop = FALSE]
  } else {
    dat <- processed_datasets[[wf]]$corrected
  }

  run_limma_assoc(dat, y_M_HOMA, block, clamp, wf)
})

names(results_M_HOMA) <- workflow_names
res_M_discordance_HOMA_all_limma <- bind_rows(results_M_HOMA)


library(dplyr)

extract_uniprot <- function(id_vec) {
  ifelse(
    grepl("^[^_]+_[^_]+_[^_]+$", id_vec),  # Olink pattern
    sub("^[^_]+_[^_]+_([^_]+)$", "\\1", id_vec),  # take last
    sub("^([^_]+)_.*$", "\\1", id_vec)            # take first (MS)
  )
}

library(dplyr)

pick_top_protein <- function(df) {
  df %>%
    mutate(Uniprot = extract_uniprot(Protein)) %>%
    group_by(Uniprot) %>%
    arrange(adj.P.Val, .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    ungroup()
}
res_M_discordance_top <- pick_top_protein(res_M_discordance_all_limma_blocked)

res_AdipoIR_discordance_top <- pick_top_protein(res_AdipoIR_discordance_all_limma_blocked)

res_HOMA_discordance_top <- pick_top_protein(res_HOMA_discordance_all_limma)

res_M_discordance_HOMA_top <- pick_top_protein(res_M_discordance_HOMA_all_limma)


extract_uniprot_by_workflow <- function(id_vec, workflow_vec) {
  mapply(function(id, wf) {
    parts <- strsplit(id, "_")[[1]]

    if (wf == "Olink_discovery_cohort") {
      tail(parts, 1)
    } else {
      parts[1]
    }
  }, id_vec, workflow_vec, USE.NAMES = FALSE)
}

extract_gene_by_workflow <- function(id_vec, workflow_vec) {
  mapply(function(id, wf) {
    parts <- strsplit(id, "_")[[1]]

    if (wf == "Olink_discovery_cohort") {
      parts[1]
    } else {
      if (length(parts) >= 2) parts[2] else NA_character_
    }
  }, id_vec, workflow_vec, USE.NAMES = FALSE)
}


add_ids <- function(df) {
  df %>%
    mutate(
      Uniprot = extract_uniprot_by_workflow(Protein, Workflow),
      Gene = extract_gene_by_workflow(Protein, Workflow)
    )
}

res_M_discordance_top <- add_ids(res_M_discordance_top)
res_AdipoIR_discordance_top <- add_ids(res_AdipoIR_discordance_top)
res_HOMA_discordance_top <- add_ids(res_HOMA_discordance_top)
res_M_discordance_HOMA_top <- add_ids(res_M_discordance_HOMA_top)


res_M_discordance_top$Is_Significant <- "No"
res_M_discordance_top[res_M_discordance_top$P.Value < 0.01,]$Is_Significant <- "Yes"
res_M_discordance_top$Direction <- ""
res_M_discordance_top[res_M_discordance_top$P.Value < 0.01 & res_M_discordance_top$logFC > 0,]$Direction <- "Up"
res_M_discordance_top[res_M_discordance_top$P.Value < 0.01 & res_M_discordance_top$logFC < 0,]$Direction <- "Down"
res_M_discordance_top_subset <- res_M_discordance_top[,11:14]

res_AdipoIR_discordance_top$Is_Significant <- "No"
res_AdipoIR_discordance_top[res_AdipoIR_discordance_top$P.Value < 0.01,]$Is_Significant <- "Yes"
res_AdipoIR_discordance_top$Direction <- ""
res_AdipoIR_discordance_top[res_AdipoIR_discordance_top$P.Value < 0.01 & res_AdipoIR_discordance_top$logFC > 0,]$Direction <- "Up"
res_AdipoIR_discordance_top[res_AdipoIR_discordance_top$P.Value < 0.01 & res_AdipoIR_discordance_top$logFC < 0,]$Direction <- "Down"
res_AdipoIR_discordance_top_subset <- res_AdipoIR_discordance_top[,11:14]

res_HOMA_discordance_top$Is_Significant <- "No"
res_HOMA_discordance_top[res_HOMA_discordance_top$P.Value < 0.01,]$Is_Significant <- "Yes"
res_HOMA_discordance_top$Direction <- ""
res_HOMA_discordance_top[res_HOMA_discordance_top$P.Value < 0.01 & res_HOMA_discordance_top$logFC > 0,]$Direction <- "Up"
res_HOMA_discordance_top[res_HOMA_discordance_top$P.Value < 0.01 & res_HOMA_discordance_top$logFC < 0,]$Direction <- "Down"
res_HOMA_discordance_top_subset <- res_HOMA_discordance_top[,11:14]





#### GSEA ####
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

make_gsea_vector <- function(res_df, stat_col = "t") {
  # 1) take full ranked results
  df <- res_df %>%
    mutate(
      Uniprot = extract_uniprot(Protein),
      stat = .data[[stat_col]]
    ) %>%
    filter(is.finite(stat), !is.na(Uniprot), Uniprot != "")

  # 2) collapse repeated protein measurements across workflows to one UniProt
  df_uni <- df %>%
    group_by(Uniprot) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup()

  # 3) map UniProt -> Entrez
  map_df <- bitr(
    unique(df_uni$Uniprot),
    fromType = "UNIPROT",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )

  # 4) join back
  df_entrez <- df_uni %>%
    inner_join(map_df, by = c("Uniprot" = "UNIPROT")) %>%
    filter(!is.na(ENTREZID), ENTREZID != "")

  # 5) CRITICAL: collapse again at ENTREZID level
  #    because multiple UniProt IDs can map to the same Entrez gene
  df_entrez2 <- df_entrez %>%
    group_by(ENTREZID) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup()

  gene_vector <- df_entrez2$stat
  names(gene_vector) <- df_entrez2$ENTREZID
  gene_vector <- sort(gene_vector, decreasing = TRUE)

  return(gene_vector)
}

run_gobp <- function(res_df, stat_col = "t") {
  gene_vector <- make_gsea_vector(res_df, stat_col = stat_col)

  gseGO(
    geneList = gene_vector,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = FALSE,
    by = "fgsea"
  )
}

gobp_M <- run_gobp(res_M_discordance_all_limma_blocked)
gobp_Adipo <- run_gobp(res_AdipoIR_discordance_all_limma_blocked)
gobp_HOMA <- run_gobp(res_HOMA_discordance_all_limma)
gobp_M_HOMA <- run_gobp(res_M_discordance_HOMA_all_limma)

gobp_M_simpl <- simplify(
  gobp_M,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

gobp_Adipo_simpl <- simplify(
  gobp_Adipo,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

gobp_HOMA_simpl <- simplify(
  gobp_HOMA,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

gobp_M_HOMA_simpl <- simplify(
  gobp_M_HOMA,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)


gobp_compare_simpl <- bind_rows(
  #gobp_M_simpl@result %>% mutate(Analysis = "M_discordance"),
  gobp_Adipo_simpl@result %>% mutate(Analysis = "AdipoIR_discordance"),
  gobp_HOMA_simpl@result %>% mutate(Analysis = "HOMA_discordance")
  #gobp_M_HOMA_simpl@result %>% mutate(Analysis = "M_discordance_HOMA")
) %>%
  dplyr::select(ID, Description, Analysis, NES, p.adjust)

gobp_sig_simpl <- gobp_compare_simpl %>%
  group_by(ID, Description) %>%
  mutate(any_sig = any(p.adjust < 0.05, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(any_sig)

top_pathways_simpl <- gobp_sig_simpl %>%
    group_by(ID, Description) %>%
    summarise(min_padj = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
    arrange(min_padj) %>%
    slice_head(n = 20)

plot_df <- gobp_sig_simpl %>%
  semi_join(top_pathways_simpl, by = c("ID", "Description")) %>%
  mutate(sig_star = ifelse(p.adjust < 0.05, "*", ""))

pathway_order <- plot_df %>%
  group_by(Description) %>%
  summarise(max_abs_NES = max(abs(NES), na.rm = TRUE), .groups = "drop") %>%
  arrange(max_abs_NES) %>%
  pull(Description)

plot_df <- plot_df %>%
  mutate(Description = factor(Description, levels = pathway_order))

# the code for the plot used for the manuscript is at the end #
ggplot(plot_df, aes(x = Analysis, y = Description, fill = NES)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sig_star), size = 5) +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    name = "NES"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.line = element_blank()
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Simplified GO BP pathway comparison across discordance analyses"
  )


#### GSEA GOBP ####


library(dplyr)
library(ggplot2)

# Keep pathways significant in at least one comparison
gobp_sig_simpl <- gobp_compare_simpl %>%
    group_by(ID, Description) %>%
    mutate(any_sig = any(p.adjust < 0.05, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(any_sig)

# Top 10 per Analysis based on min adjusted p-value
top10_per_analysis <- gobp_sig_simpl %>%
    group_by(Analysis, ID, Description) %>%
    summarise(
        min_padj = min(p.adjust, na.rm = TRUE),
        max_abs_NES = max(abs(NES), na.rm = TRUE),
        .groups = "drop"
    ) %>%
    group_by(Analysis) %>%
    arrange(min_padj, desc(max_abs_NES)) %>%
    slice_head(n = 10) %>%
    ungroup()

plot_df <- gobp_sig_simpl %>%
    semi_join(top10_per_analysis, by = c("Analysis", "ID", "Description")) %>%
    mutate(sig_star = ifelse(p.adjust < 0.05, "*", ""))





plot_one_analysis <- function(df, title_text = NULL, n_terms = 10) {

    df_pos <- df %>%
        filter(NES > 0) %>%
        mutate(
            neglog10_padj = pmin(-log10(p.adjust), 10)
        )

    top_terms <- df_pos %>%
        group_by(ID, Description) %>%
        summarise(
            min_padj = min(p.adjust, na.rm = TRUE),
            max_NES = max(NES, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        arrange(min_padj, desc(max_NES))

    top_terms <- top_terms %>%
        dplyr::slice(seq_len(min(n_terms, nrow(top_terms))))

    plot_df2 <- df_pos %>%
        inner_join(top_terms, by = c("ID", "Description"))

    pathway_order <- top_terms %>%
        arrange(min_padj, desc(max_NES)) %>%
        pull(Description)

    plot_df2 <- plot_df2 %>%
        mutate(
            Description = factor(Description, levels = rev(pathway_order)),
            x = 1
        )

    line_df <- plot_df2 %>%
        distinct(Description) %>%
        mutate(x = 0, xend = 1)

    ggplot(plot_df2, aes(x = x, y = Description)) +
        geom_segment(
            data = line_df,
            aes(x = x, xend = xend, y = Description, yend = Description),
            inherit.aes = FALSE,
            color = "grey80",
            linewidth = 0.5
        ) +
        geom_point(
            aes(size = NES, color = neglog10_padj),
            shape = 16
        ) +
        scale_color_gradient(
            low = "lightgrey",
            high = "red",
            name = expression(-log[10](FDR))
        ) +
        scale_size_continuous(
            range = c(3, 10),
            name = "NES"
        ) +
        scale_x_continuous(
            limits = c(0, 1.15),
            expand = c(0, 0)
        ) +
        theme_classic(base_size = 14) +
        theme(
            axis.text.y = element_text(color = "black"),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line = element_blank()
        ) +
        labs(
            x = NULL,
            y = NULL,
            title = title_text
        )
}

plot_df_list <- split(gobp_compare_simpl, gobp_compare_simpl$Analysis)

p1 <- plot_one_analysis(plot_df_list[[1]], names(plot_df_list)[1], n_terms = 10)
p2 <- plot_one_analysis(plot_df_list[[2]], names(plot_df_list)[2], n_terms = 10)

p1
p2

#### Discordance proteome plotting ####
library(dplyr)
library(ggplot2)
library(ggrepel)

# M discordance: ranked protein plot
plot_M <- res_M_discordance_top %>%
    filter(!is.na(adj.P.Val), !is.na(logFC)) %>%
    arrange(adj.P.Val) %>%
    mutate(
        Rank = row_number(),
        Sig = adj.P.Val < 0.05,
        Label = ifelse(Sig, Gene, NA_character_)
    )



ggplot(plot_M, aes(x = Rank, y = -log10(adj.P.Val))) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(aes(color = Sig), size = 2) +
    geom_text_repel(
        aes(label = Label),
        size = 3,
        max.overlaps = 30,
        na.rm = TRUE
    ) +
    scale_color_manual(values = c("FALSE" = "grey75", "TRUE" = "#b2182b")) +
    theme_classic(base_size = 14) +
    theme(
        axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black")
    ) +
    labs(
        x = "Ranked proteins",
        y = "Effect size (logFC)",
        color = "adj. P < 0.05",
        title = "Proteins associated with M discordance"
    )

# AdipoIR discordance: ranked protein plot
plot_Adipo <- res_AdipoIR_discordance_top %>%
    filter(!is.na(adj.P.Val), !is.na(logFC)) %>%
    arrange(adj.P.Val) %>%
    mutate(
        Rank = row_number(),
        Sig = adj.P.Val < 0.05,
        Label = ifelse(Sig, Gene, NA_character_)
    )

plot_Adipo$DirectionSig <- ifelse(
    plot_Adipo$Sig & plot_Adipo$logFC > 0, "Positive",
    ifelse(plot_Adipo$Sig & plot_Adipo$logFC < 0, "Negative", "NS")
)
ggplot(plot_Adipo, aes(x = Rank, y = -log10(adj.P.Val))) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_point(
        aes(color = DirectionSig, size = Sig),
        alpha = 0.25,
        stroke = NA
    ) +
    geom_text_repel(
        aes(label = Label),
        size = 3,
        max.overlaps = 30,
        na.rm = TRUE
    ) +
    scale_color_manual(
        values = c(
            "NS" = "grey75",
            "Positive" = "#b2182b",  # red
            "Negative" = "#2166ac"   # blue
        )
    ) +
    scale_size_manual(values = c("FALSE" = 1.5, "TRUE" = 3)) +
    guides(size = "none") +
    theme_classic(base_size = 14) +
    theme(
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", color = NA),
        axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black")
    ) +
    labs(
        x = "Ranked proteins",
        y = "-log10(adj.P)",
        color = "Association",
        title = "AdipoIR discordance"
    )


plot_HOMA <- res_HOMA_discordance_top %>%
    filter(!is.na(adj.P.Val), !is.na(logFC)) %>%
    arrange(adj.P.Val) %>%
    mutate(
        Rank = row_number(),
        Sig = adj.P.Val < 0.05,
        Label = ifelse(Sig, Gene, NA_character_)
    )


plot_HOMA$DirectionSig <- ifelse(
    plot_HOMA$Sig & plot_M$logFC > 0, "Positive",
    ifelse(plot_HOMA$Sig & plot_HOMA$logFC < 0, "Negative", "NS")
)
ggplot(plot_HOMA, aes(x = Rank, y = -log10(adj.P.Val))) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_point(
        aes(color = DirectionSig, size = Sig),
        alpha = 0.25,
        stroke = NA
    ) +
    geom_text_repel(
        aes(label = Label),
        size = 3,
        max.overlaps = 30,
        na.rm = TRUE
    ) +
    scale_color_manual(
        values = c(
            "NS" = "grey75",
            "Positive" = "#b2182b",  # red
            "Negative" = "#2166ac"   # blue
        )
    ) +
    scale_size_manual(values = c("FALSE" = 1.5, "TRUE" = 3)) +
    guides(size = "none") +
    theme_classic(base_size = 14) +
    theme(
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", color = NA),
        axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black")
    ) +
    labs(
        x = "Ranked proteins",
        y = "-log10(adj.P)",
        color = "Association",
        title = "HOMA-IR discordance"
    )
