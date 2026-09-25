library(dplyr)
library(tidyr)
library(readr)
library(stringr)
######################################
#### Tissue network mappin ####
######################################
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

# --- 1) Membership of proteins per group
membership <- bind_rows(
    tibble(GeneID = new_lists$m_positive, group = "m_positive"),
    tibble(GeneID = new_lists$m_negative, group = "m_negative"),
    tibble(GeneID = new_lists$m_all,      group = "m_all")
) %>% distinct()
view(tissue_map)
# --- 2) Tissue map: one row per GeneID–tissue
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

# --- 3) Significant tissues (you can tweak thresholds)
padj_cut <- 0.05
min_or   <- 1

sig_tissues <- combined_tissue_results %>%
    filter(!is.na(padj), padj < padj_cut, fisher_OR > min_or) %>%
    mutate(tissue = str_trim(tolower(tissue))) %>%
    dplyr::select(group, tissue, fisher_OR, padj) %>%
    distinct()

# --- 4) Proteins per (group, tissue): done via joins (no rowwise)
sig_tissue_proteins_long <- sig_tissues %>%
    inner_join(tissue_map,  by = "tissue") %>%                # add GeneID for that tissue
    inner_join(membership,  by = c("group", "GeneID")) %>%    # keep only genes in that group
    dplyr::select(group, tissue, GeneID, Uniprot, fisher_OR, padj) %>%
    distinct() %>%
    arrange(group, tissue, GeneID)

# Optional: summary counts per (group, tissue)
sig_tissue_counts <- sig_tissue_proteins_long %>%
    count(group, tissue, name = "n_proteins") %>%
    arrange(group, desc(n_proteins))

# Save
write_csv(sig_tissue_proteins_long, "significant_tissue_proteins.csv")
write_csv(sig_tissue_counts,        "significant_tissue_protein_counts.csv")

sig_tissue_proteins_long <- read.csv("significant_tissue_proteins.csv")



# --- 1) Membership of proteins per group
membership <- bind_rows(
  tibble(GeneID = new_lists$m_positive, group = "m_positive"),
  tibble(GeneID = new_lists$m_negative, group = "m_negative"),
  tibble(GeneID = new_lists$m_all,      group = "m_all")
) %>% distinct()
view(tissue_map)
# --- 2) Tissue map: one row per GeneID–tissue
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

# --- 3) Significant tissues (you can tweak thresholds)
padj_cut <- 0.05
min_or   <- 1

sig_tissues <- combined_tissue_results %>%
  filter(!is.na(padj), padj < padj_cut, fisher_OR > min_or) %>%
  mutate(tissue = str_trim(tolower(tissue))) %>%
  dplyr::select(group, tissue, fisher_OR, padj) %>%
  distinct()

# --- 4) Proteins per (group, tissue): done via joins (no rowwise)
sig_tissue_proteins_long <- sig_tissues %>%
  inner_join(tissue_map,  by = "tissue") %>%                # add GeneID for that tissue
  inner_join(membership,  by = c("group", "GeneID")) %>%    # keep only genes in that group
  dplyr::select(group, tissue, GeneID, Uniprot, fisher_OR, padj) %>%
  distinct() %>%
  arrange(group, tissue, GeneID)

# Optional: summary counts per (group, tissue)
sig_tissue_counts <- sig_tissue_proteins_long %>%
  count(group, tissue, name = "n_proteins") %>%
  arrange(group, desc(n_proteins))

# Save
write_csv(sig_tissue_proteins_long, "significant_tissue_proteins.csv")
write_csv(sig_tissue_counts,        "significant_tissue_protein_counts.csv")

sig_tissue_proteins_long <- read.csv("significant_tissue_proteins.csv")




# =========================
# Setup & data loading
# =========================
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(readr)

Discovery_proteome_data <- readRDS("data/processed_datasets.rds")
Discovery_clinical_data <- readRDS("data/Clinical_data_filt")
Discovery_sample_data   <- readRDS("data/Sample_ID_filt.rds")

# =========================
# Detect sample IDs across proteome matrices (all workflows)
# =========================
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

# =========================
# Auto-detect metadata ID column & build M-values for PRE
# =========================
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
    M_value  = as.numeric(`M..µmol..kg.min..`)
  ) %>%
  filter(!is.na(SampleID), !is.na(M_value))

keep_cols <- intersect(all_samples, meta_with_M$SampleID)
stopifnot(length(keep_cols) > 0)

# =========================
# Subset each workflow to PRE samples by name
# =========================
subset_ms    <- function(mat) mat[, intersect(colnames(mat), keep_cols), drop = FALSE]
subset_olink <- function(df)  df[, c("ID", intersect(colnames(df), keep_cols)), drop = FALSE]

Neat_mat_pre     <- subset_ms(ms_mats$Neat)
PCA_mat_pre      <- subset_ms(ms_mats$PCA)
MagNet_mat_pre   <- subset_ms(ms_mats$MagNet)
Depleted_mat_pre <- subset_ms(ms_mats$Depleted)
Olink_df_pre     <- subset_olink(olink_df)

M_values_df <- meta_with_M

# =========================
# Helpers to parse IDs & build long tables
# =========================
parse_ms_id <- function(x) {
  tibble(
    id     = x,
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

# =========================

library(dplyr)
library(ggplot2)
library(patchwork)

library(dplyr)
library(stringr)

# Clean expr_long
expr_long <- expr_long %>%
    mutate(
        Workflow = str_trim(Workflow),
        Gene     = str_trim(Gene),
        Uniprot  = toupper(str_trim(Uniprot))
    )

# Map raw workflow labels in unique_significant to canonical names used in expr_long
sig_workflow <- unique_significant %>%
    mutate(
        Gene = str_trim(ID),
        Workflow_raw = str_trim(Workflow),
        Workflow = case_when(
            Workflow_raw %in% c("Olink", "OLINK", "Olink_discovery") ~ "Olink_discovery_cohort",
            Workflow_raw %in% c("Depleted", "DEP", "Depleted_discovery") ~ "Depleted_discovery_cohort",
            Workflow_raw %in% c("MagNet", "MAGNET") ~ "MagNet_discovery_cohort",
            Workflow_raw %in% c("PCA") ~ "PCA_discovery_cohort",
            Workflow_raw %in% c("Neat") ~ "Neat_discovery_cohort",
            TRUE ~ Workflow_raw        # if already in canonical form
        )
    ) %>%
    filter(
        !is.na(Workflow), Workflow != "",
        !is.na(Gene), Gene != ""
    ) %>%
    distinct(Gene, Workflow)

sig_workflow %>% filter(Gene == "LPL")
expr_long %>% filter(Gene == "LPL") %>% count(Workflow)



library(ggplot2)

plot_gene <- function(gene_symbol) {
    gene_symbol <- str_trim(gene_symbol)

    # All rows for this gene
    dat_gene <- expr_long %>%
        filter(Gene == gene_symbol)

    if (nrow(dat_gene) == 0) {
        warning("No data available for gene: ", gene_symbol)
        return(NULL)
    }

    # Lookup discovery workflow from sig_workflow
    wf_row <- sig_workflow %>%
        filter(Gene == gene_symbol)

    wf_data <- unique(dat_gene$Workflow)

    if (nrow(wf_row) >= 1) {
        wf_sig <- wf_row$Workflow[1]

        if (wf_sig %in% wf_data) {
            chosen_wf <- wf_sig
            src <- "sig_workflow"
        } else {
            # If the annotated workflow has no expression rows (mismatch),
            # fall back to workflow with most data for this gene
            chosen_wf <- dat_gene %>%
                count(Workflow, sort = TRUE) %>%
                slice_head(n = 1) %>%
                pull(Workflow)
            src <- "fallback_most_data_no_match"
        }
    } else {
        # No sig_workflow entry: fallback based purely on data
        chosen_wf <- dat_gene %>%
            count(Workflow, sort = TRUE) %>%
            slice_head(n = 1) %>%
            pull(Workflow)
        src <- "fallback_most_data_no_sig"
    }

    # Subset to chosen workflow
    dat <- dat_gene %>%
        filter(Workflow == chosen_wf)

    # Apply outlier removal only for Olink (where that specific rule came from)
    if (chosen_wf == "Olink_discovery_cohort" &&
        all(c("Intensity", "M_value") %in% names(dat))) {
        dat <- dat %>%
            filter(!(Intensity < -2 & M_value < 30))
    }

    if (nrow(dat) == 0) {
        warning("No usable data for ", gene_symbol, " in workflow ", chosen_wf)
        return(NULL)
    }

    message(gene_symbol, ": using ", chosen_wf,
            " [", src, "], n = ", nrow(dat))

    ggplot(dat, aes(x = M_value, y = Intensity)) +
        geom_point(size = 2, alpha = 0.3, color = "darkred") +
        geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "black") +
        labs(
            title = paste0(gene_symbol, " (", chosen_wf, ", n=", nrow(dat), ")"),
            x = expression(paste("Clamp M-value (", mu, "mol·kg"^-1,"·min"^-1,")")),
            y = "Protein abundance"
        ) +
        theme_minimal(base_size = 14) +
        theme(
            panel.grid = element_blank(),
            axis.text  = element_text(color = "black"),
            axis.title = element_text(color = "black"),
            plot.title = element_text(face = "bold"),
            axis.line  = element_line(color = "black", linewidth = 0.4),
            legend.position = "none"
        )
}


sig_workflow %>% filter(Gene == "LPL")
expr_long %>% filter(Gene == "LPL") %>% count(Workflow)




# Generate the three panels
p1 <- plot_gene("CES1")
p2 <- plot_gene("INS-CPEPTIDE")
p3 <- plot_gene("CDH15")
p4 <- plot_gene("EPCAM")
p5 <- plot_gene("LPL")
p6 <- plot_gene("IL18R1")

#intestine
#plot_gene("LGALS4")

#plot_gene("PRAP1")

#
#plot_gene("IL18R1")

#heart
#plot_gene("ENG")

#ADIPOSe tissue


# Stack vertically
library(patchwork)

combined_plot <- (p1 + p2) /
    (p3 + p4) /
    (p5 + p6) +
    plot_annotation(title = "Selected insulin sensitivity–associated proteins")

combined_plot







view(tissue_long)

## =========================================
## 0) Setup
## =========================================
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(igraph)
library(ggraph)
library(tibble)

## Parameters
min_rel_weight       <- 0.25   # min relative nTPM for tissue assignment
max_genes_per_tissue <- 25     # cap proteins per tissue
exclude_tissues      <- c("placenta")
angle_deg            <- 305    # rotation of final layout
min_gap_tissue       <- 0.5    # min dist protein–tissue
min_gap_prot         <- 0.2    # min dist protein–protein


## =========================================
## 1) Build tissue_long from HPA / nTPM
## =========================================

# Significant proteins (by UniProt)
sig_uni <- All_unique_significant %>%
    transmute(Uniprot = toupper(str_trim(UniProt_ID))) %>%
    distinct()

sig_df <- merged_data %>%
    mutate(Uniprot = toupper(str_trim(Uniprot))) %>%
    semi_join(sig_uni, by = "Uniprot") %>%
    dplyr::select(any_of(c("GeneID","Uniprot","RNA.tissue.specific.nTPM")))

tissue_long <- sig_df %>%
    filter(!is.na(RNA.tissue.specific.nTPM)) %>%
    separate_rows(RNA.tissue.specific.nTPM, sep = ";") %>%
    separate(RNA.tissue.specific.nTPM, into = c("tissue","nTPM"), sep = ": ") %>%
    mutate(
        tissue = str_trim(tolower(tissue)),
        nTPM   = suppressWarnings(as.numeric(nTPM))
    ) %>%
    filter(!is.na(nTPM)) %>%
    group_by(Uniprot) %>%
    mutate(rel_weight = nTPM / max(nTPM, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(rel_weight >= min_rel_weight) %>%
    filter(!tissue %in% tolower(exclude_tissues)) %>%
    group_by(tissue) %>%
    arrange(desc(nTPM), .by_group = TRUE) %>%
    slice_head(n = max_genes_per_tissue) %>%
    ungroup()

## =========================================
## 2) Build edges + nodes
## =========================================

edges <- tissue_long %>%
    transmute(
        from   = Uniprot,              # protein (UniProt)
        to     = str_to_title(tissue), # tissue name
        weight = rel_weight
    )

nodes_tissue <- edges %>%
    distinct(to) %>%
    dplyr::rename(node = to) %>%
    mutate(type = "tissue") %>%
    left_join(
        edges %>% group_by(to) %>% summarise(size = sum(weight), .groups = "drop"),
        by = c("node" = "to")
    )

nodes_prot <- edges %>%
    distinct(from) %>%
    dplyr::rename(node = from) %>%
    mutate(type = "protein", size = 1)

nodes <- bind_rows(nodes_tissue, nodes_prot) %>%
    distinct(node, .keep_all = TRUE)

## =========================================
## 3) Attach M-value effect sizes to proteins
## =========================================

est_df <- All_unique_significant %>%
    transmute(
        Uniprot  = toupper(str_trim(UniProt_ID)),
        pre_p    = suppressWarnings(as.numeric(`M-value_p.value_pre`)),
        post_p   = suppressWarnings(as.numeric(`M-value_p.value_post`)),
        est_pre  = suppressWarnings(as.numeric(`M-value_estimate_pre`)),
        est_post = suppressWarnings(as.numeric(`M-value_estimate_post`))
    ) %>%
    filter(!is.na(Uniprot), Uniprot != "4", str_detect(Uniprot, "[A-Z]")) %>%
    mutate(
        sig_pre  = !is.na(pre_p)  & pre_p  < 0.05,
        sig_post = !is.na(post_p) & post_p < 0.05,
        est = case_when(
            sig_pre  & !sig_post ~ est_pre,
            sig_pre  &  sig_post ~ est_pre,
            !sig_pre &  sig_post ~ est_post,
            TRUE ~ NA_real_
        )
    ) %>%
    filter(!is.na(est)) %>%
    group_by(Uniprot) %>%
    summarise(est = median(est, na.rm = TRUE), .groups = "drop") %>%
    mutate(Uniprot = toupper(str_trim(Uniprot)))

edges_f <- edges %>%
    mutate(from = toupper(from)) %>%
    left_join(est_df, by = c("from" = "Uniprot"))

cap <- quantile(abs(edges_f$est), 0.95, na.rm = TRUE)
edges_f <- edges_f %>%
    mutate(est_capped = pmax(pmin(est, cap), -cap))

## Restrict to tissues with significant enrichment (from combined_tissue_results)

tissues_to_exclude2 <- c("testis","placenta","vagina","prostate","fallopian tube",
                         "seminal vesicle","ovary","endometrium 1","epididymis","cervix")

sig_tissues <- combined_tissue_results %>%
    filter(padj < 0.05) %>%
    mutate(tissue = str_trim(tolower(tissue))) %>%
    filter(!tissue %in% tolower(tissues_to_exclude2)) %>%
    distinct(tissue) %>%
    pull(tissue)

edges_f <- edges_f %>%
    mutate(tissue_l = str_trim(tolower(to))) %>%
    filter(tissue_l %in% sig_tissues) %>%
    dplyr::select(from, to, weight, est_capped) %>%
    distinct()

nodes_f <- nodes %>%
    filter(node %in% edges_f$from | node %in% edges_f$to) %>%
    distinct(node, .keep_all = TRUE)

## =========================================
## 4) Graph + base layout (stress) + spacing
## =========================================

g_f <- graph_from_data_frame(
    d = edges_f %>% dplyr::select(from, to, weight, est_capped),
    vertices = nodes_f %>% dplyr::select(node, type, size),
    directed = FALSE
)

V(g_f)$node_class <- ifelse(V(g_f)$type == "tissue", "tissue", "protein")

rotate_layout <- function(layout_df, angle_deg = 0) {
    angle <- angle_deg * pi / 180
    layout_df %>%
        mutate(
            x_new =  x * cos(angle) - y * sin(angle),
            y_new =  x * sin(angle) + y * cos(angle)
        ) %>%
        mutate(x = x_new, y = y_new) %>%
        dplyr::select(-x_new, -y_new)
}

set.seed(42)
layout_adj <- ggraph::create_layout(g_f, layout = "stress")

# ensure node_class present
if (!"node_class" %in% names(layout_adj)) {
    layout_adj$node_class <- ifelse(layout_adj$type == "tissue", "tissue", "protein")
}

# rotate
layout_adj <- rotate_layout(layout_adj, angle_deg = angle_deg)




# push proteins away from tissues
tissue_pos <- layout_adj %>%
    filter(node_class == "tissue") %>%
    dplyr::select(name, x_t = x, y_t = y)

for (i in seq_len(nrow(layout_adj))) {
    if (layout_adj$node_class[i] == "protein") {
        dx <- layout_adj$x[i] - tissue_pos$x_t
        dy <- layout_adj$y[i] - tissue_pos$y_t
        d2 <- dx^2 + dy^2
        j  <- which.min(d2)
        dist <- sqrt(d2[j])
        if (!is.na(dist) && dist < min_gap_tissue) {
            vx <- dx[j]; vy <- dy[j]
            if (vx == 0 && vy == 0) {
                vx <- 1e-3; vy <- 0; dist <- sqrt(vx^2 + vy^2)
            }
            scale <- min_gap_tissue / dist
            layout_adj$x[i] <- tissue_pos$x_t[j] + vx * scale
            layout_adj$y[i] <- tissue_pos$y_t[j] + vy * scale
        }
    }
}

# repel protein–protein
prot_idx <- which(layout_adj$node_class == "protein")
max_iter <- 50

for (iter in seq_len(max_iter)) {
    moved <- FALSE
    px <- layout_adj$x[prot_idx]
    py <- layout_adj$y[prot_idx]
    for (i in seq_along(prot_idx)) {
        for (j in seq((i + 1), length(prot_idx))) {
            dx <- px[j] - px[i]
            dy <- py[j] - py[i]
            dist <- sqrt(dx^2 + dy^2)
            if (!is.na(dist) && dist > 0 && dist < min_gap_prot) {
                push <- (min_gap_prot - dist) / 2
                ux <- dx / dist
                uy <- dy / dist
                px[i] <- px[i] - ux * push
                py[i] <- py[i] - uy * push
                px[j] <- px[j] + ux * push
                py[j] <- py[j] + uy * push
                moved <- TRUE
            }
            if (!is.na(dist) && dist == 0) {
                px[i] <- px[i] - min_gap_prot / 2
                px[j] <- px[j] + min_gap_prot / 2
                moved <- TRUE
            }
        }
    }
    layout_adj$x[prot_idx] <- px
    layout_adj$y[prot_idx] <- py
    if (!moved) break
}















## =========================================
## 5A) PLOT WITHOUT LABELS
## =========================================
p_no_labels <- ggraph(layout_adj) +
    geom_edge_link(
        aes(width = weight, colour = est_capped),
        alpha = 0.65, lineend = "round", show.legend = TRUE
    ) +
    scale_edge_width(range = c(0.12, 1.2),
                     name  = "Tissue specificity\n(nTPM-weighted)") +
    scale_edge_colour_gradient2(
        low = "#d73027", mid = "white", high = "#4575b4",
        midpoint = 0, name = "M-value effect\n(estimate)"
    ) +
    geom_node_point(
        aes(shape = node_class, size = node_class),
        stroke = 0.35, fill = "white", show.legend = TRUE
    ) +
    scale_shape_manual(
        name   = "Node type",
        values = c(tissue = 22, protein = 21),
        labels = c("Tissue", "Protein")
    ) +
    scale_size_manual(
        name   = "Node type",
        values = c(tissue = 6, protein = 2),
        labels = c("Tissue", "Protein")
    ) +
    # tissue labels only
    geom_node_text(
        data = layout_adj %>% dplyr::filter(node_class == "tissue"),
        aes(x = x, y = y, label = name),
        fontface = "bold",
        size = 3,
        vjust = -0.8
    ) +
    guides(
        edge_width  = guide_legend(order = 1, override.aes = list(colour = "grey50")),
        edge_colour = ggraph::guide_edge_colourbar(order = 2),
        shape       = guide_legend(order = 3),
        size        = guide_legend(order = 3)
    ) +
    coord_equal() +
    theme_void() +
    labs(title = NULL)

print(p_no_labels)


## =========================================
## 5B) PLOT WITH GENE LABELS (top proteins)
## =========================================

# UniProt -> Gene (adjust columns if needed)
uniprot_gene_map <- All_unique_significant %>%
    transmute(
        Uniprot = toupper(str_trim(UniProt_ID)),
        Gene    = coalesce(str_trim(ID), str_trim(ID))
    ) %>%
    distinct() %>%
    filter(!is.na(Uniprot), !is.na(Gene), Gene != "")

# top proteins per tissue by |est_capped|
top_gene_labels <- edges_f %>%
    filter(!is.na(est_capped)) %>%
    group_by(to) %>%
    slice_max(order_by = abs(est_capped), n = 12, with_ties = FALSE) %>%
    ungroup() %>%
    pull(from) %>%
    unique() %>%
    toupper()

layout_labeled <- layout_adj %>%
    mutate(
        Uniprot = toupper(name),
        Gene = uniprot_gene_map$Gene[match(Uniprot, uniprot_gene_map$Uniprot)],
        label = case_when(
            node_class == "tissue" ~ name,
            node_class == "protein" &
                !is.na(Gene) &
                Uniprot %in% top_gene_labels ~ Gene,
            TRUE ~ ""
        ),
        is_tissue = node_class == "tissue"
    )

p_with_labels <- ggraph(layout_labeled) +
    geom_edge_link(
        aes(width = weight, colour = est_capped),
        alpha = 0.65, lineend = "round", show.legend = TRUE
    ) +
    scale_edge_width(range = c(0.12, 1.8),
                     name  = "Tissue specificity\n(nTPM-weighted)") +
    scale_edge_colour_gradient2(
        low = "#d73027", mid = "white", high = "#4575b4",
        midpoint = 0, name = "M-value effect\n(estimate)"
    ) +
    geom_node_point(
        aes(shape = node_class, size = node_class),
        stroke = 0.35, fill = "white", show.legend = TRUE
    ) +
    scale_shape_manual(
        name   = "Node type",
        values = c(tissue = 22, protein = 21),
        labels = c("Tissue", "Protein")
    ) +
    scale_size_manual(
        name   = "Node type",
        values = c(tissue = 6, protein = 3),
        labels = c("Tissue", "Protein")
    ) +
    geom_node_text(
        data = layout_labeled %>% filter(label != ""),
        aes(x = x, y = y, label = label,
            fontface = ifelse(is_tissue, "bold", "plain")),
        size = 2.4,
        vjust = -0.6,
        check_overlap = TRUE
    ) +
    guides(
        edge_width  = guide_legend(order = 1, override.aes = list(colour = "grey50")),
        edge_colour = ggraph::guide_edge_colourbar(order = 2),
        shape       = guide_legend(order = 3),
        size        = guide_legend(order = 3)
    ) +
    coord_equal() +
    theme_void() +
    labs(title = NULL)

print(p_with_labels)





