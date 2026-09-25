library(data.table)
library(ggplot2)

res_anypre_anypost <- fread("data-raw/results_anypre_anypost_10y_wBCAA_LDL_FINAL.csv")

res_anypre_anypost$disease == c("Obesity","Dementias","Type 2 diabetes","Chronic kidney disease","Sleep apnea","Heart failure",
                                "Cholangitis","Myocardial infarction [Heart attack]","Malignant neoplasm of the bladder","Chronic obstructive pulmonary disease [COPD]",
                                "Angina pectoris","Malignant neoplasm of bronchus/lung","Cerebral infarction [Ischemic stroke]","Coronary atherosclerosis [Atherosclerotic heart disease]",
                                "Malignant neoplasm of the esophagus", "Osteoarhitis","Gallstones [Cholelithiasis]","Spondylosis",
                                "Colorectal cancer","Malignant neoplasm of the breast","Urinary incontinence and enuresis",
                                "Malignant neoplasm of other and ill-defined sites", "Non-Hodgin lymphoma","Kidney stone disease")


c("Obesity","Dementias","Type 2 diabetes","Chronic kidney disease","Sleep apnea","Heart failure",
  "Cholangitis","Myocardial infarction [Heart attack]","Malignant neoplasm of the bladder","Chronic obstructive pulmonary disease [COPD]",
  "Angina pectoris","Malignant neoplasm of bronchus/lung","Cerebral infarction [Ischemic stroke]","Coronary atherosclerosis",
  "Malignant neoplasm of the esophagus", "Osteoarthritis","Gallstones","Spondylosis",
  "Colorectal cancer","Malignant neoplasm of the breast","Urinary incontinence and enuresis",
  "Malignant neoplasm of other/ill-defined sites", "Non-Hodgkin lymphoma","Kidney stone disease") %in% res_anypre_anypost$disease


res_anypre_anypost$disease == c("Non-Hodgkin lymphoma","Kidney stone disease","Urinary incontinence and enuresis")

res_anypre_anypost$predictor == "LDLadj"

plot_forest_allCI_nowrap_clean <- function(res_df, title = "",
                                           dodge_width = 0.85,
                                           cindex_min = 0.25) {
    library(data.table)
    library(ggplot2)

    dt <- as.data.table(res_df)

    diseases_keep <- c("Obesity",
                       "Dementias",
                       "Fatty liver disease (FLD)",
                       "End stage renal disease [CKD stage 5]",
                       "Type 2 diabetes",
                       "Chronic kidney disease",
                       "Sleep apnea",
                       "Heart failure",
                       "Cholangitis",
                       "Hypertension",
                       "Myocardial infarction [Heart attack]",
                       "Malignant neoplasm of the bladder",
                       "Chronic obstructive pulmonary disease [COPD]",
                       "Angina pectoris",
                       "Malignant neoplasm of bronchus/lung",
                       "Cerebral infarction [Ischemic stroke]",
                       "Coronary atherosclerosis",
                       "Malignant neoplasm of the esophagus",
                       "Osteoarthritis",
                       "Gallstones",
                       "Spondylosis",
                       "Colorectal cancer",
                       "Malignant neoplasm of other/ill-defined sites")

    drop_codes <- c("CV_406.11")

    dt <- dt[
        cases > 0 &
            !phecode %in% drop_codes &
            disease %in% diseases_keep &
            !predictor %in% c("LDLadj", "Leucine", "Isoleucine", "Valine")
    ]

    comp <- dt[predictor == "Pseudo_M_Value", .(outcome_label, C_test)]
    if (nrow(comp) == 0L) return(ggplot() + theme_void())

    ord <- comp[order(C_test)][["outcome_label"]]

    clean_base <- function(x) {
        x <- gsub("\\bcases\\s*=\\s*[0-9][0-9,\\.]*", "", x, ignore.case = TRUE)
        x <- gsub("\\bcontrols?\\s*=\\s*[0-9][0-9,\\.]*", "", x, ignore.case = TRUE)
        x <- gsub("\\(\\s*\\)", "", x)
        x <- gsub("\\s*[,;:/-]+\\s*$", "", x)
        x <- gsub("\\s{2,}", " ", x)
        trimws(x)
    }

    cases_map <- dt[, .(cases = unique(cases)[1L]), by = outcome_label]
    cases_map[, base := clean_base(outcome_label)]
    cases_map[, outcome_label_n := paste0(base, " (", format(cases, big.mark = ","), ")")]

    setkey(cases_map, outcome_label)
    dt[cases_map, outcome_label_n := i.outcome_label_n, on = "outcome_label"]

    ord_n <- cases_map[match(ord, outcome_label), outcome_label_n]
    ord_n <- ord_n[!is.na(ord_n)]

    dt <- dt[!is.na(C_test)]
    dt[, outcome := factor(outcome_label_n, levels = ord_n)]

    pred_labels <- c(
        "Pseudo_M_Value" = "Predicted M-value",
        "HbA1c"          = "HbA1c",
        "TG_HDL_ratio"   = "TG/HDL ratio",
        "LDL_0"          = "LDL"
    )

    pred_order <- c(
        "Pseudo_M_Value",
        "HbA1c",
        "TG_HDL_ratio",
        "LDL_0"
    )

    dt[, predictor_lbl := factor(
        pred_labels[predictor],
        levels = rev(pred_labels[pred_order])
    )]

    dt <- dt[!is.na(predictor_lbl)]

    cols <- c(
        "Predicted M-value" = "darkblue",
        "HbA1c"             = "#E69F00",
        "TG/HDL ratio"      = "#8DA0CB",
        "LDL"               = "#66A61E"
    )

    shapes <- c(
        "Predicted M-value" = 15,
        "HbA1c"             = 16,
        "TG/HDL ratio"      = 17,
        "LDL"               = 18
    )

    ltys <- c(
        "Predicted M-value" = "solid",
        "HbA1c"             = "dashed",
        "TG/HDL ratio"      = "dotted",
        "LDL"               = "dotdash"
    )

    bg <- data.table(
        outcome = factor(ord_n, levels = ord_n),
        outcome_i = seq_along(ord_n)
    )[outcome_i %% 2 == 1L]

    pd <- position_dodge(width = dodge_width)

    ggplot() +
        geom_rect(
            data = bg,
            aes(xmin = outcome_i - 0.5,
                xmax = outcome_i + 0.5,
                ymin = cindex_min,
                ymax = 1),
            inherit.aes = FALSE,
            fill = "grey96",
            color = NA
        ) +
        geom_errorbar(
            data = dt,
            aes(x = outcome, y = C_test,
                ymin = C_boot_lo, ymax = C_boot_hi,
                color = predictor_lbl,
                linetype = predictor_lbl),
            width = 0.025,
            linewidth = 0.45,
            position = pd,
            alpha = 0.75,
            na.rm = TRUE
        ) +
        geom_point(
            data = dt,
            aes(x = outcome, y = C_test,
                color = predictor_lbl,
                shape = predictor_lbl),
            position = pd,
            size = 2.2,
            stroke = 0.6,
            na.rm = TRUE
        ) +
        geom_hline(
            yintercept = 0.5,
            linetype = "longdash",
            linewidth = 0.4,
            color = "grey40"
        ) +
        scale_color_manual(values = cols, name = "Predictor") +
        scale_shape_manual(values = shapes, name = "Predictor") +
        scale_linetype_manual(values = ltys, name = "Predictor") +
        scale_y_continuous(
            limits = c(cindex_min, 1),
            breaks = seq(0.3, 1, by = 0.1),
            expand = expansion(mult = c(0.01, 0.03))
        ) +
        coord_flip(clip = "off") +
        labs(
            x = NULL,
            y = "C-index (test set; 95% CI)",
            title = title
        ) +
        theme_minimal(base_size = 12) +
        theme(
            axis.text.y  = element_text(size = 8.8, color = "black", lineheight = 0.95, margin = margin(r = 6)),
            axis.text.x  = element_text(size = 10, color = "black"),
            axis.title.x = element_text(color = "black"),
            axis.title.y = element_text(color = "black"),
            axis.line.x  = element_line(color = "black"),
            axis.line.y  = element_line(color = "black"),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            legend.position = "right",
            legend.text = element_text(size = 9),
            legend.title = element_text(size = 10),
            plot.title = element_text(face = "bold", hjust = 0.0, color = "black")
        ) +
        guides(
            color = guide_legend(
                reverse = TRUE,
                override.aes = list(size = 3)
            ),
            shape = guide_legend(
                reverse = TRUE,
                override.aes = list(size = 3)
            ),
            linetype = guide_legend(
                reverse = TRUE,
                override.aes = list(linewidth = 0.8)
            )
        )
}


p_anypre_anypost_nowrap <- plot_forest_allCI_nowrap_clean(
    res_anypre_anypost,
    "Prevalence = ≥1 pre code excluded; Incident = ≥1 post code (first event); 10y follow-up",
    dodge_width = 0.85,
    cindex_min = 0.3
)

p_anypre_anypost_nowrap
#export 9x10



#### Suppl. figure ####
plot_forest_Mvalue_BCAA <- function(res_df, title = "",
                                    dodge_width = 0.85,
                                    cindex_min = 0.25) {
    library(data.table)
    library(ggplot2)

    dt <- as.data.table(res_df)

    diseases_keep <- c("Obesity",
                       "Dementias",
                       "Fatty liver disease (FLD)",
                       "End stage renal disease [CKD stage 5]",
                       "Type 2 diabetes",
                       "Chronic kidney disease",
                       "Sleep apnea",
                       "Heart failure",
                       "Cholangitis",
                       "Hypertension",
                       "Myocardial infarction [Heart attack]",
                       "Malignant neoplasm of the bladder",
                       "Chronic obstructive pulmonary disease [COPD]",
                       "Angina pectoris",
                       "Malignant neoplasm of bronchus/lung",
                       "Cerebral infarction [Ischemic stroke]",
                       "Coronary atherosclerosis",
                       "Malignant neoplasm of the esophagus",
                       "Osteoarthritis",
                       "Gallstones",
                       "Spondylosis",
                       "Colorectal cancer",
                       "Malignant neoplasm of other/ill-defined sites")

    drop_codes <- c("CV_406.11")

    predictors_keep <- c(
        "Pseudo_M_Value",
        "Leucine",
        "Isoleucine",
        "Valine"
    )

    dt <- dt[
        cases > 0 &
            !phecode %in% drop_codes &
            disease %in% diseases_keep &
            predictor %in% predictors_keep
    ]

    comp <- dt[predictor == "Pseudo_M_Value", .(outcome_label, C_test)]
    if (nrow(comp) == 0L) return(ggplot() + theme_void())

    ord <- comp[order(C_test)][["outcome_label"]]

    clean_base <- function(x) {
        x <- gsub("\\bcases\\s*=\\s*[0-9][0-9,\\.]*", "", x, ignore.case = TRUE)
        x <- gsub("\\bcontrols?\\s*=\\s*[0-9][0-9,\\.]*", "", x, ignore.case = TRUE)
        x <- gsub("\\(\\s*\\)", "", x)
        x <- gsub("\\s*[,;:/-]+\\s*$", "", x)
        x <- gsub("\\s{2,}", " ", x)
        trimws(x)
    }

    cases_map <- dt[, .(cases = unique(cases)[1L]), by = outcome_label]
    cases_map[, base := clean_base(outcome_label)]
    cases_map[, outcome_label_n := paste0(base, " (", format(cases, big.mark = ","), ")")]

    setkey(cases_map, outcome_label)
    dt[cases_map, outcome_label_n := i.outcome_label_n, on = "outcome_label"]

    ord_n <- cases_map[match(ord, outcome_label), outcome_label_n]
    ord_n <- ord_n[!is.na(ord_n)]

    dt <- dt[!is.na(C_test)]
    dt[, outcome := factor(outcome_label_n, levels = ord_n)]

    pred_labels <- c(
        "Pseudo_M_Value" = "Predicted M-value",
        "Leucine"        = "Leucine",
        "Isoleucine"     = "Isoleucine",
        "Valine"         = "Valine"
    )

    pred_order <- c(
        "Pseudo_M_Value",
        "Leucine",
        "Isoleucine",
        "Valine"
    )

    dt[, predictor_lbl := factor(
        pred_labels[predictor],
        levels = rev(pred_labels[pred_order])
    )]

    dt <- dt[!is.na(predictor_lbl)]

    cols <- c(
        "Predicted M-value" = "darkblue",
        "Leucine"           = "#D95F02",
        "Isoleucine"        = "#7570B3",
        "Valine"            = "#E7298A"
    )

    shapes <- c(
        "Predicted M-value" = 15,
        "Leucine"           = 0,
        "Isoleucine"        = 1,
        "Valine"            = 2
    )

    ltys <- c(
        "Predicted M-value" = "solid",
        "Leucine"           = "solid",
        "Isoleucine"        = "dashed",
        "Valine"            = "dotted"
    )

    bg <- data.table(
        outcome = factor(ord_n, levels = ord_n),
        outcome_i = seq_along(ord_n)
    )[outcome_i %% 2 == 1L]

    pd <- position_dodge(width = dodge_width)

    ggplot() +
        geom_rect(
            data = bg,
            aes(xmin = outcome_i - 0.5,
                xmax = outcome_i + 0.5,
                ymin = cindex_min,
                ymax = 1),
            inherit.aes = FALSE,
            fill = "grey96",
            color = NA
        ) +
        geom_errorbar(
            data = dt,
            aes(x = outcome, y = C_test,
                ymin = C_boot_lo, ymax = C_boot_hi,
                color = predictor_lbl,
                linetype = predictor_lbl),
            width = 0.025,
            linewidth = 0.45,
            position = pd,
            alpha = 0.75,
            na.rm = TRUE
        ) +
        geom_point(
            data = dt,
            aes(x = outcome, y = C_test,
                color = predictor_lbl,
                shape = predictor_lbl),
            position = pd,
            size = 2.2,
            stroke = 0.6,
            na.rm = TRUE
        ) +
        geom_hline(
            yintercept = 0.5,
            linetype = "longdash",
            linewidth = 0.4,
            color = "grey40"
        ) +
        scale_color_manual(values = cols, name = "Predictor") +
        scale_shape_manual(values = shapes, name = "Predictor") +
        scale_linetype_manual(values = ltys, name = "Predictor") +
        scale_y_continuous(
            limits = c(cindex_min, 1),
            breaks = seq(0.3, 1, by = 0.1),
            expand = expansion(mult = c(0.01, 0.03))
        ) +
        coord_flip(clip = "off") +
        labs(
            x = NULL,
            y = "C-index (test set; 95% CI)",
            title = title
        ) +
        theme_minimal(base_size = 12) +
        theme(
            axis.text.y  = element_text(size = 8.8, color = "black", lineheight = 0.95, margin = margin(r = 6)),
            axis.text.x  = element_text(size = 10, color = "black"),
            axis.title.x = element_text(color = "black"),
            axis.title.y = element_text(color = "black"),
            axis.line.x  = element_line(color = "black"),
            axis.line.y  = element_line(color = "black"),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            legend.position = "right",
            legend.text = element_text(size = 9),
            legend.title = element_text(size = 10),
            plot.title = element_text(face = "bold", hjust = 0.0, color = "black")
        ) +
        guides(
            color = guide_legend(
                reverse = TRUE,
                override.aes = list(size = 3)
            ),
            shape = guide_legend(
                reverse = TRUE,
                override.aes = list(size = 3)
            ),
            linetype = guide_legend(
                reverse = TRUE,
                override.aes = list(linewidth = 0.8)
            )
        )
}


p_Mvalue_BCAA <- plot_forest_Mvalue_BCAA(
    res_anypre_anypost,
    "Predicted M-value compared with branched-chain amino acids; 10y follow-up",
    dodge_width = 0.85,
    cindex_min = 0.3
)

p_Mvalue_BCAA
#export in 9x10
