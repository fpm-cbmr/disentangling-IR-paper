library(data.table)
library(ggplot2)

# Read results
res_anypre_anypost <- fread("data-raw/Inter99_cohort/Inter99_HR_CMD_extended.csv")

res_anypre_anypost <- res_anypre_anypost[, -c(2,3)]
res_anypre_anypost <- res_anypre_anypost[res_anypre_anypost$disease != "N_Diabetes mellitus", ]
res_anypre_anypost <- res_anypre_anypost[res_anypre_anypost$nevent > 19, ]

library(data.table)

fwrite(
  res_anypre_anypost,
  "res_anypre_anypost.csv",
  sep = ";",
  dec = ","
)

plot_forest_HR_clean <- function(res_df, title = "",
                                 dodge_width = 0.75,
                                 hr_min = 0.5,
                                 hr_max = 2.5) {
  library(data.table)
  library(ggplot2)

  dt <- as.data.table(res_df)

  diseases_keep <- c(
    "Chronic kidney disease", "Type 2 diabetes", "Death", "Obesity",
    "Chronic obstructive pulmonary disease [COPD]",
    "Malignant neoplasm of the of bronchus and lung",
    "Sleep apnea", "N_Stroke", "Cerebral infarction [Ischemic stroke]",
    "Hypertension"
  )

  dt <- dt[
    disease %in% diseases_keep &
      nevent > 19 &
      .id %in% c("M_Value", "Matsuda_Index", "CRS", "HOMAIR", "TG_HDL_ratio", "HbA1c")
  ]

  if (nrow(dt) == 0L) return(ggplot() + theme_void())

  dt[, HR_plot := `exp(coef)`]
  dt[, HR_lo := `lower .95`]
  dt[, HR_hi := `upper .95`]

  dt <- dt[!is.na(HR_plot) & !is.na(HR_lo) & !is.na(HR_hi)]

  pred_labels <- c(
    "M_Value"       = "Predicted M-value",
    "Matsuda_Index" = "Matsuda index",
    "CRS"           = "CRS",
    "HOMAIR"        = "HOMA-IR",
    "TG_HDL_ratio"  = "TG/HDL ratio",
    "HbA1c"         = "HbA1c"
  )

  pred_order <- c(
    "M_Value",
    "Matsuda_Index",
    "CRS",
    "HOMAIR",
    "TG_HDL_ratio",
    "HbA1c"
  )

  dt[, predictor_lbl := factor(
    pred_labels[.id],
    levels = rev(pred_labels[pred_order])
  )]

  dt <- dt[!is.na(predictor_lbl)]

  ## Use M-value HR to order diseases
  comp <- dt[.id == "M_Value", .(disease, HR_plot, nevent)]
  ord <- comp[order(HR_plot)]$disease

  cases_map <- dt[, .(nevent = unique(nevent)[1L]), by = disease]
  cases_map[, disease_n := paste0(disease, " (", format(nevent, big.mark = ","), ")")]

  setkey(cases_map, disease)
  dt[cases_map, disease_n := i.disease_n, on = "disease"]

  ord_n <- cases_map[match(ord, disease), disease_n]
  ord_n <- rev(unique(ord_n[!is.na(ord_n)]))

  dt[, outcome := factor(disease_n, levels = ord_n)]

  ## Significance labels
  dt[, sig_label := ""]
  dt[p.sign == "p<0.05",   sig_label := "*"]
  dt[p.sign == "p<0.01",   sig_label := "**"]
  dt[p.sign == "p<0.001",  sig_label := "***"]
  dt[p.sign == "p<0.0001", sig_label := "****"]

  cols <- c(
    "Predicted M-value" = "darkblue",
    "Matsuda index"     = "#D55E00",
    "CRS"               = "#CC79A7",
    "HOMA-IR"           = "#009E73",
    "TG/HDL ratio"      = "#8DA0CB",
    "HbA1c"             = "#E69F00"
  )

  shapes <- c(
    "Predicted M-value" = 15,
    "Matsuda index"     = 16,
    "CRS"               = 17,
    "HOMA-IR"           = 18,
    "TG/HDL ratio"      = 8,
    "HbA1c"             = 4
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
          ymin = hr_min,
          ymax = hr_max),
      inherit.aes = FALSE,
      fill = "grey96",
      color = NA
    ) +
      geom_hline(
          yintercept = seq(hr_min + 0.25, hr_max - 0.25, by = 0.25),
          linewidth = 0.25,
          color = "grey88"
      ) +
      geom_hline(
          yintercept = 1,
          linetype = "longdash",
          linewidth = 0.45,
          color = "grey40"
      ) +
    geom_errorbar(
      data = dt,
      aes(
        x = outcome,
        y = HR_plot,
        ymin = HR_lo,
        ymax = HR_hi,
        color = predictor_lbl,
        group = predictor_lbl
      ),
      width = 0.025,
      linewidth = 0.45,
      position = pd,
      alpha = 0.75,
      na.rm = TRUE
    ) +
    geom_point(
      data = dt,
      aes(
        x = outcome,
        y = HR_plot,
        color = predictor_lbl,
        shape = predictor_lbl,
        group = predictor_lbl
      ),
      position = pd,
      size = 2.2,
      stroke = 0.6,
      na.rm = TRUE
    ) +
    geom_text(
      data = dt,
      aes(
        x = outcome,
        y = HR_hi,
        label = sig_label,
        group = predictor_lbl
      ),
      position = pd,
      hjust = -0.25,
      size = 2.5,
      color = "black",
      na.rm = TRUE
    ) +
      scale_color_manual(values = cols, name = "Predictor") +
    scale_shape_manual(values = shapes, name = "Predictor") +
    scale_y_continuous(
      limits = c(hr_min, hr_max),
      breaks = seq(hr_min, hr_max, by = 0.25),
      expand = expansion(mult = c(0.01, 0.08))
    ) +
    coord_flip(clip = "off") +
    labs(
      x = NULL,
      y = "Hazard ratio (95% CI)",
      title = title
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y  = element_text(size = 8.8, color = "black"),
      axis.text.x  = element_text(size = 10, color = "black"),
      axis.title.x = element_text(color = "black"),
      axis.line.x  = element_line(color = "black"),
      axis.line.y  = element_line(color = "black"),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10),
      plot.title = element_text(face = "bold", hjust = 0, color = "black"),
      plot.margin = margin(5.5, 35, 5.5, 5.5)
    ) +
    guides(
      color = guide_legend(reverse = TRUE),
      shape = guide_legend(reverse = TRUE)
    )
}

p_hr_anypre_anypost <- plot_forest_HR_clean(
  res_anypre_anypost,
  title = "Inter99",
  dodge_width = 0.75,
  hr_min = 0,
  hr_max = 2.25
)

#8x12
p_hr_anypre_anypost

