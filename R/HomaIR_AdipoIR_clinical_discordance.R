#### Muscle, Adipose, and liver discordance ONLY clinical plots. See version 2 for proteome plots ####
library(dplyr)

df_ir <- Discovery_clinical_data %>%
  dplyr::mutate(
    AdipoIR_log2 = log2(AdipoIR)
  )

fit_M <- lm(
  M..µmol..kg.min.. ~ AdipoIR_log2 + Gender + Age + BMI..kg.m2.,
  data = df_ir,
  na.action = na.exclude
)

df_ir$M_discordance <- resid(fit_M)

fit_Adipo <- lm(
  AdipoIR_log2 ~ M..µmol..kg.min.. + Gender + Age + BMI..kg.m2.,
  data = df_ir,
  na.action = na.exclude
)

df_ir$AdipoIR_discordance <- resid(fit_Adipo)

Discovery_clinical_data$M_discordance <- df_ir$M_discordance
Discovery_clinical_data$AdipoIR_discordance <- df_ir$AdipoIR_discordance


library(dplyr)
library(ggplot2)

df_ir <- Discovery_clinical_data %>%
  mutate(
    AdipoIR_log2 = log2(AdipoIR),
    HomaIR_log2 = log2(HOMA1.IR)
  )

fit_M <- lm(
  M..µmol..kg.min.. ~ AdipoIR_log2 + Gender + Age + BMI..kg.m2.,
  data = df_ir,
  na.action = na.exclude
)

fit_Adipo <- lm(
  AdipoIR_log2 ~ M..µmol..kg.min.. + Gender + Age + BMI..kg.m2.,
  data = df_ir,
  na.action = na.exclude
)

fit_HOMA <- lm(
    HomaIR_log2 ~ M..µmol..kg.min.. + Gender + Age + BMI..kg.m2.,
    data = df_ir,
    na.action = na.exclude
)

df_ir <- df_ir %>%
  mutate(
    M_discordance = resid(fit_M),
    AdipoIR_discordance = resid(fit_Adipo),
    HomaIR_discordance = resid(fit_HOMA),
    M_fitted = fitted(fit_M),
    AdipoIR_fitted = fitted(fit_Adipo),
    HomaIR_fitted = fitted(fit_HOMA)
  )


ggplot(df_ir, aes(x = AdipoIR_log2, y = M..µmol..kg.min..)) +
  geom_segment(
    aes(
      x = AdipoIR_fitted,
      xend = AdipoIR_log2,
      y = M..µmol..kg.min..,
      yend = M..µmol..kg.min..,
      color = AdipoIR_discordance
    ),
    linewidth = 0.1,
    alpha = 0.7
  ) +
  geom_point(aes(color = AdipoIR_discordance), size = 3) +
  scale_color_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "AdipoIR (log2)",
    y = "M-value (µmol/kg/min)",
    color = "AdipoIR discordance"
  )

df_ir <- df_ir[Discovery_sample_data$Clamp == "Pre",]

df_ir_ranked <- df_ir %>%
    dplyr::arrange(AdipoIR_discordance) %>%
    dplyr::mutate(rank = dplyr::row_number())


ggplot(df_ir_ranked, aes(x = rank, y = AdipoIR_discordance)) +
    geom_col(aes(fill = AdipoIR_discordance), width = 1) +

    scale_fill_gradient2(
        low = "#2166ac",
        mid = "white",
        high = "#b2182b",
        midpoint = 0
    ) +

    theme_classic(base_size = 14) +
    theme(
        axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black")
    ) +

    labs(
        x = "Individuals (ranked by discordance)",
        y = "AdipoIR discordance",
        fill = "Discordance"
    )

df_ir_ranked_M <- df_ir %>%
    dplyr::arrange(M_discordance) %>%
    dplyr::mutate(rank = dplyr::row_number())




df_ir_ranked_HOMA <- df_ir %>%
    dplyr::arrange(HomaIR_discordance) %>%
    dplyr::mutate(rank = dplyr::row_number())


ggplot(df_ir_ranked_HOMA, aes(x = rank, y = HomaIR_discordance)) +
    geom_col(aes(fill = HomaIR_discordance), width = 1) +

    scale_fill_gradient2(
        low = "#2166ac",
        mid = "white",
        high = "#b2182b",
        midpoint = 0
    ) +

    theme_classic(base_size = 14) +
    theme(
        axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black")
    ) +

    labs(
        x = "Individuals",
        y = "log2(HOMA-IR) discordance",
        fill = "Discordance"
    )



ggplot(df_ir, aes(x = AdipoIR_log2, y = M..µmol..kg.min..)) +
    geom_segment(
        aes(
            x = AdipoIR_fitted,
            xend = AdipoIR_log2,
            y = M..µmol..kg.min..,
            yend = M..µmol..kg.min..,
            color = AdipoIR_discordance
        ),
        linewidth = 0.1,
        alpha = 0.7
    ) +
    geom_point(aes(color = AdipoIR_discordance), size = 3) +
    scale_color_gradient2(
        low = "#2166ac",
        mid = "#F2F2F2",
        high = "#b2182b",
        midpoint = 0
    ) +
    theme_classic(base_size = 14) +
    labs(
        x = "AdipoIR (log2)",
        y = "M-value (µmol/kg/min)",
        color = "AdipoIR discordance"
    )




ggplot(df_ir, aes(x = HomaIR_log2, y = M..µmol..kg.min..)) +
    geom_segment(
        aes(
            x = HomaIR_fitted,
            xend = HomaIR_log2,
            y = M..µmol..kg.min..,
            yend = M..µmol..kg.min..,
            color = HomaIR_discordance
        ),
        linewidth = 0.1,
        alpha = 0.7
    ) +
    geom_point(aes(color = HomaIR_discordance), size = 3) +
    scale_color_gradient2(
        low = "#2166ac",
        mid = "#F2F2F2",
        high = "#b2182b",
        midpoint = 0
    ) +
    theme_classic(base_size = 14) +
    labs(
        x = "HOMA-IR (log2)",
        y = "M-value (µmol/kg/min)",
        color = "HOMA-IR discordance"
    )
