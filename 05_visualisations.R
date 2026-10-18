# =============================================================================
# Script 05: All Manuscript Figures (Figures 1-8)
# Paper: Application of Data Science for Zoonotic Disease Modelling in
#         West Africa: A Hybrid ARIMA-ANN Approach
# Authors: Nyamekye KA, et al.
# Journal: Health Science Reports
# =============================================================================
# Generates all eight figures that appear in the manuscript.
# Run after sourcing scripts 01-04 in order.
#
# Objects used from saved .RData files:
#   lassa_preprocessed.RData : lassa_ts, lassa_train, lassa_test, n_train, n_test
#   arima_outputs.RData      : arima_fit, arima_pred, arima_residuals
#   ann_outputs.RData        : ann_fit, ann_residual_pred
#   hybrid_outputs.RData     : hybrid_pred, metrics_table, comparison
# =============================================================================

library(forecast)
library(ggplot2)
library(dplyr)
library(patchwork)

dir.create("plots", showWarnings = FALSE)

# Load all model outputs
load("lassa_preprocessed.RData")
load("arima_outputs.RData")
load("ann_outputs.RData")
load("hybrid_outputs.RData")

actual_test <- as.numeric(lassa_test)
test_time   <- as.numeric(time(lassa_test))
full_time   <- as.numeric(time(lassa_ts))

# =============================================================================
# FIGURE 1: Cumulative confirmed Lassa fever cases over the study period
# =============================================================================
cat("Generating Figure 1: Cumulative cases over time...\n")

fig1_data <- data.frame(
  year    = full_time,
  cumcase = cumsum(as.numeric(lassa_ts))
)

fig1 <- ggplot(fig1_data, aes(x = year, y = cumcase)) +
  geom_line(colour = "#D62728", linewidth = 1.1) +
  labs(
    title   = "Figure 1: Cumulative Confirmed Lassa Fever Cases Over Time\n(Nigeria, 2012–2022)",
    x       = "Year",
    y       = "Cumulative Cases",
    caption = "Source: Nigeria Centre for Disease Control (NCDC) weekly situation reports; authors' computation using R 4.3.1."
  ) +
  theme_classic() +
  theme(plot.caption = element_text(hjust = 0, size = 8))

ggsave("plots/Figure1_cumulative_cases.png", fig1,
       width = 10, height = 5, dpi = 150)
cat("  Saved: plots/Figure1_cumulative_cases.png\n")

# =============================================================================
# FIGURE 2: Monthly confirmed case counts — full time series
# =============================================================================
cat("Generating Figure 2: Monthly case counts (full series)...\n")

fig2_data <- data.frame(
  year  = full_time,
  cases = as.numeric(lassa_ts)
)

fig2 <- ggplot(fig2_data, aes(x = year, y = cases)) +
  geom_line(colour = "#1F77B4", linewidth = 0.7) +
  geom_vline(xintercept = as.numeric(time(lassa_train))[n_train],
             linetype = "dashed", colour = "grey40") +
  annotate("text", x = as.numeric(time(lassa_train))[n_train] - 0.3,
           y = max(as.numeric(lassa_ts)) * 0.9,
           label = "Train | Test", hjust = 1, size = 3, colour = "grey30") +
  labs(
    title   = "Figure 2: Monthly Confirmed Lassa Fever Case Counts\n(Nigeria, 2012–2022)",
    x       = "Year",
    y       = "Confirmed Cases",
    caption = "Source: NCDC weekly situation reports; authors' computation using R 4.3.1 (ggplot2)."
  ) +
  theme_classic() +
  theme(plot.caption = element_text(hjust = 0, size = 8))

ggsave("plots/Figure2_monthly_cases.png", fig2,
       width = 10, height = 5, dpi = 150)
cat("  Saved: plots/Figure2_monthly_cases.png\n")

# =============================================================================
# FIGURE 3: STL seasonal decomposition
# =============================================================================
cat("Generating Figure 3: STL decomposition...\n")

stl_fit <- stl(lassa_ts, s.window = "periodic", robust = TRUE)

png("plots/Figure3_STL_decomposition.png",
    width = 900, height = 700, res = 120)
plot(stl_fit,
     main = "Figure 3: STL Seasonal Decomposition of Monthly Lassa Fever Cases\n(Nigeria, 2012–2022)")
mtext("Source: NCDC weekly situation reports; authors' computation using R 4.3.1 (forecast).",
      side = 1, line = 4, cex = 0.7, adj = 0)
dev.off()
cat("  Saved: plots/Figure3_STL_decomposition.png\n")

# =============================================================================
# FIGURE 4: ARIMA residual diagnostics
# =============================================================================
cat("Generating Figure 4: ARIMA residual diagnostics...\n")

png("plots/Figure4_ARIMA_residual_diagnostics.png",
    width = 900, height = 700, res = 120)
checkresiduals(arima_fit,
               main = paste0("Figure 4: ARIMA Residual Diagnostics\nBest-Fitting Model: ",
                             "ARIMA(", paste(arimaorder(arima_fit), collapse = ","), ")\n",
                             "Source: NCDC weekly situation reports; authors' computation using R 4.3.1 (forecast)."))
dev.off()
cat("  Saved: plots/Figure4_ARIMA_residual_diagnostics.png\n")

# =============================================================================
# FIGURE 5: Hybrid forecast vs ARIMA vs actual — test period
# =============================================================================
cat("Generating Figure 5: Hybrid forecast vs actual (test period)...\n")

fig5_data <- data.frame(
  year   = test_time,
  Actual = actual_test,
  ARIMA  = arima_pred,
  Hybrid = hybrid_pred
)

fig5 <- ggplot(fig5_data, aes(x = year)) +
  geom_line(aes(y = Actual, colour = "Actual"),          linewidth = 1.1) +
  geom_line(aes(y = ARIMA,  colour = "ARIMA"),           linewidth = 0.9, linetype = "dashed") +
  geom_line(aes(y = Hybrid, colour = "Hybrid ARIMA-ANN"), linewidth = 0.9) +
  scale_colour_manual(
    values = c("Actual" = "black", "ARIMA" = "steelblue",
               "Hybrid ARIMA-ANN" = "firebrick"),
    name = NULL
  ) +
  labs(
    title   = "Figure 5: Hybrid ARIMA-ANN Out-of-Sample Forecast versus Actual\nConfirmed Lassa Fever Cases (Test Period: 2020–2022)",
    x       = "Year",
    y       = "Confirmed Cases",
    caption = "Source: NCDC weekly situation reports; authors' computation using R 4.3.1."
  ) +
  theme_classic() +
  theme(legend.position = "bottom",
        plot.caption    = element_text(hjust = 0, size = 8))

ggsave("plots/Figure5_hybrid_forecast_test.png", fig5,
       width = 10, height = 5, dpi = 150)
cat("  Saved: plots/Figure5_hybrid_forecast_test.png\n")

# =============================================================================
# FIGURE 6: 12-month prospective forecast (refitted on full series)
# =============================================================================
cat("Generating Figure 6: 12-month prospective forecast...\n")

arima_full  <- auto.arima(lassa_ts, seasonal = TRUE, stepwise = TRUE,
                           approximation = TRUE, ic = "aic")
fcast_arima <- forecast(arima_full, h = 12)
arima_12    <- as.numeric(fcast_arima$mean)
arima_lo    <- as.numeric(fcast_arima$lower[, 2])   # 95% lower
arima_hi    <- as.numeric(fcast_arima$upper[, 2])   # 95% upper

res_full <- residuals(arima_full)
set.seed(2024)
ann_full <- nnetar(res_full, p = 12, P = 1, size = 7, repeats = 20)
ann_12   <- as.numeric(forecast(ann_full, h = 12)$mean)

hybrid_12 <- arima_12 + ann_12

proj_time <- seq(from   = as.numeric(time(lassa_ts))[length(lassa_ts)] + 1/12,
                 by     = 1/12,
                 length.out = 12)

fig6_hist <- data.frame(
  year  = full_time,
  cases = as.numeric(lassa_ts)
)

fig6_proj <- data.frame(
  year      = proj_time,
  hybrid    = hybrid_12,
  arima_lo  = arima_lo + ann_12,
  arima_hi  = arima_hi + ann_12
)

fig6 <- ggplot() +
  geom_line(data = fig6_hist, aes(x = year, y = cases),
            colour = "grey60", linewidth = 0.7) +
  geom_ribbon(data = fig6_proj,
              aes(x = year, ymin = arima_lo, ymax = arima_hi),
              fill = "firebrick", alpha = 0.15) +
  geom_line(data = fig6_proj, aes(x = year, y = hybrid),
            colour = "firebrick", linewidth = 1.1) +
  geom_point(data = fig6_proj, aes(x = year, y = hybrid),
             colour = "firebrick", size = 2) +
  labs(
    title   = "Figure 6: 12-Month Prospective Hybrid ARIMA-ANN Forecast\n(January 2023 – December 2023; Shaded Band = 95% CI)",
    x       = "Year",
    y       = "Projected Confirmed Cases",
    caption = "Source: NCDC weekly situation reports; authors' computation using R 4.3.1."
  ) +
  theme_classic() +
  theme(plot.caption = element_text(hjust = 0, size = 8))

ggsave("plots/Figure6_prospective_forecast.png", fig6,
       width = 10, height = 5, dpi = 150)
cat("  Saved: plots/Figure6_prospective_forecast.png\n")

# =============================================================================
# FIGURE 7: Bar chart — RMSE comparison across models
# =============================================================================
cat("Generating Figure 7: RMSE bar chart...\n")

# Drop rows where RMSE is NA (ANN standalone not directly evaluated)
fig7_data <- metrics_table[!is.na(metrics_table$RMSE), ]

fig7 <- ggplot(fig7_data, aes(x = Unit, y = RMSE, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7),
           width = 0.6) +
  scale_fill_manual(
    values = c("ARIMA" = "steelblue", "Hybrid ARIMA-ANN" = "firebrick"),
    name = "Model"
  ) +
  labs(
    title   = "Figure 7: Root Mean Square Error (RMSE) — ARIMA vs Hybrid ARIMA-ANN\n(Test Period: 2020–2022)",
    x       = NULL,
    y       = "RMSE (confirmed cases)",
    caption = "Source: Authors' computation using R 4.3.1. Lower RMSE = better forecast accuracy."
  ) +
  theme_classic() +
  theme(legend.position = "bottom",
        plot.caption    = element_text(hjust = 0, size = 8))

ggsave("plots/Figure7_RMSE_bar_chart.png", fig7,
       width = 8, height = 5, dpi = 150)
cat("  Saved: plots/Figure7_RMSE_bar_chart.png\n")

# =============================================================================
# FIGURE 8: Bar chart — MAPE comparison across models
# =============================================================================
cat("Generating Figure 8: MAPE bar chart...\n")

fig8_data <- metrics_table[!is.na(metrics_table$MAPE), ]

fig8 <- ggplot(fig8_data, aes(x = Unit, y = MAPE, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7),
           width = 0.6) +
  scale_fill_manual(
    values = c("ARIMA" = "steelblue", "Hybrid ARIMA-ANN" = "firebrick"),
    name = "Model"
  ) +
  labs(
    title   = "Figure 8: Mean Absolute Percentage Error (MAPE, %) — ARIMA vs Hybrid ARIMA-ANN\n(Test Period: 2020–2022)",
    x       = NULL,
    y       = "MAPE (%)",
    caption = "Source: Authors' computation using R 4.3.1. Lower MAPE = better percentage accuracy."
  ) +
  theme_classic() +
  theme(legend.position = "bottom",
        plot.caption    = element_text(hjust = 0, size = 8))

ggsave("plots/Figure8_MAPE_bar_chart.png", fig8,
       width = 8, height = 5, dpi = 150)
cat("  Saved: plots/Figure8_MAPE_bar_chart.png\n")

cat("\n=== All 8 figures saved to plots/ ===\n")
cat("Figure 1: plots/Figure1_cumulative_cases.png\n")
cat("Figure 2: plots/Figure2_monthly_cases.png\n")
cat("Figure 3: plots/Figure3_STL_decomposition.png\n")
cat("Figure 4: plots/Figure4_ARIMA_residual_diagnostics.png\n")
cat("Figure 5: plots/Figure5_hybrid_forecast_test.png\n")
cat("Figure 6: plots/Figure6_prospective_forecast.png\n")
cat("Figure 7: plots/Figure7_RMSE_bar_chart.png\n")
cat("Figure 8: plots/Figure8_MAPE_bar_chart.png\n")
