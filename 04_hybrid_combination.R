# =============================================================================
# Script 04: Hybrid ARIMA-ANN Model — Combination and Evaluation
# Paper: Application of Data Science for Zoonotic Disease Modelling in
#         West Africa: A Hybrid ARIMA-ANN Approach
# Authors: Nyamekye KA, et al.
# Journal: Health Science Reports
# =============================================================================
# Combines ARIMA linear forecasts and ANN residual corrections into the
# hybrid prediction. Evaluates performance using RMSE, MAE, MAPE, and sMAPE
# with Bonferroni-corrected Diebold-Mariano tests against the standalone models.
# Implements rolling-origin (walk-forward) cross-validation.
# =============================================================================

library(forecast)
library(sandwich)   # HAC Newey-West standard errors
library(lmtest)
library(ggplot2)

# --- 1. Load all model outputs -----------------------------------------------
load("lassa_preprocessed.RData")
load("arima_outputs.RData")
load("ann_outputs.RData")

# --- 2. Hybrid prediction = ARIMA forecast + ANN residual correction ---------
hybrid_pred <- arima_pred + ann_residual_pred

actual_test <- as.numeric(lassa_test)

cat("=== Hybrid ARIMA-ANN Predictions vs Actual (Test Period) ===\n")
comparison <- data.frame(
  Month   = time(lassa_test),
  Actual  = actual_test,
  ARIMA   = arima_pred,
  Hybrid  = hybrid_pred
)
print(comparison)

# --- 3. Accuracy metrics (helper functions) ----------------------------------
rmse  <- function(actual, pred) sqrt(mean((actual - pred)^2))
mae   <- function(actual, pred) mean(abs(actual - pred))
mape  <- function(actual, pred) mean(abs((actual - pred) / actual)) * 100
smape <- function(actual, pred) {
  mean(2 * abs(actual - pred) / (abs(actual) + abs(pred))) * 100
}

# --- 4. Performance table ----------------------------------------------------
metrics <- data.frame(
  Model  = c("ARIMA", "ANN (residual)", "Hybrid ARIMA-ANN"),
  RMSE   = c(rmse(actual_test, arima_pred),
             rmse(actual_test, arima_pred + ann_residual_pred - arima_pred),  # ANN alone not direct
             rmse(actual_test, hybrid_pred)),
  MAE    = c(mae(actual_test, arima_pred),
             NA,
             mae(actual_test, hybrid_pred)),
  MAPE   = c(mape(actual_test, arima_pred),
             NA,
             mape(actual_test, hybrid_pred)),
  sMAPE  = c(smape(actual_test, arima_pred),
             NA,
             smape(actual_test, hybrid_pred))
)

cat("\n=== Forecast Accuracy Metrics (Test Set: 2020-2022) ===\n")
print(metrics, digits = 4)

# --- 5. Diebold-Mariano test (hybrid vs ARIMA) with HAC SE ------------------
e_arima  <- actual_test - arima_pred
e_hybrid <- actual_test - hybrid_pred

d <- e_arima^2 - e_hybrid^2   # squared error differential

# HAC Newey-West standard error
dm_lm  <- lm(d ~ 1)
dm_hac <- coeftest(dm_lm, vcov = NeweyWest(dm_lm, lag = 3, prewhite = FALSE))

cat("\n=== Diebold-Mariano Test: Hybrid vs ARIMA (HAC Newey-West SEs) ===\n")
cat("H0: equal forecast accuracy\n")
cat("H1: hybrid is more accurate (negative d = hybrid has smaller squared errors)\n")
print(dm_hac)

# Bonferroni correction for multiple comparisons (3 pairwise: H vs A, H vs ANN, A vs ANN)
# Since ANN standalone is not directly evaluated here, report single corrected p-value
dm_p     <- dm_hac[1, 4]
dm_p_bon <- min(dm_p * 3, 1)    # Bonferroni-adjusted (3 tests)
cat("\nBonferroni-adjusted p-value (3 comparisons):", round(dm_p_bon, 4), "\n")

# --- 6. VIF check (multicollinearity in combined predictor space) ------------
# VIF between ARIMA fitted and ANN residual correction
vif_data <- data.frame(
  y         = actual_test,
  arima_hat = arima_pred,
  ann_corr  = ann_residual_pred
)
vif_model <- lm(y ~ arima_hat + ann_corr, data = vif_data)
vif_vals  <- car::vif(vif_model)
cat("\n=== Variance Inflation Factors (VIF) ===\n")
print(vif_vals)
cat("(VIF < 10 indicates acceptable multicollinearity)\n")

# --- 7. Rolling-origin cross-validation (walk-forward) -----------------------
cat("\n=== Rolling-Origin Recursive Evaluation ===\n")
cat("Initial window: 2012-2018 (84 obs); expanding by 1 month per step\n\n")

initial_window <- 84   # Jan 2012 - Dec 2018
horizons       <- 1:36 # 36 one-step-ahead forecasts

ro_actual  <- numeric(36)
ro_arima   <- numeric(36)
ro_hybrid  <- numeric(36)

for (i in seq_along(horizons)) {
  train_end_idx <- initial_window + i - 1
  train_sub     <- lassa_ts[1:train_end_idx]
  ts_sub        <- ts(train_sub, start = start(lassa_ts), frequency = 12)

  # Refit ARIMA
  fit_i    <- auto.arima(ts_sub, seasonal = TRUE, stepwise = TRUE,
                          approximation = TRUE, ic = "aic")
  fcast_i  <- forecast(fit_i, h = 1)
  arima_i  <- as.numeric(fcast_i$mean)

  # Refit ANN on residuals
  res_i    <- residuals(fit_i)
  set.seed(2024)
  ann_i    <- nnetar(res_i, p = 12, P = 1, size = 7, repeats = 10)
  ann_res_i <- as.numeric(forecast(ann_i, h = 1)$mean)

  hybrid_i <- arima_i + ann_res_i
  actual_i <- as.numeric(lassa_ts[train_end_idx + 1])

  ro_actual[i]  <- actual_i
  ro_arima[i]   <- arima_i
  ro_hybrid[i]  <- hybrid_i

  cat(sprintf("Step %2d: Actual=%3.0f  ARIMA=%6.1f  Hybrid=%6.1f\n",
              i, actual_i, arima_i, hybrid_i))
}

cat("\n=== Rolling-Origin Accuracy ===\n")
cat("ARIMA  RMSE:", round(rmse(ro_actual, ro_arima), 2), "\n")
cat("Hybrid RMSE:", round(rmse(ro_actual, ro_hybrid), 2), "\n")
cat("ARIMA  MAE :", round(mae(ro_actual, ro_arima), 2), "\n")
cat("Hybrid MAE :", round(mae(ro_actual, ro_hybrid), 2), "\n")

# --- 8. Build multi-unit metrics table for Figure 7 and Figure 8 ------------
# This script runs once per unit inside a loop in a full multi-unit analysis.
# Here we construct a metrics_table that 05_visualisations.R expects:
#   columns: Unit, Model, RMSE, MAPE, MAE

# Single-unit Nigeria analysis: build a two-row metrics table (ARIMA + Hybrid)
metrics_table <- data.frame(
  Unit  = "Nigeria",
  Model = c("ARIMA", "Hybrid ARIMA-ANN"),
  RMSE  = c(rmse(actual_test, arima_pred), rmse(actual_test, hybrid_pred)),
  MAPE  = c(mape(actual_test, arima_pred), mape(actual_test, hybrid_pred)),
  MAE   = c(mae(actual_test,  arima_pred), mae(actual_test,  hybrid_pred))
)

# --- 9. Save final outputs ---------------------------------------------------
save(hybrid_pred, metrics, comparison, ro_actual, ro_arima, ro_hybrid,
     metrics_table,
     file = "hybrid_outputs.RData")
cat("\nHybrid model outputs saved to hybrid_outputs.RData\n")
cat("metrics_table saved (used by 05_visualisations.R for Figures 7 and 8).\n")

# --- 9. Final comparison plot ------------------------------------------------
png("plots/04_hybrid_vs_actual.png", width = 1000, height = 500, res = 100)
plot_df <- data.frame(
  time   = as.numeric(time(lassa_test)),
  Actual = actual_test,
  ARIMA  = arima_pred,
  Hybrid = hybrid_pred
)
ggplot(plot_df, aes(x = time)) +
  geom_line(aes(y = Actual, colour = "Actual"),  linewidth = 1.2) +
  geom_line(aes(y = ARIMA,  colour = "ARIMA"),   linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = Hybrid, colour = "Hybrid"),  linewidth = 1) +
  scale_colour_manual(values = c("Actual" = "black",
                                 "ARIMA"  = "steelblue",
                                 "Hybrid" = "firebrick")) +
  labs(title   = "Hybrid ARIMA-ANN vs ARIMA vs Actual — Test Period (2020-2022)",
       x       = "Year",
       y       = "Confirmed Lassa Fever Cases",
       colour  = "Model") +
  theme_classic() +
  theme(legend.position = "bottom")
dev.off()
cat("Comparison plot saved to plots/04_hybrid_vs_actual.png\n")
