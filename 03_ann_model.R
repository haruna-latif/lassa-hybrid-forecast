# =============================================================================
# Script 03: Artificial Neural Network (ANN) on ARIMA Residuals
# Paper: Application of Data Science for Zoonotic Disease Modelling in
#         West Africa: A Hybrid ARIMA-ANN Approach
# Authors: Nyamekye KA, et al.
# Journal: Health Science Reports
# =============================================================================
# Models the nonlinear structure remaining in the ARIMA residuals using a
# feed-forward ANN (nnetar). The ANN captures patterns the linear ARIMA
# component cannot, forming the second stage of the hybrid model.
# =============================================================================

library(forecast)
library(ggplot2)

# --- 1. Load ARIMA outputs ---------------------------------------------------
load("lassa_preprocessed.RData")
load("arima_outputs.RData")

# --- 2. Prepare residual series for ANN input --------------------------------
# The ANN is trained on ARIMA residuals (training window only)
residual_train <- arima_residuals

cat("ARIMA residuals summary (training period):\n")
print(summary(as.numeric(residual_train)))
cat("\n")

# --- 3. Fit nnetar ANN to residuals ------------------------------------------
# nnetar: neural network autoregression
#   p  = number of lagged inputs (chosen by AIC/CV internally)
#   P  = seasonal lagged inputs (12-month seasonality)
#   size = number of hidden nodes (set to ceiling((p + P + 1) / 2) per Zhang 2003)
#   repeats = 20 random initialisations (best retained)
set.seed(2024)  # reproducibility

cat("Fitting ANN to ARIMA residuals...\n")
ann_fit <- nnetar(
  residual_train,
  p       = 12,    # 12 lagged inputs (one full seasonal cycle)
  P       = 1,     # one seasonal lag
  size    = 7,     # hidden nodes: ceiling((12+1+1)/2)
  repeats = 20
)

cat("\nANN specification:\n")
print(ann_fit)

# --- 4. ANN fitted values (on training residuals) ----------------------------
ann_fitted_residuals <- fitted(ann_fit)

# Remove leading NAs introduced by lagged inputs
valid_idx <- !is.na(ann_fitted_residuals)
cat("\nValid fitted values:", sum(valid_idx), "of", length(ann_fitted_residuals), "\n")

# --- 5. Forecast ANN residual corrections for test period --------------------
cat("\nGenerating ANN residual forecasts for test period (36 months)...\n")
ann_residual_forecast <- forecast(ann_fit, h = n_test)
ann_residual_pred     <- as.numeric(ann_residual_forecast$mean)

cat("ANN residual forecast summary:\n")
print(summary(ann_residual_pred))

# --- 6. Diagnostics on modelled residuals ------------------------------------
modelled_residuals <- as.numeric(residual_train[valid_idx]) -
                      as.numeric(ann_fitted_residuals[valid_idx])

cat("\nLjung-Box test on ANN-modelled residuals (lag 12):\n")
lb_ann <- Box.test(modelled_residuals, lag = 12, type = "Ljung-Box")
print(lb_ann)

# --- 7. Save ANN outputs -----------------------------------------------------
save(ann_fit,
     ann_fitted_residuals,
     ann_residual_pred,
     file = "ann_outputs.RData")

cat("\nANN outputs saved to ann_outputs.RData\n")
cat("ANN residual predictions passed to Script 04 (Hybrid combination).\n")

# --- 8. Plot: ANN fitted vs actual residuals ---------------------------------
png("plots/03_ann_residuals.png", width = 900, height = 450, res = 100)
plot(residual_train,
     main = "ARIMA Residuals: Actual vs ANN-Fitted",
     ylab = "Residual", xlab = "Time",
     col  = "grey50", lwd = 1.5)
lines(ann_fitted_residuals, col = "steelblue", lwd = 2)
legend("topright",
       legend = c("ARIMA Residuals", "ANN Fitted"),
       col    = c("grey50", "steelblue"),
       lwd    = c(1.5, 2))
dev.off()
cat("Residual plot saved to plots/03_ann_residuals.png\n")
