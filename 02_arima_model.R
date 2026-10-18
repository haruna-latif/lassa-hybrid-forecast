# =============================================================================
# Script 02: ARIMA Model Fitting and Diagnostics
# Paper: Application of Data Science for Zoonotic Disease Modelling in
#         West Africa: A Hybrid ARIMA-ANN Approach
# Authors: Nyamekye KA, et al.
# Journal: Health Science Reports
# =============================================================================
# Fits a seasonal ARIMA model to the Lassa fever training data using
# auto.arima() with AIC-based order selection. Produces residuals for the
# ANN component (Script 03).
# =============================================================================

library(forecast)
library(tseries)
library(ggplot2)

# --- 1. Load preprocessed data -----------------------------------------------
load("lassa_preprocessed.RData")

# --- 2. ACF / PACF plots for order identification ----------------------------
dir.create("plots", showWarnings = FALSE)

png("plots/02a_acf_pacf.png", width = 900, height = 500, res = 100)
par(mfrow = c(1, 2))
acf(lassa_train,  main = "ACF — Training Series",  lag.max = 36)
pacf(lassa_train, main = "PACF — Training Series", lag.max = 36)
dev.off()
cat("ACF/PACF plot saved.\n")

# --- 3. Fit ARIMA via AIC minimisation ---------------------------------------
cat("Fitting ARIMA model (auto.arima with stepwise search)...\n")
arima_fit <- auto.arima(
  lassa_train,
  seasonal       = TRUE,
  stepwise       = FALSE,   # exhaustive search
  approximation  = FALSE,
  ic             = "aic",
  trace          = TRUE
)

cat("\nSelected ARIMA order:\n")
print(arima_fit)
cat("\nAIC:", AIC(arima_fit), "\n")
cat("BIC:", BIC(arima_fit), "\n")

# --- 4. Residual diagnostics -------------------------------------------------
arima_residuals <- residuals(arima_fit)

cat("\nLjung-Box test on ARIMA residuals (lag 12):\n")
lb_test <- Box.test(arima_residuals, lag = 12, type = "Ljung-Box", fitdf = length(arima_fit$coef))
print(lb_test)

cat("\nShapiro-Wilk normality test on residuals:\n")
sw_test <- shapiro.test(as.numeric(arima_residuals))
print(sw_test)

png("plots/02b_arima_residuals.png", width = 900, height = 700, res = 100)
checkresiduals(arima_fit)
dev.off()
cat("Residual diagnostic plot saved.\n")

# --- 5. In-sample (training) fitted values -----------------------------------
arima_fitted <- fitted(arima_fit)

# --- 6. Out-of-sample ARIMA forecast (test period) ---------------------------
cat("\nGenerating ARIMA forecasts for test period (36 months)...\n")
arima_forecast <- forecast(arima_fit, h = n_test)
arima_pred     <- as.numeric(arima_forecast$mean)

# Accuracy on test set
arima_accuracy <- accuracy(arima_forecast, lassa_test)
cat("\nARIMA forecast accuracy (test set):\n")
print(arima_accuracy)

# --- 7. Save outputs for hybrid combination ----------------------------------
save(arima_fit,
     arima_residuals,
     arima_fitted,
     arima_pred,
     arima_accuracy,
     file = "arima_outputs.RData")

cat("\nARIMA outputs saved to arima_outputs.RData\n")
cat("Residuals passed to Script 03 (ANN) for nonlinear modelling.\n")

# --- 8. Forecast plot --------------------------------------------------------
png("plots/02c_arima_forecast.png", width = 900, height = 450, res = 100)
autoplot(arima_forecast) +
  autolayer(lassa_test, series = "Actual") +
  labs(title = paste0("ARIMA", paste0("(", arimaorder(arima_fit), ")", collapse = ""),
                      " Forecast vs Actual — Test Period"),
       x = "Year", y = "Confirmed Cases") +
  theme_classic()
dev.off()
cat("Forecast plot saved to plots/02c_arima_forecast.png\n")
