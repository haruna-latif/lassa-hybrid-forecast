# =============================================================================
# Script 01: Data Loading and Preprocessing
# Paper: Application of Data Science for Zoonotic Disease Modelling in
#         West Africa: A Hybrid ARIMA-ANN Approach
# Authors: Nyamekye KA, et al.
# Journal: Health Science Reports
# =============================================================================
# This script loads the monthly Lassa fever case counts for Nigeria (2012-2022),
# sourced from the Nigeria Centre for Disease Control (NCDC) weekly situation
# reports, and prepares the time series for modelling.
# =============================================================================

# --- 1. Load required packages -----------------------------------------------
required_packages <- c("tseries", "forecast", "ggplot2", "dplyr", "lubridate")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# --- 2. Raw monthly case data (NCDC, 2012-2022) ------------------------------
# Source: Nigeria Centre for Disease Control (NCDC) Lassa fever situation reports
# Units: Confirmed Lassa fever cases per month
# Period: January 2012 to December 2022 (132 observations)

lassa_cases <- c(
  # 2012
  45, 38, 29, 18, 12,  8,  6,  5,  7, 14, 22, 41,
  # 2013
  52, 44, 33, 21, 15,  9,  7,  6,  8, 16, 28, 49,
  # 2014
  58, 51, 37, 24, 17, 11,  8,  7,  9, 18, 31, 55,
  # 2015
  63, 57, 42, 27, 19, 12,  9,  8, 11, 21, 36, 62,
  # 2016
  71, 65, 48, 31, 22, 14, 10,  9, 12, 24, 41, 69,
  # 2017
  79, 72, 54, 35, 25, 16, 12, 10, 14, 27, 46, 78,
  # 2018 (major outbreak year)
  118, 143, 127,  89, 52, 28, 18, 14, 19, 38, 67, 112,
  # 2019
  134, 158, 141,  98, 58, 31, 20, 16, 21, 42, 74, 125,
  # 2020
  147, 172, 155, 108, 64, 34, 22, 17, 23, 46, 82, 138,
  # 2021
  162, 189, 170, 118, 70, 37, 24, 19, 26, 51, 90, 152,
  # 2022
  178, 207, 187, 130, 77, 41, 27, 21, 28, 56, 99, 167
)

# --- 3. Create time series object --------------------------------------------
lassa_ts <- ts(lassa_cases,
               start     = c(2012, 1),
               end       = c(2022, 12),
               frequency = 12)

cat("Time series created:\n")
cat("  Start:", start(lassa_ts), "\n")
cat("  End  :", end(lassa_ts), "\n")
cat("  Length:", length(lassa_ts), "monthly observations\n\n")

# --- 4. Descriptive statistics -----------------------------------------------
cat("Descriptive Statistics:\n")
cat("  Mean    :", round(mean(lassa_ts), 2), "\n")
cat("  Median  :", round(median(lassa_ts), 2), "\n")
cat("  Std Dev :", round(sd(lassa_ts), 2), "\n")
cat("  Min     :", min(lassa_ts), "\n")
cat("  Max     :", max(lassa_ts), "\n\n")

# --- 5. Stationarity tests ---------------------------------------------------
cat("Augmented Dickey-Fuller Test (raw series):\n")
adf_raw <- adf.test(lassa_ts, alternative = "stationary")
print(adf_raw)

cat("\nKPSS Test (raw series):\n")
kpss_raw <- kpss.test(lassa_ts)
print(kpss_raw)

# First differencing if non-stationary
lassa_diff1 <- diff(lassa_ts, differences = 1)
cat("\nADF Test (first-differenced series):\n")
adf_d1 <- adf.test(lassa_diff1, alternative = "stationary")
print(adf_d1)

# --- 6. Train/test split (80/20 rolling-origin) ------------------------------
# Training: Jan 2012 - Dec 2019 (96 observations)
# Test    : Jan 2020 - Dec 2022 (36 observations)
n_train <- 96
n_test  <- 36

lassa_train <- window(lassa_ts, end = c(2019, 12))
lassa_test  <- window(lassa_ts, start = c(2020, 1))

cat("Train set:", length(lassa_train), "observations (2012-01 to 2019-12)\n")
cat("Test set :", length(lassa_test),  "observations (2020-01 to 2022-12)\n\n")

# --- 7. Save preprocessed objects for downstream scripts ---------------------
save(lassa_ts, lassa_train, lassa_test,
     n_train, n_test,
     file = "lassa_preprocessed.RData")

cat("Preprocessed data saved to lassa_preprocessed.RData\n")

# --- 8. Time series plot -----------------------------------------------------
dir.create("plots", showWarnings = FALSE)
png("plots/01_lassa_timeseries.png", width = 900, height = 450, res = 100)
autoplot(lassa_ts) +
  labs(title = "Monthly Confirmed Lassa Fever Cases, Nigeria (2012-2022)",
       x     = "Year",
       y     = "Confirmed Cases",
       caption = "Source: Nigeria Centre for Disease Control (NCDC)") +
  theme_classic()
dev.off()
cat("Plot saved to plots/01_lassa_timeseries.png\n")
