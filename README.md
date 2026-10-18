# Lassa Fever Hybrid ARIMA-ANN Forecast — R Code

This repository contains the R analysis code for:

> Nyamekye KA, et al. (2025). *Application of Data Science for Zoonotic Disease Modelling in West Africa: A Hybrid ARIMA-ANN Approach, A Retrospective Time-Series Analysis.* Health Science Reports. Manuscript ID: 1662025.

---

## Repository Structure

| Script | Description |
|---|---|
| `01_data_preprocessing.R` | Loads monthly Lassa fever case counts (NCDC, 2012–2022), runs stationarity tests (ADF, KPSS), and splits data into training (2012–2019) and test (2020–2022) sets |
| `02_arima_model.R` | Fits a seasonal ARIMA model via AIC-based automatic order selection; generates in-sample fitted values and out-of-sample forecasts; extracts residuals for the ANN stage |
| `03_ann_model.R` | Trains a feed-forward ANN (nnetar) on ARIMA residuals to capture nonlinear patterns the linear component cannot model |
| `04_hybrid_combination.R` | Combines ARIMA and ANN outputs into the hybrid prediction; evaluates RMSE, MAE, MAPE, and sMAPE; performs Diebold-Mariano tests with HAC Newey-West standard errors and Bonferroni correction; implements rolling-origin cross-validation |

---

## Data Source

Monthly confirmed Lassa fever case counts for Nigeria (January 2012 – December 2022) were obtained from the Nigeria Centre for Disease Control (NCDC) weekly situation reports, available at: https://ncdc.gov.ng/diseases/sitreps/?cat=5&name=An%20Update%20of%20Lassa%20fever%20in%20Nigeria

---

## R Packages Required

```r
install.packages(c("forecast", "tseries", "ggplot2", "dplyr",
                   "lubridate", "sandwich", "lmtest", "car"))
```

R version used: 4.3.x or later recommended.

---

## How to Run

Run the scripts in order:

```r
source("01_data_preprocessing.R")
source("02_arima_model.R")
source("03_ann_model.R")
source("04_hybrid_combination.R")
```

Each script saves intermediate `.RData` files that the next script loads. Plots are saved to a `plots/` subdirectory (created automatically).

---

## Contact

Corresponding author: Dr. Kofi Afriyie Nyamekye  
University of Professional Studies, Accra (UPSA)  
For data or code queries, contact the corresponding author via the journal.
