###############################################################
# This script analysis the cancer incidence series in 
# Europe by classification cancer, age group, and sex
# from 2007 to 2017
###############################################################

library(tidyverse) # Data manipulation and Visualizarion
library(lubridate) # Date handling
library(stringr) # to use str_to_lowcase
library(forecast)  # Time series forecasting analysis
library(tseries)   # Statistical tests
library(ggplot2)   # Advanced plotting
library(gridExtra) # Multiple plots
library(stats)     # Base statistical functions
library(seasonal)  # X-13ARIMA-SEATS decomposition
library(zoo)

directory_path = "/home/jsancheg/git_environment/Europe_Cancer_Incidence/data/processed"
source("/home/jsancheg/git_environment/Europe_Cancer_Incidence/src/models/functions_cancer_incidence.R")

df_cancer_incidence <- read.csv(file.path(directory_path, "europe_cancer_cases_long_format.csv"))

head(df_cancer_incidence)
colnames(df_cancer_incidence) <- colnames(df_cancer_incidence) %>% 
  str_to_lower()  


cat('---- Reading the dataset, adding columns day and month to the dataset ----')

df_cancer_incidence <- df_cancer_incidence %>% 
  mutate(month = 12,
         day = 31) %>%
  relocate(month, day, .before = year)

cat ('--- creating column date_reported, and preparing 
     to be used by prophet ---')

df_cancer_incidence <- df_cancer_incidence %>%
  mutate(date_reported = as.Date(paste(year, month, day, sep = "-")),
         # Prepare the dataset to use prophet:
         # Prophet requires specific column names
         # 'ds' for date and 'y' for target
         ds = date_reported,
         y = cases_by_age) %>%
  relocate(date_reported, ds, .after = age_group)  

df_cancer_incidence %>% colnames()

cat('--- Printing columns of the dataset and their type ---')

colnames(df_cancer_incidence)

str(df_cancer_incidence)

cat('--- Filtering Europe cancer incidences from 2007 - 2017 ---')
df_cancer_incidence_2007_2017 <- df_cancer_incidence %>%
  filter(year >= 2007)

colnames(df_cancer_incidence_2007_2017)
head(df_cancer_incidence_2007_2017[c("month","day","year","date_reported","ds")])

cat('--- Subseting the series to be forecasted by cancer type, sex, and age group ---')

item_id <- df_cancer_incidence_2007_2017 %>% 
  select(c("cancer_classification", "sex_character","age_group")) %>%
  unique()

europe_total_cases_year <- df_cancer_incidence_2007_2017 %>% 
  select(c(year,y)) %>%
  group_by(year) %>%
  summarise(y = sum(y))

europe_total_cases_year

europe_cases_sex_year <- df_cancer_incidence_2007_2017 %>% 
  select(c(year,sex_character,y)) %>%
  group_by(year,sex_character) %>%
  summarise(y = sum(y))


tail(europe_cases_sex_year, 12)

cat ('--- Subsetting the series ---')
europe_df <- df_cancer_incidence_2007_2017 %>% 
  select(ds,y) %>% 
  group_by(ds) %>%
  summarise(y = sum(y))


start_year <- year(min(europe_df$ds))


head(europe_df)

ts_europe <- ts(
  europe_df$y,
  start = start_year,
  frequency = 1
)

cat("Time Series Summary:\n")
cat("Start:", start(ts_europe), "\n")
cat("End:", end(ts_europe), "\n")
cat("Frequency:", frequency(ts_europe), "\n")
cat("Number of observations:", length(ts_europe), "\n\n")

#============================================================================
# 2.-Initial visualization ---------------------------------------------------
#============================================================================


p1 <- plot_basic(ts_europe, "Cancer Incidence in Europe from 2007 to 2017")
print(p1)

#============================================================================
# 3.-Descriptive Statistics -----------------------------------------------
#============================================================================

cat ("==== Descriptive Statistics ===\n")
cat ("Mean: ", mean(ts_europe, na.rm = TRUE), "\n")
cat ("Median: ", median(ts_europe, na.rm = TRUE), "\n")
cat ("SD: ", sd(ts_europe, na.rm = TRUE), "\n")
cat ("Min: ", min(ts_europe, na.rm = TRUE), "\n")
cat ("Max: ", max(ts_europe, na.rm = TRUE), "\n")
cat ("Range: ", max(ts_europe, na.rm = TRUE) - min(ts_europe, na.rm = TRUE),"\n")
cat ("CV (%): ", (sd(ts_europe, na.rm = TRUE)/ mean(ts_europe, na.rm = TRUE)) * 100, "\n\n")



#============================================================================
# 4.-Stationary Tests -----------------------------------------------
#============================================================================

cat("===== Stationary Tests =======")

# 4.1 Augmented Dickey-Fuller Test
cat('---- Augmented Dikey-Fuller Tests --- \n')
adf_test <- adf.test(ts_europe, alternative = "stationary")
cat("ADF Statistic:", adf_test$statistic, "\n")
cat("P-value:", adf_test$p.value, "\n")
cat("Interpretation:",
    ifelse(adf_test$p.value < 0.05,
           "Series is STATIONARY (reject H0)",
           "Series is NON-STATIONARY (fail to reject H0)"), "\n\n")

# 4.2 KPSS Test
cat("--- KPSS Test ---\n")
kpss_test <- kpss.test(ts_europe, null = "Trend")
cat("Kpss Statistics:", kpss_test$statistic, "\n")
cat("P-value:", kpss_test$p.value, "\n")
cat("Interpretation:", 
    ifelse(kpss_test$p.value < 0.05,
           "Series is NON-STATIONARY (reject H0)",
           "Series is STATIONARY (fail to reject H0)"), "\n\n")

# 4.3 Philips-Perron Test
cat("--- Philips-Person Test ---\n")
pp_test <- PP.test(ts_europe)
cat("PP Statistics:", pp_test$statistic, "\n")
cat("P-value:", pp_test$p.value, "\n")
cat("Interpretation:",
    ifelse(pp_test$p.value < 0.05,
           "Series in STATIONARY (reject H0)",
           "Series is NON-STATIONARY (fall to reject H0)"), "\n\n")
# 4.4 Combined Interpretaion
cat("--- Combined Stationary Assessment --- \n")
if (adf_test$p.value < 0.05 && kpss_test$p.value >= 0.05){
  cat("CONCLUSION: Series is likely STATIONARY\n\n")
} else if (adf_test$p.value >= 0.05 && kpss_test$p.value < 0.05) {
  cat("CONCLUSION: Tests disagree - series may be trend stationary\n")
  cat("RECOMMENDATION: Consider detrending rather than differencing \n\n")
} else {
  cat("CONCLUSION: Tests disagree - further investigacion needed \n\n")
}

#============================================================================
# 5.-Autocorrelation ------------------------------------------------------
#============================================================================

cat("=== Autocorrelation Analysis ===\n\n")

# 5.1 ACF and PACF plots
par(mfrow = c(2,1))
acf(ts_europe, main = "Autocorrelation Function (ACF)", lag.max = 20)
pacf(ts_europe, main = "Partial Autocorrelation Function (PACF)", lag.max = 20)
par(mfrow = c(1,1))

# 5.2 Ljung-Box Test for autocorrelation
cat("--- Ljung-Box Test (lag = 10) ---\n")
lb_test <- Box.test(ts_europe, lag = 10, type = "Ljung-Box")
cat("Chi-squared:", lb_test$statistic, "\n")
cat("P-value:", lb_test$p.value, "\n")
cat("Interpretation:", 
    ifelse(lb_test$p.value < 0.05,
           "Significant autocorrelation detected",
           "No significant autocorrelation"), "\n\n")

#============================================================================
# 6.- Decomposition Methods -----------------------------------------------
#============================================================================

cat("=== Time Series Decomposition ===\n\n")

# 6.1 Classical Decomposition (Additive)
# Norte: Requires at least 2 complete cycles (2 years for annual, 
# 24 months for monthly)

if (length(ts_europe) >= 2 * frequency(ts_europe)){
  
  cat("--- Classical Additive Decomposition ---\n")
  decomp_additive <- decompose(ts_europe, type = "additive")
  
  # Plot decomposition
  plot(decomp_additive)
  
  # Extract components
  trend_add <- decomp_additive$trend
  seasonal_add <- decomp_additive$seasonal
  random_add <- decomp_additive$random
  
  # Summary statistics of components
  cat("\nTrend component statistics:\n")
  print(Summary(as.numeric(na.omit(trend_europe))))
  
  cat("\nSeasonal component statisics:\n")
  print(summary(as.numeric(na.omit(seasonal_add))))
  
  cat("\nRandom component statistics:\n")
  print(summary(as.numeric(na.omit(random_add))))
  
  # Strength of trend and seasonality
  var_random <- var(na.omit(random_add))
  var_detrend <- var(na.omit(ts_europe - trend_add))
  var_deseas <- var(na.omit(ts_europe - seasonal_add))


strength_trend <- max(0, 1 - var_random / var_deseas)
strength_seasonal <- max(0, 1 - var_random / var_detrend)

cat("\nStrength of Trent:", round(strength_trend, 3), "\n")
cat("Strength of Seasonality:", round(strength_seasonal, 3), "\n\n")

# 6.2 Classical Decomposition (Multiplicative)
cat("--- Classical Multiplicacion Decomposition ---\n")

# Check for non-positive values
if(all(ts_europe > 0)){
  decomp_mult <- decompose(ts_europe, type = "Multiplicative")
  plot(decomp_mult)
} else {
  cat("Cannot perform multiplicative decomposition: data contains non-positive values\n")
}

} else {
  cat("Insufficient data for statistical decomposition (need at least 2 cycles")
}

# 6.3 STL Decomposition (Seasonal and Trend decomposition using Loess)
# More robust and flexible than classical decomposition

if(frequency(ts_europe) > 1) {
  cat("\n--- STL Decomposition ---\n")
  
  stl_decomp <- stl(ts_europe, s.window = "periodic")
  plot(stl_decomp)
  
  # Extract components
  stl_trend <- stl_decomp$time_series[, "trend"]
  stl_seasonal <- stl_decomp$time_series[, "seasonal"]
  stl_remainder <- stl_decomp$time_series[, "remainder"]
  
  cat("\nSTL Trend component statistics:\n")
  print(summary(as.numeric(stl_trend)))
  
  cat("\nSTL Seasonal component statistics:\n")
  print(summary(as.numeric(stl_seasonal)))
  
  cat("\nSTL Remainder component statistics:\n")
  print(summary(as.numeric(stl_remainder)))
  
}else {
  cat("\nSTL decomposition skipped: no seasonality in annual data\n")
}

# 6.4 Moving Average Decomposition (for trend extraction)
cat("\n--- Moving Average Trend Extraccion ---\n")

# Choose appropriate window (e.g., 3-year moving average for annual data)
ma_order <- min(5, floor(length(ts_europe) / 3 )) # Adaptive order
ma_trend <- ma(ts_europe, order = ma_order, centre = TRUE)

# Plot original vs trend
plot_ma <- autoplot(ts_europe) +
  autolayer(ma_trend, series = "MA Trend", linewidth = 1.2) +
  labs(
    title = paste0("Time Series with ", ma_order, 
                   "-Period Moving Average Trend"),
    x = "Time",
    y = "Cases",
    color = "Series",
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(plot_ma)

# Detrended series
detrended <- ts_europe - ma_trend

cat("\nDetrended series statistics:\n")
print(summary(as.numeric(na.omit(detrended))))

#============================================================================
# 7.- Rolling Statistics --------------------------------------------------
#============================================================================

cat("\n=== Rolling Statistics Analysis ===\n")

# Calculate rolling mean and standard deviation
window_size <- min(5, floor(length(ts_europe)/3))

rolling_mean <- rollmean(ts_europe, k = window_size, 
                         fill = NA, align = "right")

rolling_sd <- rollapply(ts_europe, width = window_size, 
                        FUN = sd, fill = NA,
                        align = "right")

# Create dataframe for plotting
rolling_df <- data.frame(
  time = time(ts_europe),
  original = as.numeric(ts_europe),
  rolling_mean = as.numeric(rolling_mean),
  rolling_sd = as.numeric(rolling_sd)
)


# Plot rolling statistics
p_rolling <- ggplot(rolling_df, aes(x = time)) +
  geom_line(aes(y = original, color = "Original"),
            linewidth = 0.8) +
  geom_line(aes(y = rolling_mean, color = "Rolling Mean"), 
            linewidth = 1.2) +
  labs(
    title = paste0("Time Series with ", 
                   window_size, "-Period Rolling Mean"),
    x = "Time",
    y = "Cases",
    color = "Series"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p_rolling)

# Plot rolling standard deviation
p_rolling_sd <- ggplot(rolling_df, aes(x = time, y = rolling_sd)) +
  geom_line(color = "darkred", linewidth = 1) +
  labs(
    title = paste0(window_size, "-Period Rolling Standard Deviation"),
    x = "Time",
    y = "Standard Deviation"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p_rolling_sd)

# Test for changing variance
cat("\nVariance stability check:\n")
first_half_var <- var(
  as.numeric(window(ts_europe, 
                    end = time(ts_europe)[floor(length(ts_europe)/2)] )))
second_half_var <- var(
  as.numeric(window(ts_europe, 
                    start = time(ts_europe)[ceiling(length(ts_europe)/2)])))

cat("First half variance:", first_half_var, "\n")
cat("Second half variance:", second_half_var, "\n")
cat("Variance ratio:", second_half_var / first_half_var, "\n")

if(abs(log(second_half_var / first_half_var)) > 0.5){
  cat("WARNING: Variance appears to be changing over time \n")
  cat("RECOMMENDATION: Consider log transformation or variance stabilization\n")
}

#============================================================================
# 8. Transformation Analysis ----------------------------------------------
#============================================================================

cat("\n=== Transformation Analysis ===\n")

# 8.1 Box-Cox Transformation
if(all(ts_europe > 0)) {
  cat("\n--- Box-Cox Transformation ---\n")
  
  lambda <- BoxCox.lambda(ts_europe)
  cat("Optimal Lambda:", round(lambda, 2), "\n")
  
  # Interpretation of Lambda
  if(abs(lambda) < 0.1) {
    cat("Interpretation: Log Transformation recommended (lambda = 0)\n")
  }else if (abs(lambda - 05) < 0.1) {
    cat("Interpretation: Square root transformation recommended (lambda = 0.5)\n")
  }else if (abs(lambda - 1) < 0.1) {
    cat("Interpretation: No transformation needed (lambda = 1)\n")
  }else {
    cat("Interpretation: Power transformation with lambda =", 
        round(lambda, 3), "\n")
  }
  

# Apply transformation
  ts_transformed <- BoxCox(ts_europe, lambda)
  
# Plot comparison
  par(mfrow = c(2,1))
  plot(ts_europe, main = "Original Series", ylab = "Cases")
  plot(ts_transformed, main = paste0("Box-Cox Transformed (
                                     λ=", round(lambda,3),")"),
       ylab = "Transformed")
  par(mfrow = c(1,1))


} else {
  cat("Box-Cox transformation skipped: data contains non-positive values\n")
}  

# 8.2 Log Transformation (if aplicable)
if(all(ts_europe > 0)) {
  ts_log <- log(ts_europe)
  
  cat("\n--- Log Transformation ---\n")
  cat("Original CV:", sd(ts_europe) / mean(ts_europe), "\n")
  cat("Log-transformed CV:", sd(ts_log) / mean(ts_log), "\n")  
  
  # Plot
  par(mfrow = c(2,1))
  plot(ts_europe, main = "Original Series", ylab = "Cases")
  plot(ts_log, main = "Log transformed Series", ylab = "Log(Cases)")
  par(mfrow = c(1,1))
}

# 8.3 First Difference
ts_diff1 <- diff(ts_europe, differences = 1)

cat("\n--- First Differencing ---\n")
cat("Original Series - ADF p-value:", adf.test(ts_europe)$p.value, "\n")
cat("First differenced - ADF p-value:", adf.test(ts_diff1)$p.value, "\n")

if(adf.test(ts_diff1)$p.value < 0.05) {
  cat("First differencing achieves stationary\n")
}

# Plot of differenced series
acf(ts_diff1, main = "ACF of First Differenced Series")


serie_df <- filtering_item(df_cancer_incidence,item_id[1,])



# Checking  for missing dates and completeness
cat("Date range: ", min(serie_df$ds), " to ", max(serie_df$ds), "\n")
cat("Number of observations: ", nrow(serie_df), "\n")
cat("Missing values in y: ", sum(is.na(serie_df$y)), "\n")

#============================================================================
# 9.- Trend Analysis ------------------------------------------------------
#============================================================================

cat("\n=== Trend Analysis ===\n")

# 9.1 Linear Trend
time_index <- 1:length(ts_europe)
linear_model <- lm(as.numeric(ts_europe) ~ time_index)

cat("\n--- Linear Trend Model ---\n")
print(summary(linear_model))

# Extract coefficients
intercept <- coef(linear_model)[1]
slope <- coef(linear_model)[2]

# Extract p-values correctly
p_value_intercept <- t(as.matrix(summary(linear_model)$coefficients))[4,1] 
p_value_slope <- t(as.matrix(summary(linear_model)$coefficients))[4,2]

cat("\nTrend Interpretation:\n")
cat("Starting level:", round(intercept, 2), "\n")
cat("Annual change:", round(slope,2), "cases per year\n")
cat("Percentage change per year:", round((slope/intercept)*100, 2), "%\n")
cat("\nP-values:\n")
cat("Intercept p-value:", round(p_value_intercept, 4), "\n")
cat("Slope p-value:", round(p_value_slope,4), "\n\n")

if(p_value_slope < 0.05){
  if(slope > 0){
    cat("Significant increasing trend detected (p < 0.05)\n")
  }else {
    cat("Significant decreasing tren detected (p < 0.05)\n")
  }
} else {
  cat("No significant trend detected (p >= 0.05)\n")
}

# 9.2 Polynomial trend (quadratic)
  time_index_sq <- time_index^2
  poly_model <- lm(as.numeric(ts_europe) ~ time_index + time_index_sq)
  
  cat("\n--- Quadratic Trend Model ---\n")
  print(summary(poly_model))
  
# Compare models
cat("\nModel Comparison:\n")
cat("Linear R²:", summary(linear_model)$r.squared, "\n")
cat("Quadratic R²", summary(poly_model)$r.squared, "\n")

# AIC comparison
cat("Linear AIC:", AIC(linear_model), "\n")
cat("Quadratic AIC:", AIC(poly_model), "\n")

# 9.3 Plot trends
fitted_lineal <- fitted(linear_model)
fitted_poly <- fitted(poly_model)

plot_trends <- ggplot(data.frame(
  time = time(ts_europe),
  original = as.numeric(ts_europe),
  linear = fitted_lineal,
  quadratic = fitted_poly),
  aes(x = time)) + 
  geom_point(aes(y = original, color = "Original"), size = 2) +
  geom_line(aes(y = original, color = "Original"), linewidth = 0.5, alpha = 0.5) +
  geom_line(aes(y = linear, color = "Linear Trend"), linewidth = 1.2) +
  geom_line(aes(y = quadratic, color = "Quadratic Trend"), linewidth = 1.2, linetype = "dashed") +
  labs(
    title = "Trend Analysis: Linear vs Quadratic",
    x = "Time",
    y = "Cases",
    color = "Series"
  ) +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(plot_trends)


#============================================================================
# 10.- Seasonality Analysis (if applicable) -------------------------------
#============================================================================

if (frequency(ts_europe) > 1) {
  cat("\n=== Seasonality Analysis ===\n")
  
  # Seasonal subseries plot
  seamsomplot(ts_europe, year.labels = TRUE, main = "Seasonal Subseries Plot")
  
  # Month plot (for monthly data)
  monthplot(ts_europe, main = "Seasonal Deviation Plot")
  
  # Test for seasonality
  cat("\n--- Test for Seasonality ---\n")
  
  # QS test (if available)
  # This test for the presence of seasonality
} else {
  cat("\n=== No seasonality (Annual Data) ===\n")
}


#============================================================================
# 11.- Outlier Detection --------------------------------------------------
#============================================================================

cat("\n=== Outlier Detection ===\n")

# 11.1 Statistical outliers (using IQR method)
Q1 <- quantile(ts_europe, 0.25)
Q3 <- quantile(ts_europe, 0.75)
IQR <- Q3 - Q1

lower_bound <- Q1 - 1.5 * IQR
upper_bound <- Q3 + 1.5 * IQR

outliers_iqr <- which(ts_europe < lower_bound | ts_europe > upper_bound)

cat("IQR Method:\n")
cat("Lower bound:", lower_bound, "\n")
cat("Upper bound:", upper_bound, "\n")
cat("Number of outliers:", length(outliers_iqr), "\n")

if(length(outliers_iqr) > 0){
  cat("Outlier position:", outliers_iqr, "\n")
  cat("Outlier values:", ts_europe[outliers_iqr], "\n")
}

# 11.2 Z-scored method
z_scores <- scale(ts_europe)
outliers_z <- which(abs(z_scores) > 3)

cat("\nZ-score Method (|z| > 3):\n")
cat("Number of outliers:", length(outliers_z), "\n")

if (length(outliers_z) > 0) {
  cat("Outlier positions:", outliers_z, "\n")
  cat("Outlier values:", ts_europe[outliers_z], "\n")
  cat("Z-scores:", z_scores[outliers_z], "\n")
}

# Plot with outliers highlighted
outlier_df <- data.frame(
  time = time(ts_europe),
  value = as.numeric(ts_europe),
  is_outlier = 1:length(ts_europe) %in% c(outliers_iqr, outliers_z)
)

p_outliers <- ggplot(outlier_df, aes(x = time, y = value)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(aes(color = is_outlier, size = is_outlier)) +
  scale_color_manual(values = c("steelblue", "red"), labels = c("Normal", "Outlier")) +
  scale_size_manual(values = c(2,4), labels = c("Normal", "Outlier")) +
  geom_hline(yintercept = c(lower_bound, upper_bound),
             linetype = "dashed", color = "red", alpha = 0.5) +
  labs(
    title = "Time Series with Outliers Highlighted",
    x = "Time",
    y = "Cases",
    color = "Point Type",
    size = "Point Type"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p_outliers)


#============================================================================
# 12.- Change Point Detection ---------------------------------------------
#============================================================================

cat("\n=== Change Point Detection ===\n")

# Simple CuSum approach
cumsum_vals <- cumsum(ts_europe - mean(ts_europe))

plot(time(ts_europe), cumsum_vals, type = "l",
     main = "CUSUM Chart for Changepoint Detection",
     xlab = "Time", ylab = "CUSUM",
     col = "darkblue", lwd = 2)
abline(h = 0, col = "red", lty = 2)

# Identify potential checkpoints (local maxima/minima)
cumsum_diff <- diff(sign(diff(cumsum_vals)))
potential_changes <- which(abs(cumsum_diff) == 2) + 1

cat("Potential changepoints at positions:", potential_changes, "\n")
if(length(potential_changes) > 0){
  cat("Potential changepoint years:", time(ts_europe)[potential_changes], "\n")
}

#============================================================================
# 13.- Summary Visualizations ---------------------------------------------
#============================================================================

cat("\n=== Generating Summary Visualizations ===\n")

# Create a comprehensive multi-panel plot
par(mfrow = c(3,2))

# Panel 1: Original series
plot(ts_europe, main = "Original Time Series", ylab = "Cases",
     col = "steelblue", lwd = 2)

# Panel 2: ACF
acf(ts_europe, main = "ACF", lag.max = 20)

# Panel 3: PACF
pacf(ts_europe, main = "PACF", lag.max = 20)

# Panel 4: Histogram
hist(ts_europe, main = "Distribution of Values", 
     xlab = "Cases", col = "lightblue", border = "white")

# Panel 5: First difference
plot(ts_europe, main = "First Difference", ylab = "Diff",
     col = "darkgreen", lwd = 2)

# Panel 6: QQ plot
qqnorm(as.numeric(ts_europe), main = "Q-Q Plot")
qqline(as.numeric(ts_europe), col = "red", lwd = 2)

par(mfrow=c(1,1))

#============================================================================
# 14.- Comprehensive Summary Report ---------------------------------------
#============================================================================

cat("\n")
cat("="=rep("=", 70), "\n")
cat("                   Comprehensive Time Series Analysis Summary\n")
cat("="=rep("=", 70), "\n\n")

cat("Data Characteristics:\n")
cat("    - Period:", start(ts_europe), "to", end(ts_europe), "\n")
cat("    - Observations:", length(ts_europe), "\n")
cat("    - Frequency:", frequency(ts_europe), "\n")
cat("    - Mean:", round(mean(ts_europe), 2), "\n")
cat("    - Standard Deviation:", round(sd(ts_europe), 2), "\n")
cat("    - Coefficient of Variation:", round(sd(ts_europe)/mean(ts_europe)*100, 2), "%\n\n")

cat("Trend:\n")
cat("    - Linear trend slope:", round(slope,3), "cases/year\n")
cat("    - Trend significance: p=", round(summary(linear_model)$coefficients[c(7,8)], 4), "\n")
cat("    - Trend is:", ifelse(summary(linear_model)$coefficients[c(7,8)] < 0.05, "SIGNIFICANT", "NO SIGNIFICANT"), "\n\n")

cat("Autocorrelation:\n")
cat("    - Ljung-Box Test: p =", round(lb_test$p.value,4), "\n")
cat("    - Autocorrelation:", ifelse(lb_test$p.value < 0.05, "PRESENT", "NO SIGNIFICANT"), "\n\n")

cat("Outliers:\n")
cat("    - IQR method:", length(outliers_iqr), "outliers_detected\n")
cat("    - Z-score method:", length(outliers_z), "outliers detected\n\n")

cat("Recommendations:\n")

if(adf_test$p.value >= 0.05){
  cat("   1. Series is non-statinary - consider differencing or detrending\n")
}

if(kpss_test$p.value < 0.05 && adf_test$p.value >= 0.05){
  cat("   2. Apply first differencing before modelling\n")
}

if(abs(log(second_half_var / first_half_var)) > 0.5) {
  cat("   3. Variance is not constant - consider transformation (log or Box-Cox)\n")
}

if(lb_test$p.value < 0.05){
  cat("   4. Significant autocorrelation present - use ARIMA or similar models\n")
}

if(length(outliers_iqr) > 0 || length(outliers_z) > 0){
  cat("   5. Outliers detected - investigate and consider robust methods\n")
}

if(summary(linear_model)$coefficients[8] < 0.05) {
  cat("   6. Significant trend present - incorporate trend in forecasting models\n")
}

cat("\n")
cat("="=rep("=", 70), "\n")
cat("                      Analysis Complete\n")
cat("="=rep("=",70), "\n")

#============================================================================
# 15.- Save Results -------------------------------------------------------
#============================================================================



