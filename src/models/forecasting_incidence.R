###############################################################
# This script forecast the next 5 years of cancer incidence in 
# Europe by classification cancer, age group, and sex
###############################################################

library(tidyverse)
library(prophet)
library(stringr) # to use str_to_lowcase
library(lubridate) 

directory_path = "/home/jsancheg/git_environment/Europe_Cancer_Incidence/data/processed"


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

cat('--- Subseting the series to be forecasted by cancer type, sex, and age group ---')
item_id <- df_cancer_incidence_2007_2017 %>% 
  select(c("cancer_classification", "sex_character","age_group")) %>%
unique()


filtering_item <- function(df, item_row)
{
  ##################################################################
  # Filter the cases of cancer and return the records
  # df: a dataset variable containing the dataset with variables 
  #     cancer case, sex, and age group
  # item_row: a dataset containing the unique set of variables
  #          that identify each item in the series
  ##################################################################

    
  selected_series <- df %>% filter(cancer_classification == item_row[1,][[1]] &
                sex_character == item_row[1,][[2]] &
                age_group == item_row[1,][[3]])
  
  selected_series
  return(selected_series)
}
colnames(item_id)
colnames(df_cancer_incidence_2007_2017)
cat ('--- Subsetting the series ---')
serie_df <- filtering_item(df_cancer_incidence_2007_2017,
                           item_id[1,]) %>%
  select(ds, y)

print(serie_df)
# Checking  for missing dates and completeness
cat("Date range: ", min(serie_df$ds), " to ", max(serie_df$ds), "\n")
cat("Number of observations: ", nrow(serie_df), "\n")
cat("Missing values in y: ", sum(is.na(serie_df$y)), "\n")


#====================================================
# 3. Exploratory Data Analysis --------------------------------------------
#====================================================

# Plot the historical data
ggplot(serie_df, aes(x = ds, y = y)) +
  geom_line(color = "steelblue", linewidth  = 1) +
  geom_point(color = "steelblue", size = 1) +
  labs(
    title = "Historical Cancer Cases (Males, Age 8-4",
    x = "Year",
    y = "Number of Cases"
  ) +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

summary(serie_df$y)

#====================================================
# 4. Prophet model setup
#====================================================

model <- prophet(
  growth = 'linear',       # Can be 'linear' or 'logistic'
  yearly.seasonality = FALSE, # Set to TRUE if you suspect yearly patterns
  weekly.seasonality = FALSE, # False for annual data
  daily.seasonality = FALSE, # False for annual data
  seasonality.mode = 'additive', # 'additive' or 'multiplicative'
  changepoint.prior.scale = 0.05, # Flexibility of trend (0.001-0.5, default 0.05)'
                                  # Lower = less flexible, higher = more flexible
  seasonality.prior.scale = 10, # Flexibility of seasonality
  interval.width = 0.95,        # 95% prediction intervals
  uncertainty.samples = 1000    # Number of samples for uncertainty estimation
)

# Optimal: Add custom changepoints if you know when interventions occurred
# (e.g., screening program introduction, treatment advances)
# model <- prophet(
# changepoints = c("2025-01-01', '2010-01-01'),
# ---
# )

# Optional: Add regressors (external variables)
# if you have demographic or environmental factors:
# model <- add_regressor(model, 'population_size')
# model <- add_regressor(model, 'screening_rate')

#====================================================
# 5. Fit the model --------------------------------------------------------
#====================================================


cat("\nFitting Prophet model....\n")
model <- fit.prophet(model, serie_df)


#====================================================
# 6. Generate Foreast -----------------------------------------------------
# Create future dataframe for 5 years (5 observations for annual data)
#====================================================


future <- make_future_dataframe(
  model,
  periods = 5,  # 5 years ahead
  freq = 'year' # Annual frequency
)


# If you added regressors, you need to provide future values:
# future$population_size <- c(historical_values, future_values)
# future$screening_rate ,- c(historical_values, future_values)

# Generate predictions
cat("Generating forecasts ...\n")
forecast <- predict(model, future)

#====================================================
# 7. Examine Results ------------------------------------------------------
#====================================================


# View forecast future periods
future_forecast <- forecast %>%
  tail(5) %>%
  select(ds, yhat, yhat_lower, yhat_upper)

print("5-year Forecast")
print(future_forecast)

# Calculate prediction intervals
future_forecast <- future_forecast  %>%
  mutate(
    prediction_interval = yhat_upper - yhat_lower,
    relative_uncertainty = (yhat_upper - yhat_lower)/ yhat * 100
  )

print("\nForecast with Uncertainty Metrics")
print(future_forecast)

#====================================================
# 8. Visualization --------------------------------------------------------
#====================================================

# Plot 1: Forecast with components
plot(model, forecast) +
  labs(
    title = "Cancer Incidence Forecast (Males, Age 0 -4)",
    x = "Year",
    y = "Number of Cases"
  ) +
  theme_minimal()

# Plot 2: Decomposition (trend, seasonality if applicable)
prophet_plot_components(model, forecast)

