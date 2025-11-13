


filtering_item <- function(df, item_row)
{
  ##################################################################
  # Filter the cases of cancer and return the records
  # df: a dataset variable containing the dataset with variables 
  #     cancer case, sex, and age group
  # item_row: a dataset containing the unique set of variables
  #          that identify each item in the series
  ##################################################################
  
  
  selected_series <- df %>% filter(cancer_classification == item_row[[1]] &
                                     sex_character == item_row[[2]] &
                                     age_group == item_row[[3]])
  
  selected_series
  return(selected_series)
}

# Basic time series plot
plot_basic <- function(ts_obj, title = "Time Series Plot"){
  autoplot(ts_obj) + 
    geom_point(color = "steelblue", size = 2) +
    geom_line(color = "steelblue", linewidth  = 1) +
    labs(
      title = title,
      x = "Time",
      y = "Cases"
    ) +
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}
