import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.tsa.stattools import adfuller, kpss
from statsmodels.graphics.tsaplots import plot_acf, plot_pacf
from statsmodels.stats.diagnostic import acorr_ljungbox
from io import StringIO




def plot_serie(ts_data, title, xlabel, ylabel):
    plt.figure(figsize=(10, 4))
    ts_data.plot(title=title, linewidth=2, marker='o')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.show()



def descriptive_statistics(ts_data):
    ds_dict = {'Mean': [ts_data.mean()],
               'Median': [ts_data.median()],
               'SD': [ts_data.std()],
               'Min': [ts_data.min()],
               'Max': [ts_data.max()],
               'Range': [ts_data.max() - ts_data.min()],
               'CV': [(ts_data.std() / ts_data.mean() ) *100]
    }
    print ("==== Descriptive Statistics ===")

    return ds_dict

def chequear_estacionaridad(serie):
    print("===== Stationary Tests =======")
    adf_result = adfuller(serie, autolag='AIC')
    adf_stat = adf_result[0]
    adf_pvalue = adf_result[1]
    print(f"ADF Statistic: {adf_stat:.4f}")
    print(f"P-value: {adf_pvalue:.4f}")
    adf_interpretation = "Series is STATIONARY (reject H0)" if adf_pvalue < 0.05 else "Series is NON-STATIONARY (fail to reject H0)"
    
    print("--- KPSS Test ---")
    # R's kpss.test(ts, null="Trend") corresponds to regression='ct' in Python.
    # R's default null="Level" corresponds to regression='c' in Python.
    # The R script uses 'Trend', so we'll use 'ct'
    kpss_result = kpss(serie, regression='ct', nlags='auto')
    kpss_stat = kpss_result[0]
    kpss_pvalue = kpss_result[1] # Note: Python's KPSS returns a truncated p-value
    kpss_critical = kpss_result[3]
    print(f"Kpss Statistics: {kpss_stat:.4f}")
    # The p-value for KPSS is often not a precise value, but a range.
    # We'll compare to a common significance level (e.g., 5%)
    is_kpss_stationary = kpss_stat < kpss_critical['5%']
    kpss_interpretation = "Series is STATIONARY (fail to reject H0)" if is_kpss_stationary else "Series is NON-STATIONARY (reject H0)"
    print(f"P-value: approx. {kpss_pvalue}")
    print(f"Interpretation: {kpss_interpretation}\n")

    resumen = { "test_name": ["dikeyfuller", "kpss"],
                    "estadistico": [adf_stat, kpss_stat],
                   "valorp": [adf_pvalue, kpss_pvalue],
                   "interpretacion": [adf_interpretation, kpss_interpretation] }
    
    return resumen    

def chequear_autocorrelacion(serie, alpha = 0.05, number_lags = 10):
    print("=== Autocorrelation Analysis ===\n")

    # 5.1 ACF and PACF plots
    #
    fig, axes = plt.subplots(2, 1, figsize=(10, 6))
    plot_acf(serie, lags=min(20, len(serie)//2 - 1), ax=axes[0], title="Autocorrelation Function (ACF)")
    plot_pacf(serie, lags=min(20, len(serie)//2 - 1), ax=axes[1], title="Partial Autocorrelation Function (PACF)")
    plt.tight_layout()
    #plt.show()

    # 5.2 Ljung-Box Test for autocorrelation
    print("--- Ljung-Box Test (lag = 10) ---")
    # statsmodels' acorr_ljungbox returns a DataFrame, lag is the number of lags tested.
    lb_result = acorr_ljungbox(serie, lags=[number_lags])
    lb_stat = lb_result.loc[10, 'lb_stat']
    lb_pvalue = lb_result.loc[10, 'lb_pvalue']
    print(f"Chi-squared: {lb_stat:.4f}")
    print(f"P-value: {lb_pvalue:.4f}")
    lb_interpretation = "Significant autocorrelation detected" if lb_pvalue < 0.05 else "No significant autocorrelation"
    print(f"Interpretation: {lb_interpretation}\n")
    lb_result['interpretacion'] = lb_interpretation
    
    return fig, axes, lb_result

def calculate_rolling_statistics(serie, window_size = None):
    print("\n=== Rolling Statistics Analysis ===")

    warning =  ""
    recommendation = ""
    variance_check_result = ""

    # Calculate rolling mean and standard deviation
    if window_size is None:
        window_size = min(5, len(serie) // 3)

    rolling_mean = serie.rolling(window = window_size).mean()
    rolling_sd = serie.rolling(window = window_size).std()


    fig, axes = plt.subplots(2,1, figsize = (10, 8))

    # Plot 1: Time Series and Rolling Mean
    serie.plot(ax = axes[0],
               label = 'Original',
               linewidth = 0.8,
               color = 'steelblue',
               marker = 'o'
               )
    rolling_mean.plot(ax = axes[0],
                      label = 'Rolling Mean',
                      linewidth = 1.2,
                      color = 'darkorange'
                      )
    axes[0].set_title(f"Time series with {window_size}--Period Rolling Mean")
    axes[0].set_xlabel('Time')
    axes[0].set_ylabel('Cases')
    axes[0].legend(title = 'Series')
    axes[0].grid(True, linestyle = '--', alpha = 0.6)

    # Plot 2: Rolling Standard Deviation
    rolling_sd.plot(ax = axes[1],
                    linewidth = 1,
                    color = 'darkred'
                    )
    axes[1].set_title(f"{window_size}--Period Standard Deviation")
    axes[1].set_xlabel('Time')
    axes[1].set_ylabel('Standard Deviation')
    axes[1].grid(True, linestyle = '--', alpha = 0.6)

    plt.tight_layout()

    # Test for changing variance
    print("\nVariance stability check:")
    half_len = len(serie) // 2
    first_half_var = serie.iloc[:half_len].var()
    second_half_var = serie.iloc[half_len:].var()
    if first_half_var > 0:
        variance_ratio = second_half_var / first_half_var
        print(f"Variance ratio: {variance_ratio:.2f}")
        if np.abs(np.log(variance_ratio)) > 0.5:
            warning = "WARNING: Variance appears to be changing over time"
            print(f"{warning}\n")
            recommendation = "RECOMMENDATION: Consider log transformation or variance stabilization"
            print(f"{recommendation}")
        else:
            variance_check_result = "Variance ratio is within acceptable limits."
    else:
        variance_check_result = "Cannot compute variance ratio (First half variance is zero or negative"

    if variance_check_result:
        print("-"*30)
        print(variance_check_result)

    

    df_rolling = pd.DataFrame({'window_size': [window_size],
                               'half_len': [half_len],
                               'first_half_var': [first_half_var],
                               'second_half_var':[second_half_var],
                               'variance_ratio': [variance_ratio],
                               'warning':[warning],
                               'recommendation':[recommendation],
                               'variance_check_result': [variance_check_result]})

    return fig, axes, rolling_mean, rolling_sd, df_rolling







    