import matplotlib.pyplot as plt
from statsmodels.tsa.stattools import adfuller, kpss
from statsmodels.graphics.tsaplots import plot_acf, plot_pacf
from statsmodels.stats.diagnostic import acorr_ljungbox


# Custom function to replace the missing R function
def filtering_item(df, item):
    """
    Simulates the filtering_item function by filtering the DataFrame
    based on 'cancer_classification', 'sex_character', and 'age_group'.
    """
    # Assuming item is a pandas Series or single-row DataFrame from item_id
    filtered_df = df[
        (df['cancer_classification'] == item['cancer_classification'].iloc[0]) &
        (df['sex_character'] == item['sex_character'].iloc[0]) &
        (df['age_group'] == item['age_group'].iloc[0])
    ].copy()
    return filtered_df



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

def chequear_autocorrelacion(serie, alpha = 0.05):
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
    lb_result = acorr_ljungbox(serie, lags=[10])
    lb_stat = lb_result.loc[10, 'lb_stat']
    lb_pvalue = lb_result.loc[10, 'lb_pvalue']
    print(f"Chi-squared: {lb_stat:.4f}")
    print(f"P-value: {lb_pvalue:.4f}")
    lb_interpretation = "Significant autocorrelation detected" if lb_pvalue < 0.05 else "No significant autocorrelation"
    print(f"Interpretation: {lb_interpretation}\n")
    lb_result['interpretacion'] = lb_interpretation
    return fig, axes, lb_result