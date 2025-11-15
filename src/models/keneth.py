from statsmodels.tsa.stattools import adfuller, kpss

def chequear_estacionaridad2(serie):
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



