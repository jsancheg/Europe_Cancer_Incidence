
from statsmodels.tsa.stattools import adfuller, kpss
def chequear_estacionaridad1(serie):
    print("===== Stationary Tests =======")
    adf_result = adfuller(serie, autolag='AIC')
    adf_stat = adf_result[0]
    adf_pvalue = adf_result[1]
    print(f"ADF Statistic: {adf_stat:.4f}")
    print(f"P-value: {adf_pvalue:.4f}")
    adf_interpretation = "Series is STATIONARY (reject H0)" if adf_pvalue < 0.05 else "Series is NON-STATIONARY (fail to reject H0)"
    
    adf_resumen = {"adf_estadistico": [adf_stat],
                   "adf_valorp": [adf_pvalue],
                   "adf_interpretacion": adf_interpretation}
    
    return adf_resumen    


def 

