import yfinance as yf
import pandas as pd
import statsmodels.api as sm
from statsmodels.tsa.stattools import adfuller
from itertools import combinations
import warnings

warnings.filterwarnings("ignore")

def scan_for_cointegration(tickers, start_date, end_date, p_value_threshold=0.05):
    print(f"Downloading data for {len(tickers)} ETFs...")
    data = yf.download(tickers, start=start_date, end=end_date)['Close']
    data = data.ffill().dropna()
    
    pairs = list(combinations(data.columns, 2))
    print(f"Testing {len(pairs)} unique ETF pairs...\n")
    
    results = []
    
    for asset1, asset2 in pairs:
        y = data[asset1]
        x = data[asset2]
        
        try:
            x_with_constant = sm.add_constant(x)
            model = sm.OLS(y, x_with_constant).fit()
            hedge_ratio = model.params.iloc[1]
            
            spread = y - (hedge_ratio * x)
            p_value = adfuller(spread)[1]
            
            if p_value < p_value_threshold:
                results.append({
                    'Asset 1': asset1,
                    'Asset 2': asset2,
                    'P-Value': p_value,
                    'Hedge Ratio': hedge_ratio
                })
        except Exception:
            continue
            
    results_df = pd.DataFrame(results)
    if not results_df.empty:
        results_df = results_df.sort_values(by='P-Value').reset_index(drop=True)
    
    return results_df

if __name__ == "__main__":
    # The "Easy to Borrow" Macro Universe
    etf_universe = [
        # --- US Major Indices ---
        'SPY',  # S&P 500
        'QQQ',  # Nasdaq 100 (Tech heavy)
        'DIA',  # Dow Jones Industrial Average
        'IWM',  # Russell 2000 (Small Caps)
        
        # --- International Indices ---
        'EWU',  # United Kingdom (FTSE Proxy)
        'EWG',  # Germany (DAX Proxy)
        'EWJ',  # Japan (Nikkei Proxy)
        'EEM',  # Emerging Markets Broad
        'EFA',  # Developed Markets ex-US/Canada
        
        # --- Commodities (Metals & Energy) ---
        'GLD',  # Physical Gold
        'SLV',  # Physical Silver
        'USO',  # US Crude Oil
        'DBC',  # Broad Commodity Index
        
        # --- US Treasury Bonds (Interest Rates) ---
        'TLT',  # 20+ Year Treasury Bonds
        'IEF',  # 7-10 Year Treasury Bonds
        'SHY'   # 1-3 Year Treasury Bonds
    ]
    
    found_pairs = scan_for_cointegration(
        tickers=etf_universe, 
        start_date="2020-01-01", 
        end_date="2024-01-01", 
        p_value_threshold=0.05
    )
    
    print("--- TOP COINTEGRATED MACRO PAIRS ---")
    if found_pairs.empty:
        print("No cointegrated pairs found.")
    else:
        found_pairs['P-Value'] = found_pairs['P-Value'].apply(lambda x: f"{x:.5f}")
        print(found_pairs.to_string())