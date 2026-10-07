import yfinance as yf
import pandas as pd
import statsmodels.api as sm
from statsmodels.tsa.stattools import adfuller
from itertools import combinations
import warnings

# Suppress warnings from statsmodels to keep the console clean
warnings.filterwarnings("ignore")

def scan_for_cointegration(tickers, start_date, end_date, p_value_threshold=0.05):
    print(f"Downloading data for {len(tickers)} assets...")
    # Download all data in bulk to save time
    data = yf.download(tickers, start=start_date, end=end_date)['Close']
    
    # Forward-fill and drop remaining NaNs to ensure all series are the same length
    data = data.ffill().dropna()
    
    # Generate all unique pairs (n * (n-1) / 2 combinations)
    pairs = list(combinations(data.columns, 2))
    print(f"Testing {len(pairs)} unique pairs for cointegration...\n")
    
    results = []
    
    for asset1, asset2 in pairs:
        y = data[asset1]
        x = data[asset2]
        
        try:
            # 1. Calculate the Hedge Ratio using OLS
            x_with_constant = sm.add_constant(x)
            model = sm.OLS(y, x_with_constant).fit()
            hedge_ratio = model.params.iloc[1]
            
            # 2. Construct the stationary Spread
            spread = y - (hedge_ratio * x)
            
            # 3. Run the ADF Test
            adf_result = adfuller(spread)
            p_value = adfuller(spread)[1]
            
            # 4. Store the results if it passes our threshold
            if p_value < p_value_threshold:
                results.append({
                    'Asset 1': asset1,
                    'Asset 2': asset2,
                    'P-Value': p_value,
                    'Hedge Ratio': hedge_ratio
                })
        except Exception as e:
            # Skip pairs that cause calculation errors (e.g., flat lines, insufficient data)
            continue
            
    # Convert results to a DataFrame for easy sorting and viewing
    results_df = pd.DataFrame(results)
    
    if not results_df.empty:
        # Sort by the lowest P-Value (strongest evidence of cointegration)
        results_df = results_df.sort_values(by='P-Value').reset_index(drop=True)
    
    return results_df

if __name__ == "__main__":
    # Define your universe of assets
    # Mix of Tech, Semiconductors, Finance, and Indices
    my_universe = [
        'AAPL', 'MSFT', 'GOOG', 'META', 'AMZN', 
        'NVDA', 'AMD', 'INTC', 'TSM',
        'V', 'MA', 'JPM', 'BAC', 'GS',
        'SPY', 'QQQ', 'DIA'
    ]
    
    # Run the scanner
    # Note: Scanning over a 3-4 year period is standard to capture full business cycles
    found_pairs = scan_for_cointegration(
        tickers=my_universe, 
        start_date="2020-01-01", 
        end_date="2024-01-01", 
        p_value_threshold=0.05
    )
    
    print("--- TOP COINTEGRATED PAIRS ---")
    if found_pairs.empty:
        print("No cointegrated pairs found in this universe.")
    else:
        # Format the P-Value for cleaner reading
        found_pairs['P-Value'] = found_pairs['P-Value'].apply(lambda x: f"{x:.5f}")
        print(found_pairs.to_string())
        
    print("\nWarning: Always check for a fundamental relationship!")