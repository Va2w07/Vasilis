import yfinance as yf
import pandas as pd
import statsmodels.api as sm
from statsmodels.tsa.stattools import adfuller
import matplotlib.pyplot as plt

# 1. Define parameters and fetch historical data
asset1 = "PLD"
asset2 = "PPLT"
print(f"Downloading data for {asset1} and {asset2}...")
data = yf.download([asset1, asset2], start="2022-01-01", end="2026-06-01")

# Extract the 'Adjusted Close' prices
prices = data['Close'].dropna()
y = prices[asset1]
x = prices[asset2]

# 2. Calculate the Hedge Ratio using Ordinary Least Squares (OLS) regression
# We want to find the relationship: Asset1 = (Hedge_Ratio * Asset2) + Intercept + Error
x_with_constant = sm.add_constant(x)
model = sm.OLS(y, x_with_constant).fit()
hedge_ratio = model.params.iloc[1]

print(f"\nCalculated Hedge Ratio: {hedge_ratio:.4f}")
print(f"For every 1 share of {asset1} bought, you short {hedge_ratio:.4f} shares of {asset2}.")

# 3. Construct the stationary Spread
spread = y - (hedge_ratio * x)

# 4. Run the Augmented Dickey-Fuller (ADF) Test
# The Null Hypothesis is that the spread is a random walk (non-stationary).
# A p-value < 0.05 means we reject the null hypothesis and assume it is cointegrated.
adf_result = adfuller(spread)
p_value = adf_result[1]
test_statistic = adf_result[0]

print("\n--- Augmented Dickey-Fuller Test Results ---")
print(f"ADF Test Statistic: {test_statistic:.4f}")
print(f"P-Value: {p_value:.4f}")

if p_value < 0.05:
    print("Result: STRONG evidence of cointegration. The spread is mean-reverting.")
else:
    print("Result: WEAK evidence of cointegration. The spread may be a random walk.")

# 5. Calculate the Z-Score for trading signals
# This normalizes the spread so we can set fixed standard deviation thresholds (+2 / -2)
z_score = (spread - spread.mean()) / spread.std()

# 6. Data Visualization
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

# Plot 1: The Raw Spread
ax1.plot(spread.index, spread, color='midnightblue', label='Price Spread')
ax1.axhline(spread.mean(), color='red', linestyle='--', label='Mean')
ax1.set_title(f'Raw Spread between {asset1} and {asset2}')
ax1.set_ylabel('Spread Value ($)')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Plot 2: The Normalized Z-Score with Trading Bands
ax2.plot(z_score.index, z_score, color='teal', label='Z-Score')
ax2.axhline(0, color='black', linestyle='-', alpha=0.5)
ax2.axhline(2.0, color='red', linestyle='--', label='+2.0 Std Dev (Short Spread)')
ax2.axhline(-2.0, color='green', linestyle='--', label='-2.0 Std Dev (Buy Spread)')
ax2.set_title('Normalized Z-Score (Trading Signals)')
ax2.set_ylabel('Standard Deviations')
ax2.set_xlabel('Date')
ax2.legend()
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()