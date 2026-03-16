import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# --- Nd:YAG Laser Rate Equations Simulation ---

# 1. Define the System Parameters (Nd:YAG Typical Values)
# These values are characteristic of a continuously pumped Nd:YAG laser
# operating above threshold (Class B laser).

# Physical constants and typical laser parameters
c = 3.0e8      # Speed of light in vacuum (m/s)
n_g = 1.82     # Group refractive index of YAG
L = 0.3        # Cavity length (m)
A_eff = 1.0e-7 # Effective mode cross-sectional area (m^2)
sigma = 2.8e-23 # Stimulated emission cross-section (m^2)

# Time constants
tau_f = 230e-6 # Upper level fluorescence lifetime (N2 decay) in seconds (Nd:YAG)
tau_p = 0.5e-7 # Photon lifetime in the cavity (seconds). 
               # Note: tau_f >> tau_p is the condition for oscillations.

# Normalized Pumping Rate (R_p)
# R_p = W_p * N_total, where W_p is the pumping probability (s^-1)
# R_th is the minimum rate required to reach threshold
R_th = 1.0 / tau_f + (c / (n_g * L)) * (1/tau_p) * (1 / (sigma * L))
R_p = 1.5 * R_th # Pumping rate (R_p) set to 1.5 times the threshold rate

# 2. Define the Coupled Rate Equations
# These two differential equations describe the change in N2 (population inversion)
# and phi (photon density) over time.

def rate_equations(t, y):
    """
    Coupled Rate Equations for a Four-Level Laser (Nd:YAG).
    
    y[0] = N2: Population density of the upper laser level (inversion).
    y[1] = phi: Photon density in the cavity.
    
    Returns: [dN2/dt, dphi/dt]
    """
    N2 = y[0]
    phi = y[1]

    # Equation 1: Change in Upper Laser Level Population (N2)
    # dN2/dt = (Pumping) - (Spontaneous Decay) - (Stimulated Emission)
    # R_p is the effective pumping rate into the N2 level.
    # The inversion N is approximated by N2 for a four-level system (N1 ~ 0).
    dN2_dt = R_p - (N2 / tau_f) - (c * sigma * N2 * phi / n_g)

    # Equation 2: Change in Photon Density (phi)
    # dphi/dt = (Stimulated Emission) + (Spontaneous Emission) - (Cavity Loss)
    # The term 1/tau_p includes all cavity losses (output coupling, scattering, etc.)
    # The spontaneous emission factor (N2 / (tau_f * V_mode)) is often ignored or simplified.
    # We use a simplified model focusing only on the two key terms for RO:
    # dphi/dt = (Gain - Loss) * phi
    
    # We use a slightly more accurate standard form that includes the fraction of 
    # spontaneous emission coupling into the lasing mode (beta * N2 / tau_f).
    beta = 1e-4 # Spontaneous emission factor (fraction coupled into the mode)
    dphi_dt = (c * sigma * N2 / n_g - 1 / tau_p) * phi + (beta * N2 / tau_f)
    
    return [dN2_dt, dphi_dt]

# 3. Setup the Solver and Initial Conditions
N2_initial = 1.0e10  # Start with a small population inversion
phi_initial = 1.0e-3 # Start with a few noise photons
initial_conditions = [N2_initial, phi_initial]

# Time span for simulation (1.5 milliseconds)
t_span = [0, 1.5e-3]
time_points = np.linspace(t_span[0], t_span[1], 1000)

# Solve the differential equations
# Method 'RK45' (Runge-Kutta 4th and 5th order) is a good general purpose solver.
solution = solve_ivp(
    rate_equations, 
    t_span, 
    initial_conditions, 
    t_eval=time_points,
    method='RK45'
)

# Extract results
time = solution.t
N2 = solution.y[0]
phi = solution.y[1]

# Calculate Laser Output Power (proportional to phi)
# P_out is proportional to the photon density phi * (1/tau_p) * (hbar * omega) * V_mode
# We will plot the normalized photon density for simplicity, as it shows the oscillation directly.
laser_output_power = phi / np.max(phi) 


# 4. Plot the Results
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
fig.suptitle(f'Nd:YAG Laser Relaxation Oscillations (Pump $1.5 \times$ Threshold)', fontsize=14)

# --- Plot 1: Population Inversion (N2) ---
ax1.plot(time * 1e3, N2, label='$N_2$ (Upper Level Population)', color='#4169E1', linewidth=2)
# Add the steady-state threshold line
N_th = n_g / (c * sigma * tau_p)
ax1.axhline(N_th, color='r', linestyle='--', label='$N_{th}$ (Threshold Inversion)')
ax1.set_ylabel('Population Inversion ($N_2$) [a.u.]', color='#4169E1')
ax1.tick_params(axis='y', labelcolor='#4169E1')
ax1.grid(True, linestyle=':', alpha=0.6)
ax1.legend(loc='upper right')

# --- Plot 2: Photon Population (Phi) / Output Power ---
ax2.plot(time * 1e3, phi * 1e-12, label='(Photon Density)', linewidth=2)
ax2.set_xlabel('Time (ms)')
ax2.set_ylabel('Photon Density ')
ax2.tick_params(axis='y')
ax2.set_ylim(bottom=0)
ax2.grid(True, linestyle=':', alpha=0.6)
ax2.legend(loc='upper right')

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()

# 5. Summary of Key Observation Points
print("\n--- Key Simulation Results ---")
print(f"Upper Level Lifetime (tau_f): {tau_f*1e6:.1f} µs")
print(f"Photon Lifetime (tau_p): {tau_p*1e9:.1f} ns")
print(f"Pumping Rate R_p is {R_p/R_th:.1f} times the threshold rate (R_th).")
print("\nLook at the graph to see the following:")
print("1. Overshoot: N2 rises above N_th before the photons kick in.")
print("2. Spiking: Phi rapidly rises and then crashes, pulling N2 down below N_th.")
print("3. Damping: The oscillations in both N2 and Phi decrease over time until they reach a steady-state value.")