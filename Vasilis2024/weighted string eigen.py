import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.linalg import eigh # for symmetric eigenvalue problem

# Parameters
L = 1.0  # Length of string
Tension = 200.0  # String tension (N)
mu = 0.05     # Mass per unit length (kg/m)
N = 100      # Number of segments
dx = L / N
c = np.sqrt(Tension / mu) # Wave speed
dt = 0.001
total_time = 0.2

# Discretization
x = np.linspace(0, L, N+1)

# Mass matrix
m = mu * dx
M = np.eye(N+1) * m

# Add point masses at 1/4, 1/2, 3/4
extra_mass = 0.005 # in kg
for pos in [L / 4, L / 2, 3 * L / 4]:
  index = np.argmin(np.abs(x - pos))
  M[index, index] += extra_mass

# Stiffness matrix
K = np.zeros((N+1, N+1))
k = Tension / dx
for i in range(1, N):
  K[i, i] = 2 * k
  K[i, i-1] = -k
  K[i, i+1] = -k

# Fixed ends (Dirichlet boundary conditions)
K[0, :] = K[-1, :] = 0
K[:, 0] = K[:, -1] = 0
K[0, 0] = K[-1, -1] = 1e10 # effectively pinned

# Solve eigenvalue problem
eigvals, eigvecs = eigh(K, M)

# Keep first few modes
n_modes = 5
freqs = np.sqrt(eigvals[:n_modes]) / (2 * np.pi)
print (freqs)
modes = eigvecs[:, :n_modes]

# Initial condition: pluck shape (triangle)
pluck_pos = 0.5
y0 = np.where(x < pluck_pos,
       x / pluck_pos,
       (L - x) / (L - pluck_pos))

# Project initial shape onto modes
a_n = modes.T @ (M @ y0)

# Time evolution
t_vals = np.linspace(0, total_time, 500)

def compute_mode_sum(t):
  y = np.zeros_like(x)
  for n in range(n_modes):
    omega_n = np.sqrt(eigvals[n])
    y += a_n[n] * modes[:, n] * np.cos(omega_n * t)
  return y

# Plotting
fig, ax = plt.subplots()
line, = ax.plot(x, compute_mode_sum(0))
for wp in [L / 4, L / 2, 3 * L / 4]:
  ax.axvline(x=wp, color='red', linestyle='--', alpha=0.5, label='Weight' if wp == L/4 else None)
ax.set_ylim(-1.2, 1.2)
ax.set_title("Vibration of String with Point Masses")
ax.set_xlabel("Position along string")
ax.set_ylabel("Displacement")
plt.legend()

def animate(i):
  y = compute_mode_sum(t_vals[i])
  line.set_ydata(y)
  return line,

ani = animation.FuncAnimation(fig, animate, frames=len(t_vals), interval=20, blit=True)
plt.show()

# Υπολογισμός ζευγών (k, ω) από τα modes
n_vals = np.arange(1, n_modes+1)
k_vals = n_vals * np.pi / L
omega_vals = 2 * np.pi * freqs  # από το eigensolver

# Σχεδίαση διασποράς
plt.figure()
plt.plot(k_vals, omega_vals, 'o-', label="Με σημειακές μάζες (προσομοίωση)")
plt.plot(k_vals, c * k_vals, '--', label="Ιδανική χορδή (ω = c k)")
plt.xlabel("Κυματικός αριθμός k [rad/m]")
plt.ylabel("Γωνιακή συχνότητα ω [rad/s]")
plt.title("Διασπορά χορδής με σημειακές μάζες")
plt.legend()
plt.grid(True)
plt.show()
# --- Υπολογισμός σήματος σε ένα σημείο της χορδής ---
point_index = N // 4   # στη μέση της χορδής
signal = [compute_mode_sum(t)[point_index] for t in t_vals]

# --- FFT ---
fft_vals = np.fft.fft(signal)
fft_freqs = np.fft.fftfreq(len(signal), d=(t_vals[1]-t_vals[0]))

# Κρατάμε μόνο τις θετικές συχνότητες
pos_mask = fft_freqs > 0
fft_freqs = fft_freqs[pos_mask]
fft_vals = np.abs(fft_vals[pos_mask])

# --- Plot FFT ---
plt.figure()
plt.plot(fft_freqs, fft_vals)
plt.xlim(0, 500)  # περιορίζουμε το εύρος για καθαρότητα
plt.title("Φάσμα Συχνοτήτων (FFT) στο μέσο της χορδής")
plt.xlabel("Συχνότητα (Hz)")
plt.ylabel("Ένταση")
plt.grid(True)
plt.show()

# --- Υπολογισμός v_phase και v_group για τους πρώτους n_modes ---
# Χρησιμοποιούμε: omega_n = 2π f_n και k_n = n π / L (δεικτικός κυματικός αριθμός ανά τρόπο)
n_vals = np.arange(1, n_modes+1)
k_vals = n_vals * np.pi / L
omega_vals = 2 * np.pi * freqs  # rad/s

# Phase velocity: v_phase = ω / k
v_phase = omega_vals / k_vals

# Group velocity: v_group ≈ dω/dk (με πεπερασμένες διαφορές)
v_group = np.gradient(omega_vals, k_vals)

# --- Διάγραμμα ταχυτήτων ---
plt.figure()
plt.plot(k_vals, v_phase, 'o-', label='v_phase (ω/k)')
plt.plot(k_vals, v_group, 's-', label='v_group (dω/dk)')
plt.axhline(np.sqrt(Tension/mu), color='gray', linestyle='--', alpha=0.6, label='c ιδανικής χορδής')
plt.xlabel("Κυματικός αριθμός k (rad/m)")
plt.ylabel("Ταχύτητα (m/s)")
plt.title("Ταχύτητα φάσης και ομάδας ανά τρόπο")
plt.grid(True)
plt.legend()
plt.show()
