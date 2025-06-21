import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.linalg import eigh  # for symmetric eigenvalue problem

# Parameters
L = 1.0           # Length of string
Tension = 200.0   # String tension (N)
mu = 0.005         # Mass per unit length (kg/m)
N = 100           # Number of segments
dx = L / N
c = np.sqrt(Tension / mu)  # Wave speed
dt = 0.001
total_time = 0.2

# Discretization
x = np.linspace(0, L, N+1)

# Mass matrix
m = mu * dx
M = np.eye(N+1) * m

# Add point masses at 1/4, 1/2, 3/4
extra_mass = 0.005  # in kg
for pos in [L/6, 2*L/6, 3*L/6, 4*L/6, 5*L/6, L]:
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
K[0, 0] = K[-1, -1] = 1e10  # effectively pinned

# Solve eigenvalue problem
eigvals, eigvecs = eigh(K, M)

# Keep first few modes
n_modes = 10
freqs = np.sqrt(eigvals[:n_modes]) / (2 * np.pi)
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
for wp in [L/6, 2*L/6, 3*L/6, 4*L/6, 5*L/6]:
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