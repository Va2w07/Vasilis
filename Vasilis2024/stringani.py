import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# String properties
L = 1.0  # Length of string (meters)
T = 0.02  # Total duration of simulation (seconds)
c = 100.0  # Wave speed (m/s)
n_modes = 8  # Number of harmonics

# Discretization
x = np.linspace(0, L, 500)
t = np.linspace(0, T, 500)

# Initial pluck shape (triangle pluck)
pluck_position = 0.3  # Where the string is plucked (as a fraction of L)
y0 = np.where(x < pluck_position * L,
              x / (pluck_position * L),
              (L - x) / (L - pluck_position * L))

# Fourier sine series coefficients for the pluck shape
def fourier_coeffs(y0, n_terms):
    a_n = []
    for n in range(1, n_terms + 1):
        integrand = y0 * np.sin(n * np.pi * x / L)
        a = 2 / L * np.trapz(integrand, x)
        a_n.append(a)
    return np.array(a_n)

a_n = fourier_coeffs(y0, n_modes)

# Create figure
fig, ax = plt.subplots()
line, = ax.plot(x, y0)
ax.set_ylim(-1.2, 1.2)
ax.set_xlabel("Position along string")
ax.set_ylabel("Displacement")
ax.set_title("Guitar String Vibration")

# Animation function
def animate(i):
    t_i = t[i]
    y = np.zeros_like(x)
    for n in range(1, n_modes + 1):
        y += a_n[n-1] * np.sin(n * np.pi * x / L) * np.cos(n * np.pi * c * t_i / L)
    line.set_ydata(y)
    return line,

ani = animation.FuncAnimation(fig, animate, frames=len(t), interval=20, blit=True)

plt.show()