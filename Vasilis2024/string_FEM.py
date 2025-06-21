import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Physical parameters
L = 1.0                # Length of string (m)
T = 100.0              # Tension (N)
mu = 0.01              # Mass per unit length (kg/m)
N = 100                # Number of elements (N+1 nodes)
dx = L / N             # Element length
c = np.sqrt(T / mu)    # Wave speed
dt = 0.0005            # Time step
t_max = 0.5            # Total simulation time
steps = int(t_max / dt)

# Discretize space
x = np.linspace(0, L, N + 1)

# Mass and stiffness matrices (lumped mass and linear stiffness)
M = np.eye(N + 1) * mu * dx  # Lumped mass matrix
K = np.zeros((N + 1, N + 1))

# Build stiffness matrix (tridiagonal)
for i in range(1, N):
    K[i, i] += 2 * T / dx
    K[i, i - 1] -= T / dx
    K[i, i + 1] -= T / dx

# Apply boundary conditions (fixed ends)
K[0, :] = K[-1, :] = 0
K[:, 0] = K[:, -1] = 0
K[0, 0] = K[-1, -1] = 1e12  # Large number to pin ends
M[0, 0] = M[-1, -1] = 1e12  # Prevent division by zero

# Initial conditions
u = np.zeros(N + 1)       # Displacement
v = np.zeros(N + 1)       # Velocity
a = np.zeros(N + 1)       # Acceleration

# Driving force parameters
f_drive = 120.0                 # Frequency in Hz
omega_drive = 2 * np.pi * f_drive
F_amp = 5.0                     # Force amplitude
f_ext = np.zeros(N + 1)

drive_index = N // 2            # Middle of the string

# Precompute inverse of M for efficiency
M_inv = np.linalg.inv(M)

# Store results for animation
history = []

for step_i in range(steps):
    t = step_i * dt

    # External force: harmonic at center
    f_ext[:] = 0
    f_ext[drive_index] = F_amp * np.sin(omega_drive * t)

    # Central difference time integration
    # a = M^(-1) * (F - K * u)
    a = M_inv @ (f_ext - K @ u)

    # Velocity-Verlet style update
    u_new = u + dt * v + 0.5 * dt ** 2 * a
    a_new = M_inv @ (f_ext - K @ u_new)
    v_new = v + 0.5 * dt * (a + a_new)

    u = u_new
    v = v_new

    if step_i % 10 == 0:
        history.append(u.copy())

# Animate result
fig, ax = plt.subplots()
line, = ax.plot(x, history[0])
ax.set_ylim(-1.5, 1.5)
ax.set_xlabel("Position along string")
ax.set_ylabel("Displacement")
ax.set_title("Driven Vibrating String (FEM)")

def animate(i):
    line.set_ydata(history[i])
    return line,

ani = animation.FuncAnimation(fig, animate, frames=len(history), interval=20, blit=True)
plt.show()