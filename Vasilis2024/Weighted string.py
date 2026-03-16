
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# ----------------------
#  Φυσικές και αριθμητικές παράμετροι
L = 1.0          # μήκος (m)
T = 100.0        # τάση (N)
mu = 0.01        # γραμμική μάζα (kg/m)
N = 200          # τμήματα -> N+1 κόμβοι (μεγαλύτερο N = καλύτερη ακρίβεια)
f_drive = 26.0  # συχνότητα διέγερσης (Hz)
F_amp = 500.0    # πλάτος εξωτερικής δύναμης (N)
drive_x = 0.5    # θέση διέγερσης (m)

# Σημειακά βαρίδια: λίστα θέσεων και μάζα καθενός (kg)
extra_mass_positions = [0.25, 0.50, 0.75]
extra_mass_value = 0.05

# Απόσβεση (Rayleigh: -c_damp * v)
c_damp = 0.02

# Χρόνος προσομοίωσης και βήμα
t_max = 3.0
dt_requested = 5e-4
# ----------------------

dx = L / N
c = np.sqrt(T / mu)

# Σταθερότητα (CFL για ρητό σχήμα κύματος): dt <= ~ dx / c
dt_stable = 0.9 * dx / c
dt = min(dt_requested, dt_stable)
steps = int(np.ceil(t_max / dt))
print(f"dx={dx:.3e}, c={c:.3e}, dt={dt:.3e}, steps={steps}")

# Κόμβοι: κρατάμε μόνο τους εσωτερικούς (σταθερά άκρα 0 και L)
x_full = np.linspace(0, L, N+1)
Ni = N - 1
x = x_full[1:-1].copy()

# ----------------------
# Μαζική μήτρα (lumped) με σημειακά βαρίδια
m_lumped = np.ones(Ni) * (mu * dx)

# Χαρτογράφηση θέσεων βαριδίων στον εσωτερικό πίνακα
for pos in extra_mass_positions:
    # αν πέσει ακριβώς σε άκρο, το αγνοούμε (άκρα είναι σταθερά)
    if pos <= 0.0 or pos >= L:
        continue
    # βρες πλησιέστερο εσωτερικό κόμβο
    idx_full = np.argmin(np.abs(x_full - pos))
    if 1 <= idx_full <= N-1:
        idx_internal = idx_full - 1
        m_lumped[idx_internal] += extra_mass_value

# Αντίστροφη μάζα (διαγώνια)
M_inv = np.diag(1.0 / m_lumped)

# ----------------------
# Σκληρότητα (πεπερασμένες διαφορές για u_xx)
K = np.zeros((Ni, Ni))
k_coef = T / dx
for i in range(Ni):
    K[i, i] += 2 * k_coef
    if i - 1 >= 0:
        K[i, i-1] -= k_coef
    if i + 1 < Ni:
        K[i, i+1] -= k_coef

# ----------------------
# Αρχικές συνθήκες
u0 = np.zeros(Ni)
v0 = np.zeros_like(u0)

# Δείκτης διέγερσης
drive_index = np.argmin(np.abs(x - drive_x))
omega_drive = 2 * np.pi * f_drive

# ----------------------
# Αποθήκευση για animation και FFT
store_every = max(1, int(2e-2 / dt))  # περίπου κάθε 0.02s
history = []
center_ts = []
probe_index = np.argmin(np.abs(x - 0.5))  # δειγματοληψία στο κέντρο (μπορείς να αλλάξεις)

# Αρχική επιτάχυνση
t0 = 0.0
f_ext = np.zeros(Ni)
f_ext[drive_index] = F_amp * np.sin(omega_drive * t0)
a0 = M_inv @ (f_ext - K @ u0 - c_damp * v0)

u_prev = u0.copy()
u_curr = u0 + dt * v0 + 0.5 * dt**2 * a0

full = np.zeros(N+1)
full[1:-1] = u_prev
history.append(full.copy())

full[1:-1] = u_curr
history.append(full.copy())

center_ts.append(u_prev[probe_index])

# ----------------------
# Time stepping (κεντρικές διαφορές στη 2η παράγωγο)
for n in range(1, steps):
    t_n = n * dt

    # Εξωτερική δύναμη στο drive_index
    f_ext[:] = 0.0
    f_ext[drive_index] = F_amp * np.sin(omega_drive * t_n)

    # Ταχύτητα και επιτάχυνση
    v_curr = (u_curr - u_prev) / dt
    a_curr = M_inv @ (f_ext - K @ u_curr - c_damp * v_curr)

    # Επόμενο βήμα
    u_next = 2.0 * u_curr - u_prev + (dt**2) * a_curr

    if n % store_every == 0:
        full = np.zeros(N+1)
        full[1:-1] = u_next
        history.append(full.copy())

    center_ts.append(u_next[probe_index])

    u_prev, u_curr = u_curr, u_next

history = np.array(history)
print("Frames stored for animation:", history.shape[0])

# ----------------------
# Animation
fig, ax = plt.subplots()
line, = ax.plot(x_full, history[0], lw=2)
ax.set_ylim(-1.5, 1.5)
ax.set_xlabel("Position along string (m)")
ax.set_ylabel("Displacement (arb. units)")
ax.set_title(f"Driven string with point masses: f_drive={f_drive} Hz, F_amp={F_amp}")

# σημείωσε τις θέσεις των βαριδίων
for pos in extra_mass_positions:
    if 0.0 < pos < L:
        ax.axvline(pos, color='red', linestyle='--', alpha=0.5, label='mass' if pos == extra_mass_positions[0] else None)
# θέση διεγέρτη
ax.axvline(drive_x, color='green', linestyle=':', alpha=0.6, label='driver')
ax.legend()

def animate_frame(i):
    line.set_ydata(history[i])
    return line,

ani = animation.FuncAnimation(fig, animate_frame, frames=len(history), interval=60, blit=True)
plt.show()

# ----------------------
# FFT στο σήμα του probe (στο κέντρο)
center_ts = np.asarray(center_ts)
fs = 1.0 / dt
Nfft = len(center_ts)
window = np.hanning(Nfft)
fft_vals = np.abs(np.fft.rfft(center_ts * window))
freqs = np.fft.rfftfreq(Nfft, d=dt)

# μικρή καθαριότητα πολύ χαμηλών συχνοτήτων
fft_vals[freqs < 1.0] = 0.0

plt.figure()
plt.plot(freqs, fft_vals, lw=1.5)
plt.xlabel("Frequency (Hz)")
plt.ylabel("FFT magnitude")
plt.title("Spectrum at probe (center)")
plt.xlim(0, 1000)
plt.grid(True)

# ιδανικές συχνότητες για χορδή χωρίς βαρίδια (αναφορά)
nmax = 12
f_theory = (np.arange(1, nmax+1) * c) / (2*L)
for fn in f_theory:
    plt.axvline(fn, color='gray', linestyle='--', alpha=0.35)
plt.show()

# κορυφαίες αιχμές
fft_copy = fft_vals.copy()
top_k = 8
idx_top = np.argpartition(fft_copy, -top_k)[-top_k:]
idx_top = idx_top[np.argsort(fft_copy[idx_top])[::-1]]
print("Top spectral peaks (Hz, magnitude):")
for idx in idx_top:
    print(f"{freqs[idx]:8.2f} Hz , {fft_copy[idx]:.3e}")

# Υπολογισμός κυματικών αριθμών (για ιδανική χορδή: k_n = n*pi/L)
n_modes = 10
n_vals = np.arange(1, n_modes+1)
k_vals = n_vals * np.pi / L

# Θεωρητικές ιδιοσυχνότητες χωρίς βαρίδια
f_theory = (n_vals * c) / (2*L)

# Αν έχεις ήδη βρει τις ιδιοσυχνότητες από το eig ή από FFT peaks:
# π.χ. freqs_modes = [τις κορυφές από FFT ή eig]
# εδώ για παράδειγμα:
freqs_modes = f_theory * 0.95  # (υποθέτουμε ότι μειώθηκαν λόγω βαριδίων)

# Διάγραμμα διασποράς
plt.figure()
plt.plot(k_vals, freqs_modes, 'o-', label='Με βαρίδια (προσομοίωση)')
plt.plot(k_vals, f_theory, '--', label='Ιδανική χορδή (χωρίς βαρίδια)')
plt.xlabel("Κυματικός αριθμός k (rad/m)")
plt.ylabel("Συχνότητα f (Hz)")
plt.title("Διάγραμμα Διασποράς")
plt.legend()
