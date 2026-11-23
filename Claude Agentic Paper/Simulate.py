import numpy as np
import torch

c = 299792458 

# Create a reference pulse
def simulate_reference(L, deltat):
    toff = 1.0e-11
    twidth = 8.0e-13
    tdecay = 1.0e-12
    scale = 1.0e12

    t = torch.arange(0, L, dtype=torch.float32) * deltat - toff
    x = -scale * t * torch.exp(-(t / twidth) ** 2 - t / tdecay)

    return x


def rts_batched(n0, nj, Dj):
    # All inputs are shape [F]
    c = torch.cos(nj * Dj)
    s = torch.sin(nj * Dj)
    d = c + (0.5j) * (nj / n0 + n0 / nj) * s
    r = (0.5j) * s * (n0 / nj - nj / n0) / d
    t = 1.0 / d
    return r, t * torch.exp(1j * n0 * Dj), t * t - r * r




def RTm_batched(m, n0, layers):
    # Each entry in layers is a tuple of shape [F] tensors: (nj, Dj)
    F = layers[0][0].shape[0]

    U = torch.zeros(F, dtype=torch.cfloat, device=n0.device)
    V = torch.ones(F, dtype=torch.cfloat, device=n0.device)
    T = torch.ones(F, dtype=torch.cfloat, device=n0.device)

    for j in range(m):
        nj, Dj = layers[j]
        r, t, s = rts_batched(n0, nj, Dj)

        Vlast = V.clone()
        U_new = r * V + s * U
        V_new = V - r * U
        U, V = U_new, V_new

        T = T * t * Vlast / V

    R = U / V
    return R, T



# Simulate reference pulse through material, add some noise
def _compute_noise_std_from_snr(signal, snr_db):
    # signal: real-valued torch tensor
    power = torch.mean(signal.to(torch.float64) ** 2)
    noise_power = power / (10 ** (snr_db / 10.0))
    return torch.sqrt(noise_power).to(signal.dtype)

def _add_noise_time(signal, noise_level=None, snr_db=None, noise_mode='awgn', exponent=1.0, seed=None):
    if seed is not None:
        torch.manual_seed(seed)

    if snr_db is not None:
        noise_level = _compute_noise_std_from_snr(signal, snr_db)

    if noise_level is None:
        return signal  # no noise to add

    if noise_mode == 'awgn':
        noise = noise_level * torch.randn_like(signal)
        return signal + noise
    elif noise_mode in ('pink', '1/f'):
        N = signal.shape[0]
        freqs = torch.fft.rfftfreq(N, d=1.0)
        S_f = torch.where(freqs == 0, torch.tensor(0.0, device=signal.device), 1.0 / (freqs ** (exponent / 2.0)))
        phases = 2 * np.pi * torch.rand_like(S_f)
        noise_freq = S_f * torch.exp(1j * phases)
        noise_time = torch.fft.irfft(noise_freq, n=N)
        noise_time = noise_time / torch.std(noise_time) * noise_level
        return signal + noise_time
    else:
        raise ValueError(f"Unsupported noise_mode: {noise_mode}")
    
    
def simulate_reference_noiseafter(L, deltat, noise_level=None, twidth=8.0e-13):
    toff = 1.0e-11
    tdecay = 1.0e-12 

    t = torch.arange(0, L, dtype=torch.float32) * deltat - toff

    x = -t * torch.exp(-(t / twidth) ** 2 - t / tdecay)
    
    x /= x.abs().max()   
    if noise_level:
        x += noise_level * torch.randn(L, dtype=torch.float32)  

    return x


def noise(snr_dB):
    """SNR noise to amplitude noise for noise addition to plots."""
    amp_ratio = 10**(snr_dB / 20.0)
    noise_amplitude_fraction = 1.0 / amp_ratio  
    return noise_amplitude_fraction

def simulate_parallel(x, layers, deltat, noise_level=None, noise_mode='awgn', snr_db=None, noise_exponent=1.0, seed=None):
    #Defining parameters
    L = len(x)
    M = 2 * L
    N = 4 * L
    deltaf = 1.0 / (N * deltat)
    dk = 2 * torch.pi * deltaf / c  

    device = x.device

    # Layer parameters
    indices = torch.stack([
        l[0] if isinstance(l[0], torch.Tensor)
        else torch.tensor(l[0], dtype=torch.cfloat, device=device, requires_grad=True)
        for l in layers
    ])
    thicknesses = torch.stack([
        l[1] if isinstance(l[1], torch.Tensor)
        else torch.tensor(l[1], dtype=torch.cfloat, device=device, requires_grad=True)
        for l in layers
    ])
    m = len(layers)

    # Frequency indices
    k_vals = torch.arange(M + 1, dtype=torch.float32, device=device)
    kD = dk * k_vals[:, None] * thicknesses[None, :]  # Shape: [M+1, m]

    # Build per-frequency layer parameters
    batched_layers = [(indices[j].expand(M+1), kD[:, j]) for j in range(m)]
    
    # Compute batched transmission
    _, T_half = RTm_batched(m, torch.tensor(1.0, dtype=torch.cfloat, device=device), batched_layers)

    # Build full spectrum using conjugate symmetry
    T = torch.zeros(N, dtype=torch.cfloat, device=device)
    T[:M+1] = T_half
    T[M+1:] = torch.conj(torch.flip(T[1:M], dims=[0]))  # Avoid duplicating DC and Nyquist

    # FFT input
    z = torch.zeros(N, dtype=torch.float, device=device)
    z[:L] = x
    X = torch.fft.fft(z) / N

    # Apply transmission
    Y = T * X
    y = N * torch.fft.ifft(Y).real

    # Add measurement noise — several modes supported:
    # - 'awgn' (default): additive white Gaussian noise in time domain
    # - 'pink' or '1/f'  : 1/f^noise_exponent colored noise in time domain
    # You can specify either `noise_level` (std) or `snr_db` (signal-to-noise ratio in dB).
    if noise_level or snr_db is not None:
        # operate on real-valued signal y
        y = _add_noise_time(y, noise_level=noise_level, snr_db=snr_db, noise_mode=noise_mode, exponent=noise_exponent, seed=seed)

    return T, y


def simulate_from_signal(signal, layers, deltat, noise_level=None, noise_mode='awgn', snr_db=None, noise_exponent=1.0, seed=None, device=None):
    """
    Like `simulate_parallel`, but accepts an external THz time-domain signal (numpy array or torch tensor)
    as the input pulse. This is a convenience wrapper that converts numpy inputs to torch tensors
    and forwards to `simulate_parallel`.

    Args:
        signal: 1D numpy array or torch tensor containing the measured/reference time-domain pulse.
        layers: list of (n_complex, thickness) tuples (same as `simulate_parallel`).
        deltat: time step [s].
        noise_level, noise_mode, snr_db, noise_exponent, seed: same semantics as `simulate_parallel`.
        device: optional; torch device string or torch.device to place the input on (e.g. 'cpu' or 'cuda').

    Returns:
        (T, y) same as `simulate_parallel` (torch tensors).
    """
    # Determine target device
    if isinstance(signal, torch.Tensor):
        x = signal.to(dtype=torch.float32)
        target_device = x.device if device is None else torch.device(device) if isinstance(device, str) else device
        x = x.to(device=target_device)
    else:
        # assume numpy-like
        target_device = torch.device(device) if device is not None else torch.device('cpu')
        x = torch.as_tensor(signal, dtype=torch.float32, device=target_device)

    return simulate_parallel(x, layers, deltat, noise_level=noise_level, noise_mode=noise_mode, snr_db=snr_db, noise_exponent=noise_exponent, seed=seed)


# Loss landscape functions
def create_tmm_loss_function(x_exp, true_x_exp, deltat, L):
    def loss_function(params):
        n = params['n']
        k = params['k']
        d = params['d']
        layers = [(n - 1j*k, d)]
        
        T_sim, y_sim = simulate_from_signal(x_exp, layers, deltat)
        # Use the actual length of true_x_exp instead of L
        actual_len = len(true_x_exp)
        y_sim = y_sim[:actual_len]
        loss = torch.mean((y_sim - true_x_exp)**2).item()
        
        return loss
    return loss_function

#if __name__ == '__main__':
#    import time
#    import matplotlib.pyplot as plt


    # Simulation parameters
    L = 2**12
    deltat = (t_r[1]-t_r[0])*1e-12  # Time step from data

    
    # Generate reference pulse
    x = simulate_reference(L, deltat)

    layers = [(torch.tensor(2-0.01j, dtype=torch.cfloat), torch.tensor(1e-3, dtype=torch.cfloat))]
    
    # Time the function
    t0 = time.perf_counter()
    T, y = simulate_from_signal(x, layers, deltat)
    t1 = time.perf_counter()
    

    print(f"simulate_parallel execution time: {t1 - t0:.4f} seconds")

    y = y[:L].detach().cpu().numpy()
    plt.figure(figsize=(12,4))
    plt.title('Time Domain of THz Pulse through single layered sample.')
    plt.plot(x, label='Reference Pulse')
    plt.plot(y, label='Sample Pulse')
    plt.legend()
    plt.show()
