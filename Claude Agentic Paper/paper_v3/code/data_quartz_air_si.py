"""Quartz (1 mm) / air gap / Si (0.4 mm) measurements, Menlo ScanControl exports.

Both files of a run share the same time axis (3486 points, dt = 33 fs, 115 ps window, 3000-waveform average).

    run2 : ref_10_min_2.txt            / quartz_air_silicon_10_min_2.txt    (2026-10-20)
    run3 : ref_3.txt                   / quartz_air_silicon_15_min_3.txt    (2026-10-21, same stack re-measured)

`load_run(tag)` returns (t_r_raw, A_r, A_s); the module-level arrays below are run2 (backward compatible).
"""
import numpy as np

RUNS = {
    'run2': ("./data/ref_10_min_2.txt", "./data/quartz_air_silicon_10_min_2.txt"),
    'run3': ("./data/ref_3.txt",        "./data/quartz_air_silicon_15_min_3.txt"),
}


def load_run(tag='run2'):
    ref_file, sample_file = RUNS[tag]
    data_ref = np.loadtxt(ref_file)
    data_sample = np.loadtxt(sample_file)
    assert np.allclose(data_ref[:, 0], data_sample[:, 0]), 'reference and sample time axes differ'
    return data_ref[:, 0], data_ref[:, 1], data_sample[:, 1]


t_r_raw, A_r, A_s = load_run('run2')
t_s_raw = t_r_raw

# Nominal thicknesses (starting values for the fit)
d1 = 1e-3      # quartz
d_air = 50e-6  # air gap
d2 = 0.4e-3    # silicon
