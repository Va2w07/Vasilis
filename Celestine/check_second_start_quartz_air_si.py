"""Second start for the quartz/air/Si fit: begin from the previous bilayer solution (quartz n=1.967, d=994 um;
Si n=3.435, d=382.7 um) with a 50 um gap, same two-stage protocol. Compares final loss / parameters with the
notebook run (start from nominal 1.982 / 1 mm, 3.40 / 0.4 mm)."""
import numpy as np, torch, sys
from dispersive_multilayer import Grid, DispersiveMultiLayer, fit
from data_quartz_air_si import t_r_raw, A_r, A_s

torch.set_default_dtype(torch.float32)
L = len(A_r); dt = (t_r_raw[1] - t_r_raw[0]) * 1e-12
grid = Grid(L, dt)
x = np.asarray(A_r, np.float32); y = np.asarray(A_s, np.float32)

tag = sys.argv[1] if len(sys.argv) > 1 else 'prev'
starts = {
    'prev': [dict(name='quartz', n=1.967, k=0.0043, d=994.0e-6),
             dict(name='air', n=1.0, k=0.0, d=50e-6, fit_n=False, fit_k=False),
             dict(name='Si', n=3.435, k=1e-4, d=382.7e-6)],
    'lit':  [dict(name='quartz', n=1.96, k=0.005, d=1.0e-3),
             dict(name='air', n=1.0, k=0.0, d=50e-6, fit_n=False, fit_k=False),
             dict(name='Si', n=3.418, k=1e-4, d=0.38e-3)],
}
m = DispersiveMultiLayer(grid, starts[tag])
m.set_dispersion(False)
tr1 = fit(m, x, y, n_iters=15000, lr=1e-4, print_step=2500, label=f'{tag} flat', progress=False)
m.report()
m.load_state(m.state()); m.set_dispersion(True)
tr2 = fit(m, x, y, n_iters=30000, lr=1e-4, print_step=5000, label=f'{tag} disp', progress=False)
m.report()
st = m.state()
np.savez(f'second_start_{tag}.npz', n=st['n'].numpy(), k=st['k'].numpy(), d=st['d'].numpy(),
         loss_flat=tr1['loss'], loss_disp=tr2['loss'], d_hist=np.concatenate([tr1['d'], tr2['d']]))
print('final loss', tr2['loss'][-1])
