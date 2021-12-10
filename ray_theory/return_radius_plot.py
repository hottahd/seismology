import numpy as np
import matplotlib.pyplot as plt
d = np.load('data/return_radius.npz')
l = d['l']
rt = d['rt']
om = d['om']

d = np.load('data/models_data.npz')

rsun = d['rsun']

plt.close('all')
fig = plt.figure('return_radius',figsize=(6,5))
ax = fig.add_subplot(111)

ax.plot(l,rt/rsun,color='black')
ax.set_xscale('log')
ax.set_xlim(1,1000)
ax.set_ylim(0,1)

ax.set_xlabel(r'Spherical harmonic degree $\ell$')
ax.set_ylabel(r'$r/R_\odot$')

fig.tight_layout()
plt.pause(0.1)
fig.savefig('figs/return_radius.pdf')
