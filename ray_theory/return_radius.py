## calculate return radius
import numpy as np
import matplotlib.pyplot as plt

d = np.load('data/models_data.npz')

rr   = d['rr']
cs   = d['cs']
omb2 = d['omb2']
oma  = d['oma']
rsun = d['rsun']

nx = len(rr)

om = 3.1e-3*2*np.pi # wave frequency (free parameter)

lt_max = 900

l = np.linspace(1,lt_max,lt_max) # spherical harmonic degree
rt0 = np.zeros(lt_max)         # 
rt1 = np.zeros(lt_max)         # 

a = rr**2/(cs**2/om**2)      #
b = rr**2/((om**2 - omb2)/(om**2 - oma**2)*cs**2/om**2)

for lt in range(1,lt_max+1):
    al = a/lt/(lt+1)
    bl = b/lt/(lt+1)
    
    for n in range(0,nx):
        if al[n] > 1:
            rt0[lt-1] = rr[n-1]
            break
    for n in range(0,nx):
        if bl[n] > 1:
            rt1[lt-1] = rr[n-1]
            break

rt = rt1.copy()
np.savez('data/return_radius.npz',om=om,l=l,rt=rt)

plt.close('all')
fig = plt.figure('return_radius',figsize=(6,6))
ax = fig.add_subplot(111)

ax.plot(l,rt/rsun,color='black')

ax.set_xscale('log')
ax.set_ylim(0,1)

ax.set_xlabel('Spherical harmonic degree $\ell$')
ax.set_ylabel('$r/R_\odot$')

fig.tight_layout()
fig.savefig('figs/return_radius.pdf')
plt.pause(0.1)

