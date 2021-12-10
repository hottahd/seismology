import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pickle

d = np.load('data/ray_calc.npz')
delta_theta = d['delta_theta']
d = np.load('data/return_radius.npz')
rt = d['rt']
om = d['om']
d = np.load('data/models_data.npz')
rsun = d['rsun']

plt.clf()
plt.close('all')

fig1 = plt.figure('Distance-Depth',figsize=(8,6))
ax11 = fig1.add_subplot(111)

ax11.plot(delta_theta/np.pi*180,rt/rsun,color='black')
ax11.set_xlim(0,50)
ax11.set_ylim(0.68,1)

ax11.set_xlabel('$\Delta$ [degree]')
ax11.set_ylabel('$r_\mathrm{t}/R_\odot$')

fig1.tight_layout()
fig1.savefig('figs/distance_radius.pdf')

f = open('data/ray_calc_th_ray.npz','rb')
th_ray_list = pickle.load(f)
f.close()

f = open('data/ray_calc_rr_ray.npz','rb')
rr_ray_list = pickle.load(f)
f.close()

fig2 = plt.figure('Ray Path',figsize=(6,6))
ax21 = fig2.add_subplot(111,aspect='equal')

deltas = np.array([3.0,10.0,19.2,30.0,45])/180*np.pi

for delta in deltas:
    lt = np.argmin(abs(delta_theta - delta))
    rr_ray = rr_ray_list[lt]
    th_ray = th_ray_list[lt]

    shift = 15/180*np.pi
    x,y = rr_ray/rsun*np.cos(th_ray - th_ray[0] + shift), rr_ray/rsun*np.sin(th_ray - th_ray[0] + shift)
    ax21.plot(x, y,color='black',linestyle='--')

alpha = 0.5
c10 = patches.Circle(xy=(0, 0), radius=1.0, ec='black',fill=False)
c09 = patches.Circle(xy=(0, 0), radius=0.9, ec='black',fill=False,alpha=alpha)
c08 = patches.Circle(xy=(0, 0), radius=0.8, ec='black',fill=False,alpha=alpha)
c07 = patches.Circle(xy=(0, 0), radius=0.7, ec='black',fill=False,alpha=alpha)
ax21.add_patch(c10)
ax21.add_patch(c09)
ax21.add_patch(c08)
ax21.add_patch(c07)

ax21.set_xlim(0.25,1.15)
ax21.set_ylim(0.15,1.05)
ax21.tick_params(labelleft=False)
ax21.tick_params(labelbottom=False)

ax21.annotate(r'$\Delta=3.0^\circ$',xy=(0.96,0.31))
ax21.annotate(r'$\Delta=10.0^\circ$',xy=(0.92,0.42))
ax21.annotate(r'$\Delta=19.2^\circ$',xy=(0.84,0.56))
ax21.annotate(r'$\Delta=30.0^\circ$',xy=(0.72,0.71))
ax21.annotate(r'$\Delta=45^\circ$',xy=(0.5,0.89))

ax21.text(0.5,0.09,r'$r/R_\odot$=')
ax21.text(0.65,0.09,r'$0.7$')
ax21.text(0.75,0.09,r'$0.8$')
ax21.text(0.85,0.09,r'$0.9$')
ax21.text(0.95,0.09,r'$1.0$')

fig2.tight_layout()
fig2.savefig('figs/ray_plot.pdf')
plt.pause(0.1)
