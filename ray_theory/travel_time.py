# calculate travel time for given flow
import numpy as np
import matplotlib.pyplot as plt
import pickle

d = np.load('data/return_radius.npz')
rt = d['rt']
om = d['om']

d = np.load('data/ray_calc.npz')
delta_theta = d['delta_theta']

f = open('data/ray_calc_th_ray.npz','rb')
th_ray_list = pickle.load(f)
f.close()

f = open('data/ray_calc_rr_ray.npz','rb')
rr_ray_list = pickle.load(f)
f.close()

d = np.load('data/models_data.npz')
rsun = d['rsun']
rr = d['rr']
oma = d['oma']
cs = d['cs']
omb2 = d['omb2']


d = np.load('data/meridional_flow_jouve2007.npz')
vr = d['vr']
vt = d['vt']
rr_meri = d['rr_meri']
th_meri = d['th_meri']

mx = 64
rfa = np.linspace(0.5*rsun,0.99*rsun,mx)
tau = np.zeros(mx)

tf = 45/180*np.pi

for m in range(mx):    
    rf = rfa[m] # foucal radius
    n = np.argmin(abs(rt-rf))
    lt = n + 1

    kh2 = lt*(lt+1)/rr**2
    kh = np.sqrt(kh2)
    kr2 = (om**2 - oma**2)/cs**2 + lt*(lt+1)/rr**2*(omb2 - om**2)/om**2
    fft = kh*(1-omb2/om**2)/cs**2
    gg = kr2
    
    for i in range(0,len(rr_ray_list[n])-2):
        rrt = rr_ray_list[n][i]/rsun
        tht = th_ray_list[n][i] + tf

        tht1 = th_ray_list[n][i+1] + tf

        nr = np.argmin(abs(rrt-rr_meri))
        nt = np.argmin(abs(tht-th_meri))
        nt1 = np.argmin(abs(tht1-th_meri))

        gg1 = gg[nr]
        gg2 = gg[nr+1]
        dgg = gg2 - gg1
        ff1 = -2*fft[nr  ]*vt[nr  ,nt  ]
        ff2 = -2*fft[nr+1]*vt[nr+1,nt+1]
        dff = ff2 - ff1
        drr = rr[nr+1] - rr[nr]
        ep = dgg/gg1

        if gg1 > 0 and gg2 > 0:
            if abs(ep) >= 0.01:
                I0 = gg1**(-0.5)*2/ep*((1+ep)**(0.5) - 1)
                I1 = gg1**(-0.5)*4/3/ep**2*(1-(1-ep/2)*(1+ep)**0.5)
            else:
                I0 = gg1**(-0.5)*(1-ep/4+ep**2/8)
                I1 = gg1**(-0.5)*0.5*(1-ep/3+3*ep**2/16)
        else:
            if gg1 < 0 and gg2 > 0:
                # gg1 < 0 and gg2 > 0
                I0 = 2*gg2**0.5/dgg
                I1 = 2*(gg2-3*gg1)*gg2**0.5/3/dgg**2
            if gg1 > 0 and gg2 < 0:
                # gg1 > 0 and gg2 < 0
                I0 = -2*gg1**0.5/dgg
                I1 = 4*gg1**1.5/3/dgg**2
            if gg1 < 0 and gg2 < 0:
                I0 = 0
                I1 = 0
    
        tau[m] = tau[m] + drr*(ff1*I0 + dff*I1)
    print(m,tau[m])

plt.close('all')
fig = plt.figure('travel_time',figsize=(5,5))
ax = fig.add_subplot(111)

ax.plot(rfa/rsun,tau,color='red',linestyle='--')
ax.set_xlim(0.5,1)
ax.set_ylim(0,1.2)

ax.tick_params(labelleft=False,labelbottom=False)

plt.pause(0.1)
fig.savefig('figs/travel_time_jouve2007.pdf',transparent=True)
