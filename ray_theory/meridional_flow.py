# calculate meridional flow model suggested by Jouve & Brun, 2007
# define function F
import numpy as np
import matplotlib.pyplot as plt

d = np.load('data/models_data.npz')
rr = d['rr']
ro = d['ro']
rsun = d['rsun']
#

z0 = rr/rsun
nth = 256
thmax = np.pi
thmin = 0
dth = (thmax - thmin)/nth
th = np.zeros(nth)
th[0] = dth/2
for n in range(1,nth):
    th[n] = th[n-1] + dth

x0 = -np.cos(th)

RR, TH = np.meshgrid(rr/rsun,th,indexing='ij')
X, Y = RR*np.cos(TH), RR*np.sin(TH)

z, x = np.meshgrid(z0,x0,indexing='ij')
ro2, tmp = np.meshgrid(ro,th,indexing='ij')

k1 = 1
phx = k1*z*(z-0.6)*(500*(z-0.8)**3 - 20*(z-0.8))*(3*x**2 - 1)
rovr = phx/z**2
vr = rovr*z

phz = k1*(z-0.6)*(500*(z-0.8)**3 - 20*(z-0.8))*(x**3 - x) \
    + k1*z*(500*(z-0.8)**3 - 20*(z-0.8))*(x**3 - x) \
    + k1*z*(z-0.6)*(1500*(z-0.8)**2 - 20)*(x**3 - x) \

rovt = -phz/z/np.sqrt(1-x**2)
vt = rovt*z

vr[RR < 0.6] = 0.
vt[RR < 0.6] = 0.

fac = 1350/vt.max()
vr = vr*fac
vt = vt*fac

np.savez('data/meridional_flow_jouve2007.npz',vr=vr,vt=vt,rr_meri=z0,th_meri=th)