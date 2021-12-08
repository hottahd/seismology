## This program requires sunmodel developed by H. Hotta
## If you need, please ask him.
import sunmodel
import numpy as np

d = sunmodel.read()

rr   = np.flip(d['rr']) ## radial geometry [cm]
cs   = np.flip(d['cs']) ## speed of sound [cm/s]
ro   = np.flip(d['ro']) ## density [g/cm^3]
pr   = np.flip(d['pr']) ## pressure [dyn/cm^2]
gm   = np.flip(d['gm']) ## ratio of heat capasity
gx   = np.flip(d['gx']) ## gravitational acceleration [cm/s^2]
omb2   = gx*np.flip(d['dlog_pr_dlog_r_gm_dlog_ro_dlog_r'])/rr ## square of Brunt-Vaisala frequency
rsun = d['rsun'] ## solar radius

nx = len(rr)

hpr = pr/ro/gx # pressure scale height
hro = 1/(gx/cs**2 + omb2/gx) # density scale height

oma = cs/2/hro # acoustic cutoff frequency

np.savez('data/models_data.npz',rr=rr,cs=cs,ro=ro,pr=pr,gm=gm,gx=gx,omb2=omb2,rsun=rsun,hpr=hpr,hro=hro,oma=oma)
