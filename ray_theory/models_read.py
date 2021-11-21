## This program requires sunmodel developed by H. Hotta
## If you need, please ask him.
import sunmodel

d = sunmodel.read()

rr = d['rr'] ## radial geometry [cm]
cs = d['cs'] ## speed of sound [cm/s]
rsun = d['rsun'] ## solar radius

np.savez('models_cs.npz',rr=rr,cs=cs,rsun=rsun)