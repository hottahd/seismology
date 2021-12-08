import pickle
import numpy as np

d = np.load('data/models_data.npz')
rr = d['rr']
cs = d['cs']
omb2 = d['omb2']
oma = d['oma']
rsun = d['rsun']

nx = len(rr)

e = np.load('data/return_radius.npz')
om = e['om']
rt = e['rt']

#lt = 42 # spherical harmonic degree

delta_theta = np.zeros(len(rt))

th_ray_list = []
rr_ray_list = []

for lt in range(1,len(rt)+1):
    kh2 = lt*(lt+1)/rr**2
    kh = np.sqrt(kh2)
    kr2 = (om**2 - oma**2)/cs**2 + lt*(lt+1)/rr**2*(omb2 - om**2)/om**2
    ff = kh/rr*(1-omb2/om**2)
    gg = kr2
    rtt = rt[lt-1]
    nt = np.argmin(abs(rtt-rr))
    
    th_ray = np.zeros(nx-nt)
    rr_ray = np.zeros(nx-nt)
    for n in range(nt,nx):
        gg1 = gg[n]
        gg2 = gg[n+1]
        dgg = gg2 - gg1
        ff1 = ff[n]
        ff2 = ff[n+1]
        dff = ff2 - ff1
        drr = rr[n+1] - rr[n]
        ep = dgg/gg1

        if gg1 > 0 and gg2 > 0:
            if abs(ep) >= 0.01:
                I0 = gg1**(-0.5)*2/ep*((1+ep)**(0.5) - 1)
                I1 = gg1**(-0.5)*4/3/ep**2*(1-(1-ep/2)*(1+ep)**0.5)
            else:
                I0 = gg1**(-0.5)*(1-ep/4+ep**2/8)
                I1 = gg1**(-0.5)*0.5*(1-ep/3+3*ep**2/16)

        else:
            if gg1 < 0:
                # gg1 < 0 and gg2 > 0 
                I0 = 2*gg2**0.5/dgg
                I1 = 2*(gg2-3*gg1)*gg2**0.5/3/dgg**2
            else:
                # gg1 > 0 and gg2 < 0
                I0 = -2*gg1**0.5/dgg
                I1 = 4*gg1**1.5/3/dgg**2                        

        rr_ray[n-nt+1] = rr[n+1]
        if n == nt:
            th_ray[n-nt] = 0
            rr_ray[n-nt] = rr[n] - (rr[n+1]-rr[n])*gg1/dgg
            
        th_ray[n-nt+1] = th_ray[n-nt] + drr*(ff1*I0 + dff*I1)

        if gg1 > 0 and gg2 < 0:
            ns = n-nt+1
            break

    print(lt,th_ray.max()*2/np.pi*180,rtt/rsun)
    th_ray = th_ray[:ns+1]
    rr_ray = rr_ray[:ns+1]
    
    th_rayf = np.zeros(2*ns+1)
    rr_rayf = np.zeros(2*ns+1)

    th_rayf[:ns+1] = -np.flip(th_ray)
    rr_rayf[:ns+1] =  np.flip(rr_ray)

    th_rayf[ns+1:2*ns+1] = th_ray[1:ns+1]
    rr_rayf[ns+1:2*ns+1] = rr_ray[1:ns+1]

    th_ray_list.append(th_rayf)
    rr_ray_list.append(rr_rayf)
    delta_theta[lt-1] = th_ray.max()*2

np.savez('data/ray_calc.npz',delta_theta=delta_theta)
f = open('data/ray_calc_th_ray.npz','wb')
pickle.dump(th_ray_list,f)
f.close()

f = open('data/ray_calc_rr_ray.npz','wb')
pickle.dump(rr_ray_list,f)
f.close()

#plt.close('all')
#fig = plt.figure('frequency',figsize=(5,5))
#ax = fig.add_subplot(111,aspect='equal')

#x,y = rr_rayf/rsun*cos(th_rayf - th_rayf[0]), rr_rayf/rsun*sin(th_rayf - th_rayf[0])
#ax.plot(x, y,color='black',linestyle='--')
#c = patches.Circle(xy=(0, 0), radius=1, ec='black')
#ax.add_patch(c)

#ax11.plot(pr, omb2,color='black')
#ax11.plot(pr,-omb2,color='black',linestyle='--')

#ax11.plot(pr, oma**2,color=orange)

#ax11.set_xscale('log')
#ax11.set_yscale('log')

#np.savez('models_cs.npz',rr=rr,cs=cs,rsun=rsun)
