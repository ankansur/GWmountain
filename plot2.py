import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as colors
from matplotlib import ticker, cm

def spherical_to_cartesian2D(r, theta):
    
    #Given 1D arrays for r and theta, the function makes a spherical (r,theta)
    #grid and then transforms it to cartesian coordinates. It outputs two 2D
    #arrays, one for x and one for z.
    
    theta_matrix, radius_matrix = np.meshgrid(theta, r)
    x = radius_matrix*np.sin(theta_matrix)
    z = radius_matrix*np.cos(theta_matrix)
    return x, z

A = np.loadtxt('A.txt')
rho = np.loadtxt('rho.txt')/1e7
data = np.genfromtxt('Source.txt',skip_header=1)


lmin =  np.amin(rho)
lmax = np.amax(rho)
ls = np.linspace(0,lmax,50)
Nr = len(A)

r = np.zeros(Nr)
Rc = 1e6
rc = data[3]
rout = data[0]
dr1 = rc/(Nr/2.0);
dr2 = (rout-rc)/(Nr/2.0 - 1.0);

for i in range(Nr):
    if (i<=Nr/2.0):
        r[i] = (i*dr1)
    else:
        r[i] = (rc + i*dr2)

thetas = np.linspace(data[1],data[2],len(A[0]))

th, rr = np.meshgrid(thetas,r)

rho = rho.astype('float')
rho[rho==0]='nan'

#np.savetxt('rho_new.txt',rho)


f = plt.figure(figsize=(12,8))
ax = f.add_subplot(111)
Bfieldl = ax.contour(th, rr, A, 12, linewidths=2.0, linestyles='solid',colors='k') 
rhomag = ax.contourf(th, rr, rho, levels=ls, cmap='jet')
cbar = plt.colorbar(rhomag)
ax.tick_params(axis="x", labelsize=18) 
ax.tick_params(axis="y", labelsize=18)
ax.set_xlabel(r'colatitude [deg]',size=22)
ax.set_ylabel(r'radius [cm]',size=22)
ticklabs = cbar.ax.get_yticklabels()
cbar.set_label(r'$\rho/10^{7} \rm gm \,cm^{-3}) $',fontsize=20)
cbar.ax.set_yticklabels(ticklabs, fontsize=12)
cbar.ax.tick_params(labelsize=12)
ax.set_ylim(1,2*rc)
ax.set_xlim(0,14.4)
plt.savefig('fieldlines.pdf',bbox_to_inches='tight')
plt.show()
