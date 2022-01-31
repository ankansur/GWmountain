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

A = np.loadtxt('Atemp.txt')
rho = np.loadtxt('rhotemp.txt')/1e7
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
        r[i] = (0.0001+i*dr1)
    else:
        r[i] = (0.0001+rc + i*dr2)

thetas = np.linspace(data[1],data[2],len(A[0]))

print (thetas)
print (r)

th, rr = np.meshgrid(thetas,r)

#rho = rho.astype('float')
#rho[rho==0]='nan'

#np.savetxt('rho_new.txt',rho)

x, z = spherical_to_cartesian2D(r, thetas)

f = plt.figure(figsize=(14,14))
ax = f.add_subplot(111)
Bfieldl = ax.contour(x, z, A, 12, linewidths=0.7, linestyles='solid',colors='k') 
rhomag = ax.contourf(x, z, rho, levels=ls, cmap='RdPu')
cbar = plt.colorbar(rhomag)
ax.tick_params(axis="x", labelsize=18) 
ax.tick_params(axis="y", labelsize=18)
ax.set_xlabel(r'x (km) ',size=22)
ax.set_ylabel(r'z (km)',size=22)
ticklabs = cbar.ax.get_yticklabels()
cbar.set_label(r'$\rho/ \rm gm \,cm^{-3}) $',fontsize=20)
cbar.ax.set_yticklabels(ticklabs, fontsize=12)
cbar.ax.tick_params(labelsize=12)
ax.set_ylim(0,10.0)
ax.set_xlim(0,10)
plt.savefig('fieldlines.pdf',bbox_to_inches='tight')
plt.show()
