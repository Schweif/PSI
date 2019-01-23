import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
#imports for Gauss fit
import pylab as plb
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp


fig = plt.figure()
ax = fig.gca(projection='3d')

DataRaw =  '/home/just/Documents/msShutterNew.dat'
da= np.genfromtxt(DataRaw,skip_header=0, usecols=(0, 1, 2))
#da= np.genfromtxt(DataRaw,skip_header=2, usecols=(0, 1, 2, 3,4,5,6))

xdata= da[:,0]
ydata= da[:,1]
Pdata= da[:,2]
y_prime_data= da[:,1]

ax.scatter3D(xdata, ydata, Pdata, c=Pdata, cmap='Greens')
#plt.show()
plt.xlabel('y, vertical position [mm]')
plt.ylabel('x, horizontal position [mm]')
#plt.zlabel("P, power density on shutter surface [W/mm^2]")
'''
surf = ax.plot_surface(xdata, ydata, Pdata, cmap=cm.coolwarm, linewidth=0, antialiased=False)


# Customize the z axis.
ax.set_zlim(0,11)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
'''
ax.scatter3D(xdata, ydata, Pdata, c=Pdata, cmap='Greens')
plt.show()



