from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
#imports for Gauss fit
import pylab as plb
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp

DataRaw =  '/home/just/Documents/msShutterNew.dat'
import numpy as np #did not work on first imort doing it twice instead
da= np.genfromtxt(DataRaw,skip_header=0, usecols=(0, 1, 2))

xdata= da[:,0].reshape((-1,101))
ydata= da[:,1].reshape((-1,101))
zdata= da[:,2].reshape((-1,101))

X=xdata
Y=ydata
Z=zdata

#fig = plt.figure(figsize=plt.figaspect(0.5)*1.5) 
fig = plt.figure()
ax = plt.axes(projection='3d')

#ax.set_aspect('equal')
surf = ax.plot_surface(xdata, ydata, zdata,linewidth=0)

#keep Aspect ratio

ax.auto_scale_xyz([-40, 40], [-135, 135], [0, 11])
ax.pbaspect = [0.6, 2.0, 0.8]


#Axis Definition
ax.set_xlabel('x, position hor. [mm]')
ax.set_ylabel('y, position ver. [mm]')
ax.set_zlabel("P, power dens. [W/mm^2]")


# Set viewpoint.
ax.azim = -75
ax.elev = 15


#ax.scatter3D(xdata, ydata, zdata, c=zdata, cmap='Greens')


plt.show()
