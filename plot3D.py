import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
fig = plt.figure()
ax = fig.gca(projection='3d')
#imports for Gauss fit
import pylab as plb
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp


DataRaw =  '/home/just/Documents/powerdensitiy'
da= np.genfromtxt(DataRaw,skip_header=2, usecols=(0, 1, 2, 3,4,5,6))

xdata= da[:,0]
ydata= da[:,1]
Pdata= da[:,3]
y_prime_data= da[:,5]

#ax.scatter3D(xdata, ydata, Pdata, c=Pdata, cmap='Greens')
#plt.show()
ax.scatter3D(xdata, y_prime_data, Pdata, c=Pdata, cmap='Greens')
plt.show()


