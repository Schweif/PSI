#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

# March 2019
# Paul Scherrer Institut, PSI
# David Marco Just
# david.just@psi.ch


import numpy as np
from os import listdir, chdir
from os.path import isfile, join
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, ticker
import matplotlib.pyplot as plt
import matplotlib.mlab as ml
import scipy
from scipy import integrate
import re
from decimal import *



DataA = '/home/just/Documents/PSI/XBPM/rawData/cSAXS_Calcs/PowerDensU14At10mK1.65.dta'
DataB = '/home/just/Documents/PSI/XBPM/rawData/cSAXS_Calcs/PowerDensU19At10mK1.78.dta'

dist_from_source = 10 #m
import numpy as np #did not work on first imort doing it twice instead
da= np.genfromtxt(DataA,skip_header=10, usecols=(0, 1, 2))
db= np.genfromtxt(DataB,skip_header=10, usecols=(0, 1, 2))

xdata= da[:,0]
ydata= da[:,1]
za= da[:,2]
zb= db[:,2]

z= za-zb



def eformat(f, prec, exp_digits):
    s = "%.*e"%(prec, f)
    mantissa, exp = s.split('e')
    # add 1 to digits as 1 is taken by sign +/-
    return "%se%+0*d"%(mantissa, exp_digits+1, int(exp))


def plot3D(x, y, z):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_trisurf(x, y, z,  cmap=cm.jet, linewidth=0.2)
    ax.set_xlabel('x, position hor. [mm]')
    ax.set_ylabel('y, position ver. [mm]')
    ax.set_zlabel("Power Density, W/mm^2")
    #plt.xlim(-5,5)
    #plt.ylim(-5,5)
    plt.axis('scaled')
    plt.show()

def plot2D(x, y, z): 
    xmax= np.max(x)
    ymax= np.max(y)
    xmin= np.min(x)
    ymin= np.min(y)
    nx = len(x)
    ny= len(y)
    
    xi = np.linspace(xmin, xmax, nx)
    yi = np.linspace(ymin, ymax, ny)
    zi = ml.griddata(x, y, z, xi, yi)
    plt.contourf(xi, yi, zi, 15, colors = 'k')
    plt.pcolormesh(xi, yi, zi, cmap = plt.get_cmap('rainbow'))

    plt.colorbar() 
    #plt.scatter(x, y, marker = 'o', c = 'b', s = 5, zorder = 10)
    #plt.scatter(x, y, c = 'b', s = 5, zorder = 10)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.show()
    #  see: https://stackoverflow.com/questions/13781025/matplotlib-contour-from-xyz-data-griddata-invalid-index

def plot2D_colourMesh(x, y, z): 
    xmax= np.max(x)
    ymax= np.max(y)
    xmin= np.min(x)
    ymin= np.min(y)
    nx = len(x)
    ny= len(y)
    
    xi = np.linspace(xmin, xmax, nx)
    yi = np.linspace(ymin, ymax, ny)
    zi = ml.griddata(x, y, z, xi, yi)
    #plt.contourf(xi, yi, zi, 15, colors = 'k')
    plt.pcolormesh(xi, yi, zi, cmap = plt.get_cmap('rainbow'))

    plt.colorbar() 
    #plt.scatter(x, y, marker = 'o', c = 'b', s = 5, zorder = 10)
    #plt.scatter(x, y, c = 'b', s = 5, zorder = 10)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.show()
    #  see: https://stackoverflow.com/questions/13781025/matplotlib-contour-from-xyz-data-griddata-invalid-index


i=0
while i <len(xdata):
    print eformat(xdata[i], 5, 2) +' ' +eformat(ydata[i], 5, 2) +' ' +eformat(z[i], 3, 2)
    i = i+1

plot3D(xdata,ydata,za*1000/(dist_from_source**2))
plot3D(xdata,ydata,zb*1000/(dist_from_source**2))
plot3D(xdata,ydata,z*1000/(dist_from_source**2))

print 'Data SetA Max Vallue:  ' +str(za.max()) +' kW/mrad^2, coresponds to: ' +str(za.max()*1000/(dist_from_source**2)) + ' W/mm^2'
print 'Data SetB Max Vallue:  ' +str(zb.max()) +' kW/mrad^2, coresponds to: ' +str(zb.max()*1000/(dist_from_source**2)) + ' W/mm^2'
print 'Difference Max Vallue: ' +str(z.max()) +' kW/mrad^2, coresponds to: ' +str(z.max()*1000/(dist_from_source**2)) + ' W/mm^2'

#plot2D(xdata,ydata,za*1000/(dist_from_source**2))
#plot2D(xdata,ydata,zb*1000/(dist_from_source**2))
plot2D(xdata,ydata,z*1000/(dist_from_source**2))

print 'done'
