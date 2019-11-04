#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
#imports for 3D Ploting of collected data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.mlab as ml
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
fig = plt.figure()
ax = fig.gca(projection='3d')

#imports for Gauss fit
import pylab as plb
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp

csv=False

DataRaw =  '/home/just/Documents/PSI/rawData/06_U14_SLS1_K1.5_PowerDensity_10m.dta'
dist_from_source = 10 #m


def Gauss(x,a,x0,sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))

def FitAndPlot2DGauss(axis):
    if axis == 'x':
        sigma = sigmaX
        mean = meanX
        n = nX
        x= y0xdata
        y= y0zdata
    if axis == 'y':
        sigma = sigmaY
        mean = meanY
        n = nY
        x= x0ydata
        y= x0zdata
    if axis == 0:
        n= len(y)
        mean = sum(x * y) / sum(y)
        sigma = np.sqrt(sum(y * (x - mean)**2) / sum(y))

 
    popt,pcov = curve_fit(Gauss, x, y, p0=[max(y), mean, sigma])
    ampl  = popt[0]
    mean  = popt[1]
    sigmaN = popt[2]

    print "*********************************"
    print "* Results for "+axis +" - axis          *"
    print "* (x = horizontal, y = vertical)*"
    print "*********************************"
    print "sigma (est)= " +str(sigma) +" mm"
    print "sigma (fit)= " +str(sigmaN) +" mm"
    print "sigma = " +str(sigmaN/dist_from_source) +" mrad"
    print "3 sigma = " +str(3*sigmaN) +" mm"
    print "3 sigma = " +str(3*sigmaN/dist_from_source) +" mrad"
    print "mean = " +str(mean)
    print "************************************************"

    #print "pcov = " +str(pcov)
    plt.plot(x, y, 'b+:', label='data')
    plt.plot(x, Gauss(x, *popt), 'r-', label='fit')
    plt.legend()
    plt.title(axis+' Axis')
    plt.xlabel('length (mm)')
    plt.ylabel('PowerDens (W/mm**2)')
    plt.text(1,0,"sigma = " +str(round(sigmaN,3)) +" mm", fontsize=12, horizontalalignment='left', verticalalignment='bottom')
    plt.show()  

  
def plot2D(x, y, z,txt=''):
    if not 'title' in globals():
        title = ''
    if not 'autoSave' in globals():
        autoSave = False
    fname = title + '_2D' + txt + '.png' 
    xmax= np.max(x)
    ymax= np.max(y)
    xmin= np.min(x)
    ymin= np.min(y)
    nx = len(x)
    ny= len(y)
    
    #xi = np.linspace(-5, 5, nx)
    #yi = np.linspace(-5, 5, ny)
    xi = np.linspace(xmin, xmax, nx)
    yi = np.linspace(ymin, ymax, ny)
    zi = ml.griddata(x, y, z, xi, yi)
    CS = plt.contourf(xi, yi, zi, 15, colors = 'k')
    CS2 = plt.pcolormesh(xi, yi, zi, cmap = plt.get_cmap('rainbow'))
    plt.title(title)
    plt.xlabel('x, position hor. [mm]')
    plt.ylabel('y, position ver. [mm]')

    cbar = plt.colorbar(CS2) 
    cbar.ax.set_ylabel("Flux, (arbitary)")
    
    #plt.scatter(x, y, marker = 'o', c = 'b', s = 5, zorder = 10)
    #plt.scatter(x, y, c = 'b', s = 5, zorder = 10)
    plt.xlim(-7, 7)
    plt.ylim(-7, 7)
    #plt.xlim(xmin, xmax)
    #plt.ylim(ymin, ymax)
    if autoSave == True:
        plt.savefig(fname)
        plt.close()
        plt.clf()
    else:
        plt.show()
        plt.clf()
    #  see: https://stackoverflow.com/questions/13781025/matplotlib-contour-from-xyz-data-griddata-invalid-index

##Prepare Data
if csv== True:
    da= np.genfromtxt(DataRaw,delimiter=",",skip_header=1, usecols=(1, 2, 3))
else:
    da= np.genfromtxt(DataRaw,skip_header=10, usecols=(0, 1, 2))

xdata= da[:,0]
ydata= da[:,1]
zdata= da[:,2]


#cut throu x=0
x0zdata=[]
x0ydata=[]

for d in da:
    if d[0] == 0:
        x0ydata.append(d[1])
        x0zdata.append(d[2])

#cut throu y=0
y0zdata=[]
y0xdata=[]

for d in da:
    if d[1] == 0:
        y0xdata.append(d[0])
        y0zdata.append(d[2])

#get Numpy arrays from list
x0ydata=np.array(x0ydata)
x0zdata=np.array(x0zdata)
y0xdata=np.array(y0xdata)
y0zdata=np.array(y0zdata)

##Estimate starting vallues for fit
# for x=0
nY = len(x0zdata)
meanY = sum(x0ydata*x0zdata)/ sum(x0zdata)
sigmaY = np.sqrt(sum(x0zdata*(x0ydata-meanY)**2)/sum(x0zdata))

# for y=0
nX = len(y0zdata)
meanX = sum(y0xdata*y0zdata)/ sum(y0zdata)
sigmaX = np.sqrt(sum(y0zdata*(y0xdata-meanX)**2)/sum(y0zdata))


##Plot and Print
#3D Plot
ax = plt.axes(projection='3d')
ax.scatter3D(xdata, ydata, zdata, c=zdata, cmap='Greens')
plt.show()

#2D Colour Plot
plot2D(xdata, ydata, zdata)

#Cut Plots    
FitAndPlot2DGauss('x')
FitAndPlot2DGauss('y')
