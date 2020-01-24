#!/usr/bin/python3
# -*- coding: utf-8 -*-
#imports for 3D Ploting of collected data
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
import os

#imports for Gauss fit
from scipy.optimize import curve_fit
from scipy import exp

csv=False

DataRaw = "H:/final files/IDFE/X06SFE/5_SLS2.0/03_Calculations/02_Python/01_Flux_related/SpaceDept_FluxDens_0.1_0.1_1m_k0.4995_3499eV.dta"
dist_from_source = 1 #m
if dist_from_source == 1:
    xUnit= 'mrad'
else:
    xUnit= 'mm'
    
yLabel = "Flux Density [ph/s/mrad^2/0.1%BW]" # "Flux, (arbitary)" ; "Power Density [kW/mrad^2]" , "Flux Density [ph/s/mrad^2/0.1%BW]"

autoSave = False # Set Ture to automatically Save the Plots in DataRaw
txt=''

dir = os.path.dirname(DataRaw)
title= os.path.basename(DataRaw)[:-4]


def Gauss(x,a,x0,sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))

def FitAndPlot2DGauss(axis):
    if axis == 'x':
        fname = title + '_x' + txt + '.png' 
        sigma = sigmaX
        mean = meanX
        x= y0xdata
        y= y0zdata
        j=1
    if axis == 'y':
        fname = title + '_y' + txt + '.png' 
        sigma = sigmaY
        mean = meanY
        x= x0ydata
        y= x0zdata
        j=2
    if axis == 0:
        mean = sum(x * y) / sum(y)
        sigma = np.sqrt(sum(y * (x - mean)**2) / sum(y))

 
    popt,pcov = curve_fit(Gauss, x, y, p0=[max(y), mean, sigma])
    ampl  = popt[0]
    mean  = popt[1]
    sigmaN = popt[2]

    print( "*********************************")
    print( "* Results for "+axis +" - axis          *")
    print( "* (x = horizontal, y = vertical)*")
    print( "*********************************")
    print( "sigma (est)= " +str(sigma) +" mm")
    print( "sigma (fit)= " +str(sigmaN) +" mm")
    print( "sigma = " +str(sigmaN/dist_from_source) +" mrad")
    print( "3 sigma = " +str(3*sigmaN) +" mm")
    print( "3 sigma = " +str(3*sigmaN/dist_from_source) +" mrad")
    print( "mean = " +str(mean))
    print( "************************************************")

    #print( "pcov = " +str(pcov))
    plt.figure(j)
    plt.plot(x, y, 'b+:', label='data')
    plt.plot(x, Gauss(x, *popt), 'r-', label='fit')
    plt.legend()
    plt.title(axis+' Axis')
    plt.xlabel('length [' +xUnit+']')
    plt.ylabel(yLabel)
    plt.text(x.max(),ampl/15,"sigma = " +str(round(sigmaN,3)) +" " +xUnit, fontsize=9, horizontalalignment='right', verticalalignment='center')
    if autoSave == True:
        plt.savefig(dir+'/'+fname)
        plt.close()
        plt.clf()
    else:
        plt.show()


  
def plot2D(x, y, z,txt=''):
    fname = title + '_2D' + txt + '.png' 
    xmax= np.max(x)
    ymax= np.max(y)
    xmin= np.min(x)
    ymin= np.min(y)
    X, Y = np.meshgrid(x, y)

    N = int(len(z)**.5)
    Z = z.reshape(N, N)


    fig, ax = plt.subplots()
    im = ax.imshow(Z,
                   #norm=colors.LogNorm(vmin=Z.min(), vmax=Z.max()),
                   cmap=cm.rainbow, 
                   interpolation='bilinear',
                   #interpolation='none',
                   origin='lower', extent=[xmin, xmax, ymin, ymax],
                   vmax=z.max(), vmin=z.min())
    
    plt.title(fname,pad=25)
    plt.xlabel('x, position hor. ['+xUnit+']')
    plt.ylabel('y, position ver. ['+xUnit+']')
    cbar = plt.colorbar(im, ax=ax)
    cbar.ax.set_ylabel(yLabel)
    
    if autoSave == True:
        plt.savefig(dir+'/'+fname)
        plt.close()
        plt.clf()
    else:
        plt.show()


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
plt.figure(0)
fname = title + '_3D' + txt + '.png' 
ax = plt.axes(projection='3d')
ax.scatter3D(xdata, ydata, zdata, c=zdata, cmap='rainbow')
if autoSave == True:
    plt.savefig(dir+'/'+fname)
    plt.close()
    plt.clf()
else:
    plt.show()



#Cut Plots    
FitAndPlot2DGauss('x')
FitAndPlot2DGauss('y')

#2D Colour Plot
plot2D(xdata, ydata, zdata, txt)

