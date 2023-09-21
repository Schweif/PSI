# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 15:45:27 2023

@author: just_d
"""

#!/usr/bin/python3
# -*- coding: utf-8 -*-
#imports for 3D Ploting of collected data
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
import os
from pathlib import Path, PureWindowsPath
from termcolor import colored
import json


dataFormat= json
csv=False

DataRaw = PureWindowsPath(r"H:\00_final files\IDFE\X11MFE\5_SLS2.0\02_Berechnungen\02_XBPM\PowerDensity\01_LH\PowerDensity-10.txt")
dist_from_source = 11.492 #m
if dist_from_source == 1:
    xUnit= 'mrad'
else:
    xUnit= 'mm'
    
yLabel = "Power Density [W/mm^2]" # "Flux, (arbitary)" ; "Power Density [kW/mrad^2]" , "Flux Density [ph/s/mrad^2/0.1%BW]"
autoSave = True # Set Ture to automatically Save the Plots in DataRaw
txt=''

os.chdir(os.path.dirname(DataRaw))
#title= os.path.basename(DataRaw)[:-4]
title = 'UE36kn LH2 (Ptot=7.5 kW) @11.5 m'
maxZ=80 #Defines the maximum level of the Z axis. I.e. higher vallues than maxZ will be read in the plot, set to 0 if not used


### CONFIGURATION ENDS HERE




def plot2D(x, y, z, txt='', unit='', zMax=0):
    fname = title + txt + '.png' 
    xmax= np.max(x)
    ymax= np.max(y)
    xmin= np.min(x)
    ymin= np.min(y)
    
    X, Y = np.meshgrid(x, y)
    #Rearrange Data for Image Plot
    N = int(len(z)**.5)
    Z = z.reshape(N, N)
    
    fig, ax = plt.subplots()
    
    #This is used for isobars
    levels = np.arange(50, 90, 10) #start, end, intervall
    CS = ax.contour(Z, levels,colors='k', origin='lower', extent=[xmin, xmax, ymin, ymax])
    ax.clabel(CS, fontsize=9, inline=True,fmt='%1.0f') #fmt is used to edit the digits of the isobars label

    #this plots the back ground
    if zMax != 0:
        im = ax.imshow(Z,
                       cmap=cm.rainbow, 
                       interpolation=  'bilinear', #'none',
                       origin='lower', extent=[xmin, xmax, ymin, ymax],
                       vmax=zMax, vmin=z.min())
    else:
        im = ax.imshow(Z,
                       cmap=cm.rainbow, 
                       interpolation=  'bilinear', #'none',
                       origin='lower', extent=[xmin, xmax, ymin, ymax],
                       vmax=z.max(), vmin=z.min())                   
    
    plt.title(fname,pad=25)
    plt.xlabel('x, position hor. ' +unit)
    plt.ylabel('y, position ver. ' +unit)
    cbar = plt.colorbar(im, ax=ax)
    cbar.ax.set_ylabel(yLabel)
    if autoSave == True:
        plt.savefig(fname, dpi = 1200)
        plt.close()
        plt.clf()
    else:
        plt.show()
        plt.clf()
 
##Prepare Data


 
if maxZ != 0:
     print (colored("Warning: maxZ is not 0. Plots will be croped at " +str(maxZ) +"W/mm^2." , "yellow"))
 
#read into numpy format
if str(DataRaw).endswith('.csv'):
    da= np.genfromtxt(DataRaw,delimiter=",",skip_header=1, usecols=(1, 2, 3))
else:    
    da= np.genfromtxt(DataRaw,skip_header=2, usecols=(0, 1, 2))
    #New spectra Output has two header lines old version had 10
""" JSON import not finished    
if str(DataRaw).endswith('.json'):
    with open(DataRaw) as f:
        data = json.load(f)
        da=data['Output']['data']        
"""
plot2D(da[0:,0], da[0:,1], da[0:,2]*1000, txt= '', unit= '[' +xUnit + ']')

if maxZ != 0:
    plot2D(da[0:,0], da[0:,1], da[0:,2]*1000, txt= ' crptd@' +str(maxZ) +' W(mm)^-2', unit= xUnit, zMax= maxZ)
     
