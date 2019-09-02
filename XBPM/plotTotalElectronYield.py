#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

# March 2019
# Paul Scherrer Institut, PSI
# David Marco Just
# david.just@psi.ch


import numpy as np
from os import listdir, chdir, path
from os.path import isfile, join
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, ticker
import matplotlib.pyplot as plt
import matplotlib.mlab as ml
import scipy
from scipy import integrate
import re



CRXO = True #Slect source type for Yield data. If True, source is X-Ray Attenuation Length from CRXO website if false it is from Evaluated Nuclear Data File Libary
pathToYield = '/home/just/Documents/PSI/XBPM/rawData/CRXO_AttenuationLengths/W.txt'
pathToFluxes = '/home/just/Documents/PSI/XBPM/rawData/cSAXS_Calcs/DetailedFluxScanSLS2.0SS_U14/'
title = path.dirname(pathToFluxes)
#title = path.basename(title)
title = 'SLS21_U15_K1.87_100-30000eV'

maxEnergy = 30000.0  # eV given by flux calculations
minEnergy = 30.0  # eV given by yield table
distanceFromSource = 10  # m distance from source at which the flux was calculated
bucketSize = 40 #  eV

def prepare_yield_data_CRXO(pathToYield):
    """Reads in the electron yield file, converts the attenuation length to cross section to mm²/g, removes unneeded vallues and returs a data set with (energy [eV]:cross section [mm²/g])"""
    yieldPerEnergy = []
    yieldPerEnergy = np.genfromtxt(pathToYield, skip_header=2, usecols=(0, 1))
    correctionFactor = 5
    #  set energy to eV
    yieldPerEnergy[:, 0] = yieldPerEnergy[:, 0]
    #  transpose Attenuation Length to photoinonaization crosssection
    yieldPerEnergy[:, 1] = 1 / yieldPerEnergy[:, 1] * correctionFactor
    # remove energies lower than min and higher than max
    i=0
    for E in yieldPerEnergy[:,0]:
        if E < minEnergy or E > maxEnergy:
            yieldPerEnergy = np.delete(yieldPerEnergy, i, 0)
            i=i-1
        i=i+1
    return(yieldPerEnergy)


def prepare_yield_data(pathToYield):
    """Reads in the electron yield file, converts energies to eV instead of MeV and the cross section to mm²/g instead of cm²/g, removes unneeded vallues and returs a data set with (energy [eV]:cross section [mm²/g]"""
    yieldPerEnergy = []
    yieldPerEnergy = np.genfromtxt(pathToYield, skip_header=18, usecols=(0, 8))

    #  set energy to eV
    yieldPerEnergy[:, 0] = yieldPerEnergy[:, 0] * 1000000.0
    #  set yield to mm²
    yieldPerEnergy[:, 1] = yieldPerEnergy[:, 1] / 100
    # remove energies lower than min and higher than max
    i=0
    for E in yieldPerEnergy[:,0]:
        if E < minEnergy or E > maxEnergy:
            yieldPerEnergy = np.delete(yieldPerEnergy, i, 0)
            i=i-1
        i=i+1
    return(yieldPerEnergy)

def sorted_aphanumeric(data):
    """Ensures that the data is read in from the lowest to the highest energy. Important for the integrate function"""
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(data, key=alphanum_key)

def read_flux_data(pathToFluxes, minEnergy, maxEnergy):
    """Read in the energy value
	read in coordinates and fluxes
    returns a data set with the format (Energy[eV],[x[mm],y[mm],flux[ph/s/mr^2/0.1%]]) all data in one set
    returns the number of energies in the data set i.e. hoe many [x[mm],y[mm],flux[ph/s/mr^2/0.1%]] arrays are present
    """
    chdir(pathToFluxes)
    da = []
    i = 0
    files = listdir(pathToFluxes)
    files = sorted_aphanumeric(files)
    for f in files:
        if isfile(join(pathToFluxes, f)) and f.endswith('.dta'):
            fo= open(f, "r")
            lines = list(fo)
            Eline = lines[5]
            fo.close()
            Energy=Eline[-11:-2]
            if float(Energy) >= minEnergy and float(Energy) <= maxEnergy:
                da.append(Energy)
                da.append(np.genfromtxt(f, skip_header=10, usecols=(0, 1, 2)))
                i=i+1
    noEnergies = i
    return(da, noEnergies)

def find_nearest(array, value):
    """Compares and array to a given value V and returns the closest value inside the array to V"""
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def photons_per_energy_bucket(fluxData,bucketSize):
    """Removes the bandwidth term by multipying the flux density with the bandwidth and dividing it by the bucket size. I.e. changes the fluy density form [ph/s/mr^2/0.1%] to [ph/s/mr^2/bucket] """
    for data in fluxData:
        if type(data)== str:
            Energy = float(data)
            BW = Energy * 0.001
        else:
            data[:,2] =data[:,2] * BW / bucketSize
    return fluxData

def multiply_flux_with_yield(fluxData, yieldPerEnergyData):
    """Returns yield or weighted flux densities"""
    yieldValue = 0
    for data in fluxData:
        if type(data)== str:
            Energy = float(data)
            nearestEnergy = find_nearest(yieldPerEnergyData[:, 0], Energy)
            for yieldPerEnergy in yieldPerEnergyData:
                if yieldPerEnergy[0]== nearestEnergy:
                    yieldValue= yieldPerEnergy[1]
        else:
            data[:,2]= data[:,2]*yieldValue
    return fluxData

def summ_all_weighted_fluxes(weightedFluxes):
    """summs up all yields per coordinate over the whole energy range"""
    i=0
    summedFluxes= weightedFluxes[1]
    for fluxPerEnergy in weightedFluxes:
        if i & 1 and i > 1: #only the arrays (odds)
            summedFluxes[:,2]= fluxPerEnergy[:,2]+summedFluxes[:,2]
        i=i+1
    return summedFluxes

def integrate_all_weigthed_fluxes(weightedFluxes):
    """integrates all yields per cooridinate over the whole energy range"""
    i=0
    energies = []
    demo = weightedFluxes[1]
    #get x-axis:
    for value in weightedFluxes:
        if i % 2 == 0: #only the energy vallues (evens)
            energies.append(float(value))
        i=i+1
    energies = np.asarray(energies)
    j=0 # itteration per coordinate
    len_j =len(demo)
    len_ii= len(energies)
    len_i = len(weightedFluxes)
    yarray = np.zeros((len(demo),len(energies)))
    while j < len(demo): #go throu all coordinates
        i=0 # itteration per energyVallue and arrays        
        ii=0 # itteration per energyVallue and arrays
        for value in weightedFluxes: #gehe durch alle energien
            if i & 1: #only the arrays (odds)
                yarray[j,ii] = value[j][2]
                ii=ii+1
            i=i+1
        j=j+1
    i=0
    for coordinate in yarray:
        #demo[i][2]= coordinate.sum()
        demo[i][2]=scipy.integrate.trapz(coordinate, energies)
        i=i+1
    return demo

def flux_per_mm_sqr(weightedFluxes, distanceFromSource):
    weightedFluxes[:,2]=weightedFluxes[:,2]/ distanceFromSource**2
    return weightedFluxes


def plot3D(x, y, z):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_trisurf(x, y, z, cmap=cm.jet, linewidth=0.2)
    ax.set_title(title)
    ax.set_xlabel('x, position hor. [mm]')
    ax.set_ylabel('y, position ver. [mm]')
    ax.set_zlabel("Flux, (arbitary)")
    plt.xlim(-5,5)
    plt.ylim(-5,5)
    plt.axis('scaled')
    plt.show()

def plot_top_view(x,y,z):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_trisurf(x, y, z, cmap=cm.jet, linewidth=0.2)
    ax.set_title(title)
    ax.set_xlabel('x, position hor. [mm]')
    ax.set_ylabel('y, position ver. [mm]')
    ax.set_zlabel("Flux, (arbitary)")
    plt.xlim(-5,5)
    plt.ylim(-5,5)
    ax.azim = -90
    ax.elev = 90
    plt.axis('scaled')
    plt.show()

def plot2D_bu(x, y, z): 
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
    ax = plt.contourf(xi, yi, zi, 15, colors = 'k')
    plt.pcolormesh(xi, yi, zi, cmap = plt.get_cmap('rainbow'))
    plt.colorbar() 
    #plt.scatter(x, y, marker = 'o', c = 'b', s = 5, zorder = 10)
    #plt.scatter(x, y, c = 'b', s = 5, zorder = 10)
    plt.xlim(-7, 7)
    plt.ylim(-7, 7)
    #plt.xlim(xmin, xmax)
    #plt.ylim(ymin, ymax)
    plt.show()
    #  see: https://stackoverflow.com/questions/13781025/matplotlib-contour-from-xyz-data-griddata-invalid-index

def plot2D(x, y, z): 
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
    plt.show()
    #  see: https://stackoverflow.com/questions/13781025/matplotlib-contour-from-xyz-data-griddata-invalid-index

def normalize(data):
    norm = (data-data.min())/(data.max()-data.min())
    return norm



if __name__ == '__main__':
    if CRXO == True:
        yieldPerEnergy = prepare_yield_data_CRXO(pathToYield)
    else:
        yieldPerEnergy = prepare_yield_data(pathToYield)
    fluxData, noEnegries = read_flux_data(pathToFluxes, minEnergy, maxEnergy)
    fluxData = photons_per_energy_bucket(fluxData, bucketSize)
    fluxDataYielded = multiply_flux_with_yield(fluxData,yieldPerEnergy)
    #allFluxes= summ_all_weighted_fluxes(fluxDataYielded)
    allFluxes= integrate_all_weigthed_fluxes(fluxDataYielded)
    allFluxes = flux_per_mm_sqr(allFluxes, distanceFromSource)
    plot3D(allFluxes[:,0], allFluxes[:,1], allFluxes[:,2])
    plot3D(allFluxes[:,0], allFluxes[:,1], normalize(allFluxes[:,2]))
    #plot_top_view(allFluxes[:, 0], allFluxes[:, 1], allFluxes[:, 2])
    #plot_top_view(allFluxes[:, 0], allFluxes[:, 1], normalize(allFluxes[:,2]))
    plot2D(allFluxes[:,0], allFluxes[:,1], allFluxes[:,2])
    plot2D(allFluxes[:,0], allFluxes[:,1], normalize(allFluxes[:,2]))

'''
x = allFluxes[:,0]
y = allFluxes[:,1]
z = allFluxes[:,2]
'''
