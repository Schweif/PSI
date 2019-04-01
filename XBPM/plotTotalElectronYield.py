#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

# March 2019
# Paul Scherrer Institut, PSI
# David Marco Just
# david.just@psi.ch
#TODO: Import of FLux Tables schould be in numerical order in order to have the integration working

import numpy as np
from os import listdir, chdir
from os.path import isfile, join
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, ticker
import matplotlib.pyplot as plt
import matplotlib.mlab as ml
import scipy
from scipy import integrate




pathToYield = '/home/just/Documents/PSI/XBPM/rawData/EPDL97_74.dat'
pathToFluxes = '/home/just/Documents/PSI/XBPM/rawData/DiamondI18/'
maxEnergy = 30000.0  # eV given by flux calculations
minEnergy = 100.0  # eV given by yield table
distanceFromSource = 10  # m distance from source at which the flux was calculated
bucketSize = 40 #  eV

def prepare_yield_data(pathToYield):
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

def read_flux_data(pathToFluxes, minEnergy, maxEnergy):
    chdir(pathToFluxes)
    da = []
    i = 0
    for f in listdir(pathToFluxes):
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
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def photons_per_energy_bucket(fluxData,bucketSize):
    for data in fluxData:
        if type(data)== str:
            Energy = float(data)
            BW = Energy * 0.001
        else:
            data[:,2] =data[:,2] * BW / bucketSize
    return fluxData

def multiply_flux_with_yield(fluxData, yieldPerEnergyData):
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
    i=0
    summedFluxes= weightedFluxes[1]
    for fluxPerEnergy in weightedFluxes:
        if i & 1 and i > 1: #only the arrays (odds)
            summedFluxes[:,2]= fluxPerEnergy[:,2]+summedFluxes[:,2]
        i=i+1
    return summedFluxes

def integrate_all_weigthed_fluxes(weightedFluxes):
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
    while j < len(demo): #gehe durch alle coordinaten
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
    ax.set_xlabel('x, position hor. [mm]')
    ax.set_ylabel('y, position ver. [mm]')
    ax.set_zlabel("Flux, (arbitary)")
    plt.xlim(-5,5)
    plt.ylim(-5,5)
    ax.azim = -90
    ax.elev = 90
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
    plt.contourf(xi, yi, zi, 15, linewidths = 0.5, colors = 'k')
    plt.pcolormesh(xi, yi, zi, cmap = plt.get_cmap('rainbow'))

    plt.colorbar() 
    #plt.scatter(x, y, marker = 'o', c = 'b', s = 5, zorder = 10)
    plt.scatter(x, y, c = 'b', s = 5, zorder = 10)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.show()
    #  see: https://stackoverflow.com/questions/13781025/matplotlib-contour-from-xyz-data-griddata-invalid-index


if __name__ == '__main__':
    yieldPerEnergy = prepare_yield_data(pathToYield)
    fluxData, noEnegries = read_flux_data(pathToFluxes, minEnergy, maxEnergy)
    fluxData = photons_per_energy_bucket(fluxData, bucketSize)
    fluxDataYielded = multiply_flux_with_yield(fluxData,yieldPerEnergy)
    #allFluxes= summ_all_weighted_fluxes(fluxDataYielded)
    allFluxes= integrate_all_weigthed_fluxes(fluxDataYielded)
    allFluxes = flux_per_mm_sqr(allFluxes, distanceFromSource)
    plot3D(allFluxes[:,0], allFluxes[:,1], allFluxes[:,2])
    plot_top_view(allFluxes[:, 0], allFluxes[:, 1], allFluxes[:, 2])
    plot2D(allFluxes[:,0], allFluxes[:,1], allFluxes[:,2])

'''
x = allFluxes[:,0]
y = allFluxes[:,1]
z = allFluxes[:,2]
'''
