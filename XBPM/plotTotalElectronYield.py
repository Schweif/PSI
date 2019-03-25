#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

# March 2019
# Paul Scherrer Institut, PSI
# David Marco Just
# david.just@psi.ch

import numpy as np
from os import listdir, chdir
from os.path import isfile, join
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

pathToYield = '/home/just/Documents/PSI/XBPM/rawData/EPDL97_74.dat'
pathToFluxes = '/home/just/Documents/PSI/XBPM/rawData/X04S_flux/'
maxEnergy = 30000.0  # eV given by flux calculations
minEnergy = 100.0  # eV given by yield table


def prepare_yield_data(pathToYield):
    yieldPerEnergy = []
    yieldPerEnergy = np.genfromtxt(pathToYield, skip_header=18, usecols=(0, 8))

    #  set energy to eV
    yieldPerEnergy[:, 0] = yieldPerEnergy[:, 0] * 1000000.0
    #  set yield to mm²
    yieldPerEnergy[:, 1] = yieldPerEnergy[:, 1] * 100
    # remove energies lower than min and higher tha max
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
            data[:,1]= data[:,1]*yieldValue
    return fluxData

def summ_all_weightetd_fluxes(weightedFluxes):
    i=0
    summedFluxes= weightedFluxes[1]
    for fluxPerEnergy in weightedFluxes:
        if i & 1 and i > 1: #only the arrays (odds)
            summedFluxes[:,2]= fluxPerEnergy[:,2]+summedFluxes[:,2]
        i=i+1
    return summedFluxes

def plot3D(xdata, ydata, zdata):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    surf = ax.plot_surface(xdata, ydata, zdata, linewidth=0)
    ax.set_xlabel('x, position hor. [mm]')
    ax.set_ylabel('y, position ver. [mm]')
    ax.set_zlabel("Flux, (arbitary)")





if __name__ == '__main__':
    yieldPerEnergy = prepare_yield_data(pathToYield)
    fluxData, noEnegries = read_flux_data(pathToFluxes, minEnergy, maxEnergy)
    fluxDataYielded = multiply_flux_with_yield(fluxData,yieldPerEnergy)
    allFluxes= summ_all_weightetd_fluxes(fluxDataYielded)
    plot3D(allFluxes[:,0], allFluxes[:,1], allFluxes[:,2])
    plt.show()