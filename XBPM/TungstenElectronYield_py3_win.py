#!/usr/bin/python3.7
# -*- coding: utf-8 -*-

# Nov 2021
# Paul Scherrer Institut, PSI
# David Marco Just
# david.just@psi.ch


import numpy as np
from os import listdir, chdir, path
from os.path import isfile, join
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, ticker
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.mlab as ml
import scipy
from scipy import integrate
import re
import pandas as pd
import json


#CONFIGURATION
CRXO = True #Select source type for Yield data. If True, source is X-Ray Attenuation Length from CRXO website if false it is from Evaluated Nuclear Data File Libary
if CRXO == True:
    pathToYield = r"U:\05_InsertionDevices\BerchnungenFürXBPMs\Attenuations_e-_crosssections\W_Stiched.txt"
else:
    pathToYield = r"U:\05_InsertionDevices\BerchnungenFürXBPMs\Attenuations_e-_crosssections\EPDL97_74.dat"
pathToFluxes = r"E:\Spectra\02_U14_SLS2.0\01_RawData"
#title = path.dirname(pathToFluxes)+'SiC'
#title = path.basename(title)'
title = 'U14 @ SLS 2.0 (K1.46,N120) W_stiched'
autoSave = True # Set True to automatically save the plots in pathToFluxes
yLabel = "Tungsten response (arbitrary)"
maxEnergy = 30000.0  # eV given by flux calculations
minEnergy = 30.0  # eV given by yield table
distanceFromSource = 10  # m distance from source at which the flux was calculated
bucketSize = 1 #  eV
mrad = True #Sets X and Y axis in Plots to mrad. Set to False if X and Y axis in plots should be mm at distanceFromSource
#CONFIGURATION ENDS HERE


def prepare_yield_data_CRXO(pathToYield):
    """Reads in the electron yield file, converts the attenuation length to cross section to mm²/g, removes unneeded vallues and returns a data set with (energy [eV]:cross section [mm²/g])"""
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
    if np.isnan(yieldPerEnergy[0]).any():
        yieldPerEnergy= np.delete(yieldPerEnergy,0,0)
    return(yieldPerEnergy)


def prepare_yield_data(pathToYield):
    """Reads in the electron yield file, converts energies to eV instead of MeV and the cross section to mm²/g instead of cm²/g, removes unneeded vallues and returs a data set with (energy [eV]:cross section [mm²/g]"""
    yieldPerEnergy = []
    yieldPerEnergy = np.genfromtxt(pathToYield, skip_header=18, usecols=(0, 9))

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
    if np.isnan(yieldPerEnergy[0]).any():
        yieldPerEnergy= np.delete(yieldPerEnergy,0,0)
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
    returns the number of energies in the data set i.e. how many [x[mm],y[mm],flux[ph/s/mr^2/0.1%]] arrays are present
    """
    chdir(pathToFluxes)
    da = []
    i = 0
    files = listdir(pathToFluxes)
    files = sorted_aphanumeric(files)
    IgnoredFiles = []
    for f in files:
        if isfile(join(pathToFluxes, f)) and f.endswith(".dta"): 
            """prepares old Spectra data set"""
            print(f)
            fo= open(f, "r")
            lines = list(fo)
            Eline = lines[5]
            fo.close()
            Energy=Eline[-11:-1]
        
            if float(Energy) >= minEnergy and float(Energy) <= maxEnergy:
                da.append(Energy)
                da.append(np.genfromtxt(f, skip_header=10, usecols=(0, 1, 2)))
                i=i+1
        if  isfile(join(pathToFluxes, f)) and f.endswith(".json"):
            """prepares Spectra V11 and grater data sets"""
            with open(f) as v:
                data = json.load(v)
                if 'Output' in data:
                    print(f)
                    Energy = data['Output']['Set Value']
                    if Energy >= minEnergy and Energy <= maxEnergy:
                        da.append(str(Energy))
                        j=0
                        arr = np.empty((0,3))
                        for x in data['Output']['data'][0]:
                            for y in  data['Output']['data'][1]:
                                fluxDens = data['Output']['data'][2][j]
                                j=j+1
                                arr = np.append(arr, np.array([[x,y,fluxDens]]), axis=0)
                        da.append(arr)
                        i=i+1  
                    else:
                        IgnoredFiles.append(f)                                             
                else:
                    IgnoredFiles.append(f)              
        else:
            IgnoredFiles.append(f)       
    noEnergies = i
    print('Ignored files: ',IgnoredFiles)
    return(da, noEnergies)


def find_nearest(array, value):
    """Compares and array to a given value V and returns the closest value inside the array to V"""
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def photons_per_energy_bucket(fluxData,bucketSize):
    """Removes the bandwidth term by multipying the flux density with the bandwidth and dividing it by the bucket size. I.e. changes the fluxS density form [ph/s/mr^2/0.1%] to [ph/s/mr^2/bucket] """
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


def plot3D(x, y, z,txt=''):
    xmax= np.max(x)
    ymax= np.max(y)
    xmin= np.min(x)
    ymin= np.min(y)
    fname = title + '_3D' + txt + '.png'
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_trisurf(x, y, z, cmap=cm.jet, linewidth=0.2)
    ax.set_title(title)
    ax.set_xlabel('x, position hor. [mm]')
    ax.set_ylabel('y, position ver. [mm]')
    ax.set_zlabel("Flux, (arbitary)")
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)    
    #plt.xlim(-0.5,0.5)
    #plt.ylim(-0.5,0.5)
    #plt.axis('scaled')
    if autoSave == True:
        plt.savefig(fname)
        plt.clf()
    else:
        plt.show()
        plt.clf()


def plot_top_view(x,y,z,txt=''):
    xmax= np.max(x)
    ymax= np.max(y)
    xmin= np.min(x)
    ymin= np.min(y)
    fname = title + '_3D_top' + txt + '.png'
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_trisurf(x, y, z, cmap=cm.jet, linewidth=0.2)
    ax.set_title(title)
    ax.set_xlabel('x, position hor. [mm]')
    ax.set_ylabel('y, position ver. [mm]')
    ax.set_zlabel("Flux, (arbitary)")
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)
    #plt.xlim(-0.5,0.5)
    #plt.ylim(-0.5,0.5)
    ax.azim = -90
    ax.elev = 90
    plt.axis('scaled')
    if autoSave == True:
        plt.savefig(fname)
        plt.close()
        plt.clf()
    else:
        plt.show()
        plt.clf()


def plot2D_bu(x, y, z,txt=''):
    """Depercated Function"""
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
    ax = plt.contourf(xi, yi, zi, 15, colors = 'k')
    plt.pcolormesh(xi, yi, zi, cmap = plt.get_cmap('rainbow'))
    plt.colorbar() 
    #plt.scatter(x, y, marker = 'o', c = 'b', s = 5, zorder = 10)
    #plt.scatter(x, y, c = 'b', s = 5, zorder = 10)
    #plt.xlim(-0.5,0.5)
    #plt.ylim(-0.5, 0.5)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    if autoSave == True:
        plt.savefig(fname)
        plt.close()
        plt.clf()
    else:
        plt.show()
        plt.clf()
    #  see: https://stackoverflow.com/questions/13781025/matplotlib-contour-from-xyz-data-griddata-invalid-index


def plot2D(x, y, z, txt='', unit=''):
    fname = title + '_2D' + txt + '.png' 
    xmax= np.max(x)
    ymax= np.max(y)
    xmin= np.min(x)
    ymin= np.min(y)
    nx = len(x)
    ny= len(y)
    X, Y = np.meshgrid(x, y)
    #Rearrange Data for Image Plot
    N = int(len(z)**.5)
    Z = z.reshape(N, N)

    fig, ax = plt.subplots()
    im = ax.imshow(Z,
                   cmap=cm.rainbow, 
                   interpolation=  'bilinear', #'none',
                   origin='lower', extent=[xmin, xmax, ymin, ymax],
                   vmax=z.max(), vmin=-z.min())
    
    plt.title(fname,pad=25)
    plt.xlabel('x, position hor. ' +unit)
    plt.ylabel('y, position ver. ' +unit)
    cbar = plt.colorbar(im, ax=ax)
    cbar.ax.set_ylabel(yLabel)
    if autoSave == True:
        plt.savefig(fname)
        plt.close()
        plt.clf()
    else:
        plt.show()
        plt.clf()
 
    
def plot2D_Log(x, y, z, txt='', unit=''):
    fname = title + '_2D' + txt + '.png' 
    xmax= np.max(x)
    ymax= np.max(y)
    xmin= np.min(x)
    ymin= np.min(y)
    nx = len(x)
    ny= len(y)
    X, Y = np.meshgrid(x, y)

    N = int(len(z)**.5)
    Z = z.reshape(N, N)


    fig, ax = plt.subplots()
    im = ax.imshow(Z, 
                   norm=colors.LogNorm(vmin=Z.min(), vmax=Z.max()),
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
        plt.savefig(fname)
        plt.close()
        plt.clf()
    else:
        plt.show()
        plt.clf()
        
        
def normalize(data):
    norm = (data-data.min())/(data.max()-data.min())
    return norm

def saveToCSV(x, y, z, unit='mm'):
    d= {'x horizontal '+unit: x, 'y horizontal '+unit: y, 'z weighted flux [arb.]': z}
    df = pd.DataFrame(data=d)
    df.to_csv(title + '.csv')



if __name__ == '__main__':
    if CRXO == True:
        yieldPerEnergy = prepare_yield_data_CRXO(pathToYield)
    else:
        yieldPerEnergy = prepare_yield_data(pathToYield)
    fluxData, noEnergies = read_flux_data(pathToFluxes, minEnergy, maxEnergy)
    print( str(noEnergies) +' energy data sets read in.')
    fluxData = photons_per_energy_bucket(fluxData, bucketSize)
    fluxDataYielded = multiply_flux_with_yield(fluxData,yieldPerEnergy)
    #allFluxes= summ_all_weighted_fluxes(fluxDataYielded)
    allFluxes= integrate_all_weigthed_fluxes(fluxDataYielded)
    allFluxes = flux_per_mm_sqr(allFluxes, distanceFromSource)
    #saveToCSV(allFluxes[:,0], allFluxes[:,1], allFluxes[:,2])
    #plot3D(allFluxes[:,0], allFluxes[:,1], allFluxes[:,2])
    #plot3D(allFluxes[:,0], allFluxes[:,1], normalize(allFluxes[:,2]),'_Norm')
    #plot_top_view(allFluxes[:, 0], allFluxes[:, 1], allFluxes[:, 2])
    #plot_top_view(allFluxes[:, 0], allFluxes[:, 1], normalize(allFluxes[:,2]))
    if mrad == True:
            plot2D(allFluxes[:,0]/distanceFromSource, allFluxes[:,1]/distanceFromSource, allFluxes[:,2],'','mrad')
            plot2D(allFluxes[:,0]/distanceFromSource, allFluxes[:,1]/distanceFromSource, normalize(allFluxes[:,2]),'_Norm','mrad')
            plot2D_Log(allFluxes[:,0]/distanceFromSource, allFluxes[:,1]/distanceFromSource, allFluxes[:,2],'_Log','mrad')
            #saveToCSV(allFluxes[:,0], allFluxes[:,1], allFluxes[:,2],'mrad')
    else: 
            plot2D(allFluxes[:,0], allFluxes[:,1], allFluxes[:,2],'','mm')
            plot2D(allFluxes[:,0], allFluxes[:,1], normalize(allFluxes[:,2]),'_Norm','mm')
            plot2D_Log(allFluxes[:,0], allFluxes[:,1], allFluxes[:,2],'_Log','mm')
            #saveToCSV(allFluxes[:,0], allFluxes[:,1], allFluxes[:,2],'mm')

'''
x = allFluxes[:,0]
y = allFluxes[:,1]
z = allFluxes[:,2]
'''
