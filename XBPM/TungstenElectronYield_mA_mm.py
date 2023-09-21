#!/usr/bin/python3.7
# -*- coding: utf-8 -*-

# Nov 2021
# Paul Scherrer Institut, PSI
# David Marco Just
# david.just@psi.ch


import numpy as np
from os import listdir, chdir, path
from os.path import isfile, join
from matplotlib import cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import scipy
import re
import pandas as pd
import json
from shutil import copy

#CONFIGURATION
CRXO = False #Select source type for Yield data. If True, source is X-Ray Attenuation Length from CRXO website if false it is from Evaluated Nuclear Data File Libary
if CRXO == True:
    pathToYield = path.dirname(path.realpath(__file__))+r".\RawData\W.txt"
else:
    #pathToYield = path.dirname(path.realpath(__file__))+r".\RawData\EPDL97_74.dat"
    pathToYield = r"C:\Users\just_d\Documents\Python Scripts\XBPM\Raw Data\Yield\EPDL97_74.dat"
           
 
pathToFluxes =  r"H:\00_final files\IDFE\X11MFE\5_SLS2.0\02_Berechnungen\02_XBPM\FluxDensity\01_LH\Raw"
title = 'UE36kn Test verified calculation'
autoSave = True # Set True to automatically save the plots in pathToFluxes
yLabel = "Current density [mA/mm^2]"
maxEnergy = 30000.0  # eV given by flux calculations
minEnergy = 30.0  # eV given by yield table
distanceFromSource = False #  m distance from source at which the flux was calculated, if False distance will be extracted from SPECTRAs .json,
plotDistance = False  # m distance from source at which the result is wanted, if False distanceFromSource will be taken
bucketSize = 1 #  eV
mrad = False #Sets X and Y axis in Plots to mrad. Set to False if X and Y axis in plots should be mm at distanceFromSource
calibrationFactor = 1.8E-8 #Factor that calibrates the photo electric cross section * energy [cm2*eV/g] to yield [e-/ph]
ec = 1.602e-19 #elementary charge [C]
exampleFluxDensData = r"U:\10_Skripts\XBPM\RawData\TestDataSets\SmallDataSet\UE36kn_mm2"
maxZ=30 #Defines the maximum level of the Z axis. I.e. higher vallues than maxZ will be read in the plot, set to 0 if not used
#CONFIGURATION ENDS HERE


def prepare_yield_data_CRXO_cal(pathToYield):
    """Reads in the electron yield file, converts the attenuation length to cross section to mm²/g, removes unneeded vallues and returns a data set with (energy [eV]:cross section [mm²/g])"""
    yieldPerEnergy = []
    yieldPerEnergy = np.genfromtxt(pathToYield, skip_header=2, usecols=(0, 1))
    correctionFactor = 5*100
    #  set energy to eV
    yieldPerEnergy[:, 0] = yieldPerEnergy[:, 0]
    #  set yield mA / ph
    yieldPerEnergy[:, 1] = (1 / yieldPerEnergy[:, 1] * correctionFactor)* yieldPerEnergy[:, 0] * calibrationFactor *ec *1000
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

def prepare_yield_data_cal(pathToYield):
    """Reads in the electron yield file, converts energies to eV instead of MeV and the cross section to mm²/g instead of cm²/g, removes unneeded vallues and returs a data set with (energy [eV]:cross section [mm²/g]"""
    yieldPerEnergy = []
    yieldPerEnergy = np.genfromtxt(pathToYield, skip_header=18, usecols=(0, 9))

    #  set energy to eV
    yieldPerEnergy[:, 0] = yieldPerEnergy[:, 0] * 1000000.0
    #  set yield mA / ph
    yieldPerEnergy[:, 1] = yieldPerEnergy[:, 1] * yieldPerEnergy[:, 0] * calibrationFactor * ec * 1000
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
    if not path.exists(pathToFluxes):
            pathToFluxes = exampleFluxDensData
            print('ATTTENTION: Raw flux density data not found. Using example data set instead')
    chdir(pathToFluxes)
    da = []
    i = 0
    y = 0
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
            """prepares Spectra V11 and greater data sets"""
            with open(f) as v:
                data = json.load(v)
                if 'Output' in data:
                    if y == 0:
                        dist_from_source, mr2_or_mm2 = getSettingsFromSpectraFiles(data)
                    y = y + 1
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
    return(da, noEnergies, dist_from_source, mr2_or_mm2)


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
    weightedFluxes[:,2]=weightedFluxes[:,2]/ (plotDistance)**2
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
    ax.set_zlabel("Current density [A/mm^2] ")
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)    
    #plt.xlim(-0.5,0.5)
    #plt.ylim(-0.5,0.5)
    #plt.axis('scaled')
    if autoSave == True:
        plt.savefig(fname, dpi = 1200)
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
        plt.savefig(fname, dpi = 1200)
        plt.close()
        plt.clf()
    else:
        plt.show()
        plt.clf()


def plot2D_BackUp(x, y, z, txt='', unit=''):
    fname = title + '_2D' + txt + '.png' 
    xmax= np.max(x)
    ymax= np.max(y)
    xmin= np.min(x)
    ymin= np.min(y)

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
        plt.savefig(fname, dpi = 1200)
        plt.close()
        plt.clf()
    else:
        plt.show()
        plt.clf()


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
    if zMax != 0:
        im = ax.imshow(Z,
                       cmap=cm.rainbow, 
                       interpolation=  'bilinear', #'none',
                       origin='lower', extent=[xmin, xmax, ymin, ymax],
                       vmax=zMax, vmin=-z.min())
    else:
        im = ax.imshow(Z,
                       cmap=cm.rainbow, 
                       interpolation=  'bilinear', #'none',
                       origin='lower', extent=[xmin, xmax, ymin, ymax],
                       vmax=z.max(), vmin=-z.min())                   
    
    plt.title(fname,pad=25)
    plt.xlabel('x, position hor. ' +unit)
    plt.ylabel('y, position ver. ' +unit)
    cbar = plt.colorbar(im, ax=ax)
    #ToDO yLabel in Normalized plot is wrong
    #if txt == '_Norm':
    #    yLabel = 'Power denstiy AU'
    cbar.ax.set_ylabel(yLabel)
    if autoSave == True:
        plt.savefig(fname, dpi = 1200)
        plt.close()
        plt.clf()
    else:
        plt.show()
        plt.clf()
 
    
def plot2D_Log(x, y, z, txt='', unit=''):
    fname = title + txt + '.png' 
    xmax= np.max(x)
    ymax= np.max(y)
    xmin= np.min(x)
    ymin= np.min(y)
    X, Y = np.meshgrid(x, y)

    N = int(len(z)**.5)
    Z = z.reshape(N, N)

    fig, ax = plt.subplots()
    im = ax.imshow(Z, 
                   norm=colors.LogNorm(vmin=Z.min(), vmax=Z.max()),
                   cmap=cm.rainbow, 
                   interpolation=  'bilinear', #'none',
                   origin='lower', extent=[xmin, xmax, ymin, ymax]
                   )
    
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
        
        
def normalize(data):
    norm = (data-data.min())/(data.max()-data.min())
    return norm

def saveToCSV(x, y, z, unit='mm'):
    d= {'x horizontal '+unit: x, 'y horizontal '+unit: y, yLabel : z}
    df = pd.DataFrame(data=d)
    df.to_csv(title + '.csv')
    
def saveToXLSX(x, y, z, unit='mm'):
    d= {'x horizontal '+unit: x, 'y horizontal '+unit: y, yLabel : z}
    df = pd.DataFrame(data=d)
    df.to_excel(title + '.xlsx')

def getSettingsFromSpectraFiles(data): 
    """Extracts the distance to source at which the Spectra calculation was performed and 
    finds out wheteher the flux density was stored as ph/s/mr^2/0.1%B.W. or ph/s/mm^2/0.1%B.W. in the Spectra data files (.json)
    input being a read in .json (data = json.load(f); f = open('FluxDensity-22_0.json')"""
                                 
    dist_from_source = data['Input']['Configurations']['Distance from the Source (m)']
    mr2_or_mm2 = data['Output']['units'][2]
    if 'mm^2' in mr2_or_mm2:
        mr2_or_mm2 = 'mm2'
    elif 'mr^2' in mr2_or_mm2:
        mr2_or_mm2 = 'mr2'
    else:
       print('Warning could not determine if flux is ph/s/mr^2/0.1%B.W. or ph/s/mm^2/0.1%B.W.') 
    return(dist_from_source,mr2_or_mm2)
       
    
        

if __name__ == '__main__':
    script = path.abspath(__file__) 
    if CRXO == True:
        yieldPerEnergy = prepare_yield_data_CRXO_cal(pathToYield)
    else:
        yieldPerEnergy = prepare_yield_data_cal(pathToYield)
    fluxData, noEnergies, dist_from_source, mr2_or_mm2 = read_flux_data(pathToFluxes, minEnergy, maxEnergy)
    if not distanceFromSource:
        distanceFromSource = dist_from_source
    if not plotDistance:
        plotDistance = distanceFromSource 
    print( str(noEnergies) +' energy data sets read in.')
    fluxData = photons_per_energy_bucket(fluxData, bucketSize)
    fluxDataYielded = multiply_flux_with_yield(fluxData,yieldPerEnergy)
    allFluxes= integrate_all_weigthed_fluxes(fluxDataYielded)
    if mr2_or_mm2 == 'mr2':
        allFluxes = flux_per_mm_sqr(allFluxes, distanceFromSource)
    if mrad == True: #Z axis is allways per mm since flux_per_mm_sqr() is used
            plot2D(allFluxes[:,0]/distanceFromSource, allFluxes[:,1]/distanceFromSource, allFluxes[:,2],'','mrad')
            plot2D(allFluxes[:,0]/distanceFromSource, allFluxes[:,1]/distanceFromSource, normalize(allFluxes[:,2]),'_Norm','mrad')
            plot2D_Log(allFluxes[:,0]/distanceFromSource, allFluxes[:,1]/distanceFromSource, allFluxes[:,2],'_Log','mrad')
    else: #Z axis is allways per mm since flux_per_mm_sqr() is used
            plot2D(allFluxes[:,0]/distanceFromSource*plotDistance, allFluxes[:,1]/distanceFromSource*plotDistance, allFluxes[:,2],'','mm')
            plot2D(allFluxes[:,0]/distanceFromSource*plotDistance, allFluxes[:,1]/distanceFromSource*plotDistance, normalize(allFluxes[:,2]),'_Norm','mm')
            plot2D_Log(allFluxes[:,0]/distanceFromSource*plotDistance, allFluxes[:,1]/distanceFromSource*plotDistance, allFluxes[:,2],'_Log','mm')
            saveToCSV(allFluxes[:,0], allFluxes[:,1], allFluxes[:,2],'mm')
            saveToXLSX(allFluxes[:,0], allFluxes[:,1], allFluxes[:,2],'mm')
            if maxZ != 0:
                plot2D(allFluxes[:,0]/distanceFromSource*plotDistance, allFluxes[:,1]/distanceFromSource*plotDistance, allFluxes[:,2],' crpt @' +str(maxZ),'mm', maxZ)
    copy(script, title +'.py') #save a coppy of the script allong with the plots for quality control

