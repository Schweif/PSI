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
from itertools import combinations

pathToYield = '/home/just/Documents/PSI/XBPM/rawData/EPDL97_74.dat'
pathToFluxes = '/home/just/Documents/PSI/XBPM/rawData/X04S_flux/'
maxEnergy = 30000.0 # eV given by flux calculations
minEnergy = 100.0 # eV given by yield table


def prepare_Yield_Data(pathToYield):
    yieldPerEnergy = []
    yieldPerEnergy = np.genfromtxt(pathToYield, skip_header=18, usecols=(0, 8))

    #  set Energy to eV
    yieldPerEnergy[: , 0] = yieldPerEnergy[ : , 0] * 1000000.0
    # remove energies lower than min and higher tha max
    i=0
    for E in yieldPerEnergy[:,0]:
        if E < minEnergy or E > maxEnergy:
            yieldPerEnergy = np.delete(yieldPerEnergy, i, 0)
            i=i-1
        i=i+1
    return(yieldPerEnergy)

def preapre_Flux_Data():
    chdir(pathToAttenuatorFoils)
    i= 0
    for f in listdir(pathToFluxes):
        data = []
        da = []
        if isfile(join(pathToFluxes, f)) and f.endswith('.dta'):
            fo= open(f, "r")
            lines = list(fo)
            Eline = lines[5]
            fo.close()
            Energy=Eline[-11:-2]
            if i = 0:
                da = np.genfromtxt(f, skip_header= 10, usecols=(0,1,3))
                EArray= np.full(da.shape[0],1), Energy)



if __name__ == '__main__':
    yieldPerEnergy = prepare_Yield_Data(pathToYield)
    print yieldPerEnergy

