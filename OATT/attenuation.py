#!/usr/bin/python2.7
# -*- coding: utf-8 -*-


import numpy as np
from os import listdir, chdir
from os.path import isfile, join
import matplotlib.pyplot as plt
from itertools import combinations

pathToAttenuatorFoils = '/home/just/Documents/PSI/OATT/materials/CandSiOnly/selection20'


# append the rawData from all Foils
i = 0
da = []
for f in listdir(pathToAttenuatorFoils):
    chdir(pathToAttenuatorFoils)
    if isfile(join(pathToAttenuatorFoils, f)) and f.endswith('.dat'):
        da.append(f[:-4])
        da.append(np.genfromtxt(f, skip_header=2, usecols=(0, 1)))
        i = i + 1
noOfFoils = i

# count from 0 to nOfFoils
foils = []
i = 0
while i < noOfFoils:
    foils.append(i)
    i = i + 1

# get all possible combinations
noOfFoilsToCombine = 2
possibleCombinations = []
while noOfFoilsToCombine <= noOfFoils:
    comb = combinations(foils, noOfFoilsToCombine)
    possibleCombinations.extend(list(comb))
    noOfFoilsToCombine = noOfFoilsToCombine + 1
noCombinations = len(possibleCombinations)

for combination in possibleCombinations:
    nFoils = len(combination)
    i = 0
    nameCombination = ''
    while i < nFoils:
        # create new name for data Set
        foilNumber = combination[i]
        nameFoil = da[foilNumber * 2]
        if nameCombination == '':
            nameCombination = nameFoil
        else:
            nameCombination = nameCombination + '+' + nameFoil
        i = i + 1
    da.append(nameCombination)

    xDataFoil = []
    yDataFoil = []
    yDataCombination = []
    i = 0
    while i < nFoils:
        # create new data Set
        foilNumber = combination[i]
        xDataFoil = da[foilNumber * 2 + 1][:, 0]
        yDataFoil = da[foilNumber * 2 + 1][:, 1]
        if yDataCombination == []:
            yDataCombination = yDataFoil
        else:
            yDataCombination = yDataCombination * yDataFoil
        i = i + 1
    newData = np.column_stack([xDataFoil, yDataCombination])
    da.append(newData)

# create plot
i = 0
while i + 1 <= noOfFoils + noCombinations:
    x = da[i * 2 + 1][:, 0]
    y = da[i * 2 + 1][:, 1]
    vLabel = da[i * 2]
    plt.scatter(x, y)
    #plt.scatter(x, y, label=vLabel)
    i = i + 1

plt.semilogy(10)
plt.ylim((10e-6, 1))  # set the ylim to bottom, top
plt.xlim((250, 2000))  # set the xlim to left, right
plt.xlabel('Photon Energy [eV]')
plt.ylabel('Transmission')
plt.legend()
plt.show()
