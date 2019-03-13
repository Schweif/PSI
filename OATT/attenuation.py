#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

#March 2019
#Paul Scherrer Institut, PSI
#David Marco Just
#david.just@psi.ch


#calculates and displays all posible attenuations of a set of different filter foils to be used in the solid state attenuator at the ATHOS beam line at SwissFEL
#Maximum size of 24 foils can be corectly calculated



import numpy as np
from os import listdir, chdir
from os.path import isfile, join
import matplotlib.pyplot as plt
from itertools import combinations
import sys

pathToAttenuatorFoils = '/home/just/Documents/PSI/OATT/materials/selection18Kirsten/'
duplicateFoils = True  # specify if the foils will be duplicated for redundancy and for more possible combinations

# append the rawData from all Foils
i = 0
da = []
for f in listdir(pathToAttenuatorFoils):
    print f
    chdir(pathToAttenuatorFoils)
    if isfile(join(pathToAttenuatorFoils, f)) and f.endswith('.dat'):
        da.append(f[:-4])
        da.append(np.genfromtxt(f, skip_header=2, usecols=(0, 1)))
        i = i + 1
noOfFoils = i

if duplicateFoils == True:
    for f in listdir(pathToAttenuatorFoils):
        print f
        chdir(pathToAttenuatorFoils)
        if isfile(join(pathToAttenuatorFoils, f)) and f.endswith('.dat'):
            da.append(f[:-4])
            da.append(np.genfromtxt(f, skip_header=2, usecols=(0, 1)))
            i = i + 1
    noOfFoils = noOfFoils * 2

if noOfFoils > 23:
	sys.exit("Error the list of foils is larger than 24. The resulting combinations will be to many. Please lower the ammount of foils or addapt this script")	

# count from 0 to nOfFoils
foils = []
i = 0
while i < noOfFoils:
    foils.append(i)
    i = i + 1

# get all possible combinations
noOfFoilsToCombine = 2
possibleCombinations = []
while noOfFoilsToCombine <= noOfFoils or noOfFoilsToCombine <= 12:
    comb = combinations(foils, noOfFoilsToCombine)
    possibleCombinations.extend(list(comb))
    noOfFoilsToCombine = noOfFoilsToCombine + 1

# exclude invalid combinations:
newComb = []
for i in possibleCombinations:
    # exclude combinations on the same arm
    if set([0, 12]).issubset(i):
        continue
    elif set([1, 13]).issubset(i):
        continue
    elif set([2, 14]).issubset(i):
        continue
    elif set([3, 15]).issubset(i):
        continue
    elif set([4, 16]).issubset(i):
        continue
    elif set([5, 17]).issubset(i):
        continue
    elif set([6, 18]).issubset(i):
        continue
    elif set([7, 19]).issubset(i):
        continue
    elif set([8, 20]).issubset(i):
        continue
    elif set([9, 21]).issubset(i):
        continue
    elif set([10, 22]).issubset(i):
        continue
    elif set([11, 23]).issubset(i):
        continue
    else:
        newComb.append(i)

# display difference before and after excluding some combinations; debug only
noCombinations = len(possibleCombinations)
newNoCombinations = len(newComb)

# Update list of combinations after deleting foribiden combinations
possibleCombinations= newComb
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
plt.ylim((10e-5, 1))  # set the ylim to bottom, top
plt.xlim((250, 2000))  # set the xlim to left, right
plt.xlabel('Photon Energy [eV]')
plt.ylabel('Transmission')
plt.legend()
plt.show()
