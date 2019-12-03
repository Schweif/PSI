#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

#March 2019
#Paul Scherrer Institut, PSI
#David Marco Just
#david.just@psi.ch


#calculates and displays all posible attenuations of a set of different filter foils to be used in the solid state attenuator at the ATHOS beam line at SwissFEL
#Maximum size of 24 foils can be corectly calculated
#TODO: Funtion which depends on number of arms is not properly implemented



import numpy as np
from os import listdir, chdir, path
from os.path import isfile, join
import matplotlib.pyplot as plt
from itertools import combinations
import sys
import re

pathToAttenuatorFoils = '/home/just/Documents/PSI/OATT/materials/selection18Kirsten_2d/'
duplicateFoils = False  # specify if the foils will be duplicated for redundancy and for more possible combinations
nArms = 6   #Specifiy the ammount of arms to be equiped with foils 
autoSave = True
fname = 'Attenuation_d2.png'



def sorted_aphanumeric(data):
    """Sorts the Foils according to theitr alpha numeric order. I.e. allows to specify which foil to use where e.g. 01_Si2u.dat;02_Al_1u.dat"""
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(data, key=alphanum_key)


# append the rawData from all Foils
i = 0
da = []

files = listdir(pathToAttenuatorFoils)
files = sorted_aphanumeric(files)
for f in files:
    print f
    chdir(pathToAttenuatorFoils)
    if isfile(join(pathToAttenuatorFoils, f)) and f.endswith('.dat'):
        da.append(f[:-4])
        da.append(np.genfromtxt(f, skip_header=2, usecols=(0, 1)))
        i = i + 1
noOfFoils = i

if duplicateFoils == True:
    for f in listdir(pathToAttenuatorFoils):
        print f[:-4]
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

# exclude invalid combinations, for 6 independend arms:
newComb = []
for i in possibleCombinations:
    # exclude combinations on the same arm
    if set([0, 1]).issubset(i):
        continue
    elif set([0, 2]).issubset(i):
        continue
    elif set([1, 2]).issubset(i):
        continue
    elif set([3, 4]).issubset(i):
        continue
    elif set([3, 5]).issubset(i):
        continue
    elif set([4, 5]).issubset(i):
        continue
    elif set([6, 7]).issubset(i):
        continue
    elif set([6, 8]).issubset(i):
        continue
    elif set([7, 8]).issubset(i):
        continue
    elif set([9, 10]).issubset(i):
        continue
    elif set([9, 11]).issubset(i):
        continue
    elif set([10, 11]).issubset(i):
        continue
    elif set([12, 13]).issubset(i):
        continue
    elif set([12, 14]).issubset(i):
        continue
    elif set([13, 14]).issubset(i):
        continue
    elif set([15, 16]).issubset(i):
        continue
    elif set([15, 17]).issubset(i):
        continue
    elif set([16, 17]).issubset(i):
        continue
    else:
        newComb.append(i)

# exclude invalid combinations, for 12 independend arms:
'''
newComb = []
for i in possibleCombinations:
    # exclude combinations on the same arm
    if set([0, nArms]).issubset(i):
        continue
    elif set([1, nArms+1]).issubset(i):
        continue
    elif set([2, nArms+2]).issubset(i):
        continue
    elif set([3, nArms+3]).issubset(i):
        continue
    elif set([4, nArms+4]).issubset(i):
        continue
    elif set([5, nArms+5]).issubset(i):
        continue
    elif set([6, nArms+6]).issubset(i):
        continue
    elif set([7, nArms+7]).issubset(i):
        continue
    elif set([8, nArms+8]).issubset(i):
        continue
    elif set([9, nArms+9]).issubset(i):
        continue
    elif set([10, nArms+10]).issubset(i):
        continue
    elif set([11, nArms+11]).issubset(i):
        continue
    else:
        newComb.append(i)
'''

# display difference before and after excluding some combinations; debug only
noCombinations = len(possibleCombinations)
newNoCombinations = len(newComb)

# Update list of combinations after deleting foribiden combinations
possibleCombinations= newComb
noCombinations = len(possibleCombinations)
print 'Possible combinations using this set of foils: ' +str(noCombinations)

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
if autoSave == True:
    plt.savefig(fname)
    plt.close()
    plt.clf()
else:
    plt.show()
    plt.clf()
