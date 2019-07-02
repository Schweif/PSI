#imports for 3D Ploting of collected data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
fig = plt.figure()
ax = fig.gca(projection='3d')

#imports for Gauss fit
import pylab as plb
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp


DataRaw =  '/home/just/Documents/PSI/Flux_5005_MS-SLS1.dat'
dist_from_source = 10 #m
import numpy as np #did not work on first imort doing it twice instead
da= np.genfromtxt(DataRaw,skip_header=10, usecols=(0, 1, 2))

xdata= da[:,0]
ydata= da[:,1]
zdata= da[:,2]


ax = plt.axes(projection='3d')
ax.scatter3D(xdata, ydata, zdata, c=zdata, cmap='Greens')
plt.show()




'''
# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
'''


#cut throu x=0
x0zdata=[]
x0ydata=[]

for d in da:
    if d[0] == 0:
        x0ydata.append(d[1])
        x0zdata.append(d[2])

#cut throu y=0
y0zdata=[]
y0xdata=[]

for d in da:
    if d[1] == 0:
        y0xdata.append(d[0])
        y0zdata.append(d[2])



#get Numpy arrays from list
x0ydata=np.array(x0ydata)
x0zdata=np.array(x0zdata)
y0xdata=np.array(y0xdata)
y0zdata=np.array(y0zdata)

# for x=0
nY = len(x0zdata)
meanY = sum(x0ydata*x0zdata)/ sum(x0zdata)
sigmaY = np.sqrt(sum(x0zdata*(x0ydata-meanY)**2)/sum(x0zdata))

# for y=0
nX = len(y0zdata)
meanX = sum(y0xdata*y0zdata)/ sum(y0zdata)
sigmaX = np.sqrt(sum(y0zdata*(y0xdata-meanX)**2)/sum(y0zdata))

def Gauss(x,a,x0,sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))

def FitAndPlot2DGauss(axis):
    if axis == 'x':
        sigma = sigmaX
        mean = meanX
        n = nX
        x= y0xdata
        y= y0zdata
    if axis == 'y':
        sigma = sigmaY
        mean = meanY
        n = nY
        x= x0ydata
        y= x0zdata
    if axis == 0:
        n= len(y)
        mean = sum(x * y) / sum(y)
        sigma = np.sqrt(sum(y * (x - mean)**2) / sum(y))
    print "***********************"
    print "* Results for "+axis +" -axis *"
    print "***********************"
    print "sigma = " +str(sigma) +" mm"
    print "sigma = " +str(sigma/dist_from_source) +" mrad"
    print "3 sigma = " +str(3*sigma) +" mm"
    print "3 sigma = " +str(3*sigma/dist_from_source) +" mrad"
    print "mean = " +str(mean)
    print "************************************************"
 
    popt,pcov = curve_fit(Gauss, x, y, p0=[max(y), mean, sigma])
    #print "popt = " +str(popt)
    #print "pcov = " +str(pcov)
    plt.plot(x, y, 'b+:', label='data')
    plt.plot(x, Gauss(x, *popt), 'r-', label='fit')
    plt.legend()
    plt.title(axis+' Axis')
    plt.xlabel('length (mm)')
    plt.ylabel('PowerDens (W/mm**2)')
    plt.text(1,0,"sigma = " +str(round(sigma,3)) +" mm", fontsize=12, horizontalalignment='left', verticalalignment='bottom')
    plt.show()    

    
FitAndPlot2DGauss('x')
FitAndPlot2DGauss('y')




'''

plt.plot(x0ydata,x0zdata,'b+:',label='data')
plt.plot(x0ydata,Gauss(x0ydata,*popt),'ro:',label='fit')
plt.legend()
plt.title('X=0, Gaussfit')
plt.xlabel('Y in mm')
plt.ylabel('Flux Desity')
plt.show()
'''
'''
check:https://stackoverflow.com/questions/19206332/gaussian-fit-for-python
import numpy as np
import matplotlib.pyplot as plt

DataRaw =  
testFile = open(DataRaw,'r')
lines = testFile.readlines()
testFile.close

data = []
i =0

for line in lines[10:]:
    data.append([line[:13], line[14:26], line[32:41]])
data

dat= dat.astype(np.float)

fig = plt.figure()
ax = plt.axes(projection='3d')
da= np.genfromtxt(DataRaw,skip_header=10, usecols=(0, 1, 2))
da[:,0].min()
######################'''

