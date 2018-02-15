import os
import numpy as np
import glob
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

os.chdir('C://Users//mbgnwlr2//Documents//PhD_SynchrotronStuff//Winter18_beamtime_data//xrd//Al_SiC_1D_0_')
filenames = glob.glob('Al_SiC_1D_0__*.dat')

numfiles = len(filenames)

a = 1114
b = 1164

##Calculate wavelength from beam energy
E = 60000       #eV
h = 4.136e-15   #eV.s
c = 3e+08       #m/s
    
l = (h*c)/E     #m
#l = 1.2398419292e-11

np.D = []
np.Z = []

i = 0
z = -2.98

for f in filenames:
   
    ##pull data from the file into an array
    data = np.genfromtxt(f,
                         dtype = float,
                         skip_header=52,
                         delimiter='    '
                         )
    
    x = data[:,0]
    y = data[:,1]
    
    x1 = x[a:b]
    y1 = y[a:b]
    
    n = len(x1)
    
#    plt.plot(x1, y1)
    
    peak = max(y1)
    base = min(y1)
    mean = sum(x1*y1)/sum(y1)
    sigma = np.sqrt((sum(y1*(x1 - mean)**2))/sum(y1))
    
    ##generalised Gaussian function
    def gauss(x, a1, b1, c1, a2, b2, c2):
        return a1*np.exp(-((x-b1)/c1)**2) + a2*np.exp(-((x-b2)/c2)**2)
    
    ##fit a curve to Gaussian function, 'popt' returns optimal values for a, b and c
    popt, pcov = curve_fit(gauss, x1, y1, bounds=([base, (mean-sigma), 0, base, (mean-sigma), 0], [peak, (mean+sigma), 1, peak, (mean+sigma), 1]))
    
#    plt.plot(x1, gauss(x1, *popt))
#    plt.plot([popt[1], popt[1]], [15, 45])
    
     

    
    ##Bragg's law: 2*d*sin(theta) = n*lambda
    rad = (popt[1])/2
    deg = (360/(2*np.pi))*rad
    d = l/(2*np.sin(deg))
      
    np.D.insert(i, d)
    np.Z.insert(i, z)
    
    i = i+1
    z = z+0.02

plt.plot(np.Z, np.D)
print(np.D)
#plt.savefig('testgraph.png')
    
