import os
import numpy as np
import glob
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from lmfit import Model, Parameters

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

th = 0

##generalised 2-factor Gaussian function
def gauss2(x, a1, b1, c1, a2, b2, c2):
    return a1*np.exp(-((x-b1)/c1)**2) + a2*np.exp(-((x-b2)/c2)**2)

plt.figure(figsize=(20,10))

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
    
    plt.plot(x1, y1, 'bo')
    
    peak = max(y1)
    base = min(y1)
    mean = sum(x1*y1)/sum(y1)
    sigma = np.sqrt((sum(y1*(x1 - mean)**2))/sum(y1))
    

    ##scipy.curve_fit method
    ##fit a curve to Gaussian function, 'popt' returns optimal values for a, b and c
    popt, pcov = curve_fit(gauss2, x1, y1, bounds=([base, (mean-sigma), 0, base, (mean-sigma), 0], [peak, (mean+sigma), 1, peak, (mean+sigma), 1]))
    
    plt.plot(x1, gauss2(x1, *popt))
    plt.plot([popt[1], popt[1]], [15, 45])
    plt.plot([popt[4], popt[4]], [15, 45])
    print(popt)
    
    
    
    if popt[1] > popt[4]:
        th = popt[1]
    else:th = popt[4] 
    
    
    
#    ##lmfit.Model method
#    gmodel = Model(gauss2)
#    params = Parameters()
#    params.add_many(('a1', (peak+base)/2, True, base, peak, None, None),
#                    ('b1', mean, True, (mean-sigma), (mean+sigma), None, None),
#                    ('c1', sigma, True, 0, 1, None, None),
#                    ('a2', (peak+base)/2, True, base, peak, None, None),
#                    ('b2', mean, True, (mean-sigma), (mean+sigma), None, None),
#                    ('c2', sigma, True, 0, 1, None, None)
#                    )
#    
#    result = gmodel.fit(y1, params, x=x1)
#    print(result.fit_report())
#
#    plt.plot(x1, result.init_fit)
#    plt.plot(x1, result.best_fit)
    
    ##Bragg's law: 2*d*sin(theta) = n*lambda
    deg = (th)/2
    rad = ((2*np.pi)/360)*deg
    d = (l/(2*np.sin(rad)))*1e10
      
    np.D.insert(i, d)
    np.Z.insert(i, z)
    
    i = i+1
    z = z+0.02
    
    

plt.plot(np.Z, np.D)
plt.xlim(-3, 3)
plt.ylim(1.535, 1.54)
plt.xlabel("z (mm)")
plt.ylabel("d spacing (Angstrom)")
#print(np.D)
#plt.savefig('testgraph.png')
    
