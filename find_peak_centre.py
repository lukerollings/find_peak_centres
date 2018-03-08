import os
import numpy as np
import glob
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import lfilter

os.chdir('C://Users//mbgnwlr2//Documents//PhD_SynchrotronStuff//Winter18_beamtime_data//xrd')

sample = "13D"

loads = glob.glob('Al_SiC_'+str(sample)+'_*')



a = 1114
b = 1164
    
##Calculate wavelength from beam energy
E = 60000       #eV
h = 4.136e-15   #eV.s
c = 3e+08       #m/s
        
l = (h*c)/E     #m
#l = 1.2398419292e-11

##generalised 2-factor Gaussian function
def gauss2(x, a1, b1, c1, a2, b2, c2):
    return a1*np.exp(-((x-b1)/c1)**2) + a2*np.exp(-((x-b2)/c2)**2)

I = 0

for L in loads:
    file = str(L)+'//'+str(L)+'_*.dat'
    filenames = glob.glob(file) 
    
    np.D = []
    np.Z = []
    
    i = 0
    z = -2.98
    
    th = 0
    
    for f in filenames:
       
        ##pull data from the file into an array
        data = np.genfromtxt(f,
                             dtype = float,
                             skip_header=54,
                             delimiter='    '
                             )
        
        x = data[:,0]
        y = data[:,1]
        
        x1 = x[a:b]
        y1 = y[a:b]
        

        
        peak = max(y1)
        base = min(y1)
        mean = sum(x1*y1)/sum(y1)
        sigma = np.sqrt((sum(y1*(x1 - mean)**2))/sum(y1))
        
    
        ##scipy.curve_fit method
        ##fit a curve to Gaussian function, 'popt' returns optimal values for a, b and c
        popt, pcov = curve_fit(gauss2, x1, y1, bounds=([base, (mean-sigma), 0, base, (mean-sigma), 0], [peak, (mean+sigma), 1, peak, (mean+sigma), 1]))
        
#        a1, b1, c1, a2, b2, c2 = popt[0], popt[1], popt[2], popt[3], popt[4], popt[5]
#        
#        popt, pcov = curve_fit(gauss2, x1, y1, bounds=([base, (mean-sigma), 0, base, (mean-sigma), 0], [peak, (mean+sigma), 1, peak, (mean+sigma), 1]))
        
#        ##For single files, remove comments if you wish to display the fit
#        plt.plot(x1, y1, 'bo')
#        plt.plot(x1, gauss2(x1, *popt))
#        plt.plot([popt[1], popt[1]], [15, 45])
#        plt.plot([popt[4], popt[4]], [15, 45])
#        print(popt)
        
        ##Returns the highest of the caluclated means as the 2theta value 'th'
        if popt[1] > popt[4]:
            th = popt[1]
        else:th = popt[4] 
        
        
        ##Bragg's law: 2*d*sin(theta) = n*lambda
        deg = (th)/2
        rad = ((2*np.pi)/360)*deg
        d = (l/(2*np.sin(rad)))
        e = (d-1.54290106970688e-10)/1.54290106970688e-10       #convert from d spacing to strain using value calculated from exposed fibre
          
        np.D.insert(i, e)
        np.Z.insert(i, z)
        
        i = i+1
        z = z+0.02
        
        
    plt.figure(1, figsize=(20,10))
    
    plt.plot(np.Z, np.D, label=str(sample)+'_load_stage_'+str(I))
    plt.xlim(-3, 3)
    plt.ylim(-0.001, 0.01)
    plt.xlabel("z (mm)")
    plt.ylabel("strain, $\epsilon$")
    plt.title(str(sample))
    plt.legend(loc='upper left')
    plt.savefig(str(sample)+'_strain.png')
    
    plt.figure(2, figsize=(20,10))
    grad = np.gradient(np.D)
    
    N = 15
    B = [1.0/N]*N
    A = 1
    
    dedz = lfilter(B, A, grad) #array with values of de/dz, smoothed out with lfilter
    
    plt.plot(np.Z, dedz, label=str(sample)+'_load_stage_'+str(I))
    plt.xlim(-3, 3)
    plt.ylim(-0.0001, 0.00015)
    plt.xlabel("z (mm)")
    plt.ylabel("d$\epsilon$/dz")
    plt.title(str(sample))
    plt.legend(loc='upper left')
    plt.savefig(str(sample)+'_gradient.png')
    
    I = I+1
    


    
