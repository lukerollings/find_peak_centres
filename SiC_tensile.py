# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 20:21:46 2018

@author: mbgnwlr2
"""

import os
import numpy as np
import glob
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit


##DEFINE VARIABLES AND FUNCTIONS

#redirect to data folder
os.chdir("##INSERT ROOT DIRECTORY##")

#Enter sample ID (sample 13D gave most useful results)
sample = "13D"

#Pulls all folders labelled with the sample ID into a list of strings
loads = glob.glob("Al_SiC_"+str(sample)+"_*")

#Initiates empty arrays for storing strain data, including smoothed and raw values
strain = np.empty((300, len(loads)))
strain_raw = np.empty((300, len(loads)))

#Define the region of interest for finding a peak
#For SiC peaks, points 1114-1164 cover the (1 0 8) peak around 2theta = 7.68deg
a = 1114
b = 1164

##Calculate wavelength from beam energy
E = 60000       #eV
h = 4.136e-15   #eV.s
c = 3e+08       #m/s
        
l = (h*c)/E     #m

#Generalised 2-factor Gaussian function
#Defined for peak fitting algorithm
def gauss2(x, a1, b1, c1, a2, b2, c2):
    return a1*np.exp(-((x-b1)/c1)**2) + a2*np.exp(-((x-b2)/c2)**2)

#Iteration counter to track load count
I = 0


##LOAD CYCLE

#Read list "loads", for each item, open folder and produce an array to store individual .dat files
for L in loads:
    print(L)
    file = str(L)+"//"+str(L)+"_*.dat"
    
    #Produce a list of .dat files in each folder identified in "loads"
    #Stores these in a list names "filenames"
    filenames = glob.glob(file)
    
    #Ensures code won't crash if it encounters a folder with less than 300 files, just moves on to the next value in "loads"
    if(len(filenames)<300):
        continue
    
    #Initiate arrays for storing d-spacing values (smooth and raw) and points along the z axis
    D = []
    D_raw = []
    Z = []
    
    #Reset iteration counters for file counter and z position
    i = 0
    z = -2.98   #mm
    
    #Reset value of th
    th = 0
    
    
    ##READING CYCLE
    
    for f in filenames:
        
        #Pull values from .dat files specified in list "filenames"
        #Reads each value as floating point, starting at line 56
        data = np.genfromtxt(f,
                             dtype = float,
                             skip_header = 56,
                             delimiter= "    "
                             )
        
        #Separate read data into x and y arrays
        x = data[:,0]
        y = data[:,1]
        
        #Snip just the region of interest
        x1 = x[a:b]
        y1 = y[a:b]
        
        #read the max, min and approx mean and sigma (required for gaussian fitting)
        peak = max(y1)
        base = min(y1)
        mean = sum(x1*y1)/sum(y1)
        sigma = np.sqrt((sum(y1*(x1 - mean)**2))/sum(y1))
        
        
        ##CURVE FITTING
        #scipy.curve_fit is used to fit the data in the region of interest to the pre defined function, in this case a 2-factor gaussian
        #The fit uses the calculated max, min, mean and sigma as initial estimates
        popt, pcov = curve_fit(gauss2, x1, y1, bounds=([base, (mean-sigma), 0, base, (mean-sigma), 0], [peak, (mean+sigma), 1, peak, (mean+sigma), 1]))
        
        
        ##Returns the highest of the caluclated means as the 2theta value 'th'
        if popt[1] > popt[4]:
            th = popt[1]
        else:th = popt[4] 
        
        #conversion from 2theta to local strain using Bragg's Law
        deg = (th)/2
        rad = ((2*np.pi)/360)*deg
        d = (l/(2*np.sin(rad)))
        e = (d-1.54290106970688e-10)/1.54290106970688e-10       #convert from d spacing to strain using value calculated from exposed fibre
        
        #Populate D and D_raw arrays with calculated values of local strain, and Z with relevant z axis position
        D.insert(i, e)
        D_raw.insert(i, e)
        Z.insert(i, z)
        
        #Replaces outlying points with previous points, gives smoother curve with no noise
        #May need changing to a more scientifically accurate/reliable method
        if abs((D[i]) - (D[i-1])) > 0.0005:
            D[i]=D[i-1]
        
        #Increase count on data ID and z position counters
        i = i+1
        z = z+0.02
        
    #Define new arrays for storing strain values, store to a dataframe
    strain[:,I] = D
    strain_raw[:,I] = D_raw
    
    df_strain = pd.DataFrame(strain, columns=loads)
    df_strain_raw = pd.DataFrame(strain_raw, columns=loads)
    
    
    ##DEFINE PLOT PARAMETERS
    #Set names of each I iteration - required for an accurate legend
    if I == 0:
        name = "6MPa"
        mark = "."
    elif I == 1:
        name = "105MPa"
        mark = "s"
    elif I == 2:
        name = "105MPa"
        mark = "^"
    elif I == 3:
        name = "111MPa"
        mark = "+"
    elif I == 4:
        name = "124MPa"
        mark = "x"
    elif I == 5:
        name = "6MPa"
        mark = "D"
    
    #Plot 1 - smoothed tensile strain    
    plt.figure(1, figsize = (20,10))
    
    plt.scatter(Z, D, label = name, marker = mark)
    plt.xlim(-3, 3)
    plt.ylim(-0.005, 0.01)
    plt.xlabel("z (mm)")
    plt.ylabel("strain, $\epsilon$")
    plt.grid(True)
    plt.legend(loc='upper left')
    plt.savefig(str(sample)+'_strain.png')
    
    #Plot 2 - raw tensile strain
    plt.figure(2, figsize = (20,10))
    
    plt.scatter(Z, D_raw, label = name, marker = mark)
    plt.xlim(-3, 3)
    plt.ylim(-0.005, 0.01)
    plt.xlabel("z (mm)")
    plt.ylabel("strain, $\epsilon$")
    plt.grid(True)
    plt.legend(loc='upper left')
    plt.savefig(str(sample)+'_strain_raw.png')
    
    #Add to load count iteration counter
    I = I+1

#After all iterations, save dataframes to .csv files for further analysis
os.chdir("C://Users//mbgnwlr2//Documents//Autumn18_FibreFrag//fig4")
df_strain.to_csv(str(sample)+'_strain.csv', ',')
df_strain_raw.to_csv(str(sample)+'_strain_raw.csv', ',')
