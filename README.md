# find_peak_centres
Series of scripts for analysing localised strains from x-ray diffraction data.

The scripts read multiple directories containing a series of .dat files produced from integration of diffraction rings. In each case, the code highlights the peak in each .dat file between pre-set channels, equivalent to the relevant diffraction peaks of interest. The peak is fitted to a one ore two factor gaussian (depending on background noise)to determine the peak centre, giving a value of the crystallographic d-spacing (calculated by Bragg's law). This is used to calculate the local tensile or shear strain at any given point, and plotted into strain "maps" for each sub-directory (load stage).
Returns .png images of the smoothed and raw maps, and the data in two .csv dataframes for further analysis
