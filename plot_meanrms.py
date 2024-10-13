"""
Created on Thu Jan 18 12:51:46 2024
Finished editing on 2/21/24 at 7:07 am
@author: aidan
"""
import numpy as np
import matplotlib.pyplot as plt
import argparse

# initialize parser (copied in part from Thea's code sub_nar_updated.py)
parser = argparse.ArgumentParser(description='puts continuum fluxes into list that BRAINS can read') 
parser.add_argument('objpath', metavar='obj_path', help='path to the folder containing all original combined spectra (probabaly "original" directory)')
parser.add_argument('season', metavar='season', help='season number, i.e. "season4"')
arg = parser.parse_args()

# AJF read in parameter and wave/flux data
param_path = arg.objpath+'parameters'
parameters = open(param_path).readlines()
z = [float(i.split()[1]) for i in parameters if i.split()[0] == 'z'][0] 
data = np.loadtxt(arg.objpath+'combined.txt.meanrms_'+arg.season, usecols = (0,1))
wave, flux = data[:,0], data[:,1]
wave /= (1+z)

# AJF quick plot the whole spectra, adjust x and y lims as needed
plt.plot(wave, flux,color = 'k', linestyle = '-')
#plt.xlim(5075, 5575)
#plt.ylim(1e-14, 2e-14)
plt.show()

