# Here I will calculate the contribution of a delay-bin 
# or a range of delay-bins to the total number of objects 
# in the survey.

import numpy as np
import dtdutils
import os
import string
DTDpath = os.getenv('DTD')

#This is for RR Lyrae

ages_edges = np.array([1300., 2000., 3200., 5000., 7900., 12600, 20000.])*1.0e6
dtd = np.array([1.19, 1.95, 1.71, 2.12, 0.94, 1.03])*1.0e-5
dtd_plus = np.array([0.2, 0.13, 0.24, 0.2, 0.16, 0.09])*1.0e-5
dtd_minus = np.array([0.25, 0.11, 0.23, 0.25, 0.21, 0.07])*1.0e-5

#For Planetary Nebula

#ages_edges = np.array([5.0, 35.0, 800., 2000., 5000., 8000., 20000.])*1.0e6
#dtd = np.array([1.4, 0.6])*1.0e-6
#dtd_plus = np.array([0.8, 0.1])*1.0e-6
#dtd_minus = np.array([1.0, 0.2])*1.0e-6

#Prefixes

nCells = 808
nAgeBins = 16
N = 23459.
#Regenerating the SFH Map
sfhFileName = DTDpath + '/MC_SFH_Maps/lmc_sfh.dat'
outPathName = DTDpath + '/Output_SFH_Files/'
objClassName = 'RRLyrae'
obj_subtype = 'All'
binningScheme = 'Unbinned'

sfhMap = np.zeros((nAgeBins,nCells))
sfhMapMin = np.zeros((nAgeBins,nCells))
sfhMapMax = np.zeros((nAgeBins,nCells))
sfhMapRange = np.zeros((nAgeBins,nCells))
with open(outPathName+'LMC_SFH_Cells_'+objClassName+obj_subtype+'_'+binningScheme+'.dat', 'r') as f:
    for cell in np.arange(nCells) :
        cellLineWords = string.split(f.readline().strip())
        sfhMap[:,cell] = np.asarray(cellLineWords[2:]).astype(np.float)
        cellLineWords = string.split(f.readline().strip())
        sfhMapMin[:,cell] = np.asarray(cellLineWords).astype(np.float)
        cellLineWords = string.split(f.readline().strip())
        sfhMapMax[:,cell] = np.asarray(cellLineWords).astype(np.float)

sfhMapPlus = sfhMapMax - sfhMap
sfhMapMinus = sfhMap - sfhMapMin
print 'Total LMC stellar mass (x 10^9 M_sun) = ', np.sum(sfhMap)/1.0e9

for i in [10, 11, 12, 13, 14, 15]:
    N_i = (100./N)*dtd[i-10]*np.sum(sfhMap[i, :])
    sigma_N_i_p = (100./N)*np.sqrt(((dtd_plus[i-10]**2)*(np.sum(sfhMap[i,:]**2))) + ((dtd[i-10]**2)*(np.sum(sfhMapPlus[i,:]**2))))
    sigma_N_i_n = (100./N)*np.sqrt(((dtd_minus[i-10]**2)*(np.sum(sfhMap[i,:]**2))) + ((dtd[i-10]**2)*(np.sum(sfhMapMinus[i,:]**2))))
    print 'Contribution % for Ages ', ages_edges[i-10]/1.0e9, '-', ages_edges[i-10+1]/1.0e9, 'Gyrs = ', N_i, sigma_N_i_p, sigma_N_i_n
