# Here I will calculate the contribution of a delay-bin 
# or a range of delay-bins to the total number of objects 
# in the survey.

import numpy as np
import dtdutils
reload(dtdutils)
import os
import sys
import string
DTDpath = os.getenv('DTD')

#This is for RR Lyrae
binningScheme = 'CrazyCeph'
#ages_edges = np.array([1300., 2000., 3200., 5000., 7900., 12600, 20000.])*1.0e6
age_cent, ages_edges = dtdutils.sfh_ageBins(binningScheme, return_age_edge_array=True)

dtd = np.array([4.4091889998768743e-05, 3.5933597204133906e-05, 4.4982926113994287e-05, 4.6493106599294147e-05, 1.3046380005460337e-05, 9.8019434341299551e-06, 1.6250090650115581e-07, 5.9715772026879523e-07])

nCells = 808
nAgeBins = age_cent.size
N = 3162.
#Regenerating the SFH Map
sfhFileName = DTDpath + '/MC_SFH_Maps/lmc_sfh.dat'
outPathName = DTDpath + '/Output_SFH_Files/'
objClassName = 'Cepheids'
obj_subtype = 'All'


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

for i in range(age_cent.size):
    N_i = (100./N)*dtd[i]*np.sum(sfhMap[i, :])
    print 'Contribution % for Ages ', ages_edges[i]/1.0e9, '-', ages_edges[i+1]/1.0e9, 'Gyrs = ', N_i
