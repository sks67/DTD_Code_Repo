#==================================================================#
#==================================================================#
#
# Main SFH map reading script. Will look at how to treat SFH errors
# from here
# 
# Date Created: 07/03/2017
# Copyright: Sumit Sarbadhicary
#
#==================================================================#
#==================================================================#




import corner
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii
#import dtdutils
#reload(dtdutils)
#import dtdplotutils as plot_dtd
#reload(plot_dtd)
from scipy import interpolate
from scipy import stats
from scipy import misc
from scipy import integrate
import string
import emcee

sfhFileName = '/Users/sumits2k/Desktop/Research/SNResearch2/RadioSNRs/DTD/MC_SFH_Maps/lmc_sfh.dat'
outPathName = '/Users/sumits2k/Desktop/Research/SNResearch2/RadioSNRs/DTD/Output_SFH_Files/'
print 'Reading SFH from file '+sfhFileName
nCells = 0
nAgeBins = 16
cellNames = []
cellCentersRA = []
cellCentersDec = []
sfhMap = []
with open(sfhFileName, 'r') as f:
    for i in range(17): f.readline()
    while True:
        line = f.readline() #Hyphens
        if not line : 
            print 'End of file reached after cell ', nCells
            break #Check EOF
        
        words = f.readline().split()
        cellNames.append(words[1])  #Cell name    
        #print 'Reading cell ' + cellNames[-1]
        words = f.readline().split() #Cell center
        cellCentersRA.append(float(words[1][0:2]) + float(words[2][0:2])/60.0)
        cellCentersDec.append(float(words[3][0:3]) - float(words[4][0:2])/60.0)
        line = f.readline() #Hyphens
        sfhCell = np.zeros([3,nAgeBins])
        for bin in range(nAgeBins) :
            floats = [float(x) for x in f.readline().split()]
            floatArr = np.asarray(floats)
            sfhCell[0,nAgeBins-1-bin] = floatArr[[1,4,7,10]].sum() #Best Fit
            sfhCell[1,nAgeBins-1-bin] = floatArr[[2,5,8,11]].sum() #Lower limit
            sfhCell[2,nAgeBins-1-bin] = floatArr[[3,6,9,12]].sum() #Upper limit
            f.readline() #Blank line

        sfhMap.append(sfhCell)
        nCells += 1
        
sfhMapArr = np.asarray(sfhMap)


#Selecting SFH ages and bins.
logAgeArr = np.asarray([6.8, 7.1, 7.4, 7.7, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0, 10.2])
ageArr = 10.0**logAgeArr
logAgeLimsArr = np.asarray([6.65,6.95,7.25,7.55,7.85,8.1,8.3,8.5,8.7,8.9,9.1,9.3,9.5,9.7,9.9,10.1,10.3])
ageLimsArr = 10.0**logAgeLimsArr
ageIntervalsArr = ageLimsArr[1:]-ageLimsArr[:-1]

np.set_printoptions(suppress=False)
tot_st_mass = np.sum(np.dot(sfhMapArr[:,0,:]/1.0e6, ageIntervalsArr))/1.0e9
print '\n\nTotal LMC stellar mass (x 10^9 M_sun) = ', tot_st_mass

binningScheme = [[0,1,2],[3,4,5],[6,7,8],[9,10],[11],[12],[13],[14,15]]
nBinsScheme = len(binningScheme)
sfhMapBinned = np.zeros((nBinsScheme, nCells))
sfhMapBinned_errlow = np.zeros((nBinsScheme, nCells))
sfhMapBinned_errhigh = np.zeros((nBinsScheme, nCells))
cell = 100
for i, bins in enumerate(binningScheme):
    sfhMapBinned[i, cell] = np.dot(sfhMapArr[cell, 0, bins]/1.0e6, ageIntervalsArr[bins])
    sfhMapBinned_errlow[i, cell] = np.dot(sfhMapArr[cell, 1, bins]/1.0e6, ageIntervalsArr[bins])
    sfhMapBinned_errhigh[i, cell] = np.dot(sfhMapArr[cell, 2, bins]/1.0e6, ageIntervalsArr[bins])

print '\n\nsfhMapBinned = ', sfhMapBinned[:, cell]
print 'Lower limit sfhMapBinned = ', sfhMapBinned_errlow[:, cell]
print 'Upper limit sfhMapBinned = ', sfhMapBinned_errhigh[:, cell]
if 0:
    print '\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n'
    print 'Median LMC stellar mass (x 10^9 M_sun) in cell ', cell, ' = ', np.sum(sfhMapBinned)/1.0e9
    print 'Error LMC stellar mass (x 10^9 M_sun) in cell ', cell, ' = ', np.sum(sfhMapBinned)/1.0e9 - np.sum(sfhMapBinned_errlow)/1.0e9
    print 'Error LMC stellar mass (x 10^9 M_sun) in cell ', cell, ' = ', np.sum(sfhMapBinned_errhigh)/1.0e9 - np.sum(sfhMapBinned)/1.0e9

if 0:
    print '\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n'
    print 'Median log M in cell ', cell, ' = ', sfhMapBinned[:, cell]
    print 'Lower limit log M in cell ', cell, ' = ', sfhMapBinned[:, cell] - sfhMapBinned_errlow[:, cell]
    print 'Upper limit log M in cell ', cell, ' = ', sfhMapBinned[:, cell] - sfhMapBinned_errhigh[:, cell]

