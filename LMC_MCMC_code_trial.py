import numpy as np
import matplotlib.pyplot as plt
import sys
import dtdutils as dtd
from scipy.stats import truncnorm
import dtdplotutils as dtdplot
reload(dtdplot)

nAgeBins = 16
nCells = 20
nIterations = 200

checkSFH = raw_input('Check if random SFHs are consistent with HZ09? [y/n]: ')

DTDpath = '/Users/sumits2k/Desktop/Research/SNResearch2/RadioSNRs/DTD/'
sfhMap = np.loadtxt(DTDpath+'DTD_Code_Repo/sfHMap_output/sfhMap.txt')
sfhMapMin = np.loadtxt(DTDpath+'DTD_Code_Repo/sfHMap_output/sfhMapMin.txt')
sfhMapMax = np.loadtxt(DTDpath+'DTD_Code_Repo/sfHMap_output/sfhMapMax.txt')

log_sfhMap = np.log10(sfhMap)
log_sfhMapMin = np.log10(sfhMapMin)
log_sfhMapMax = np.log10(sfhMapMax)

log_sfhMap[np.isneginf(log_sfhMap)] = 0
log_sfhMapMin[np.isneginf(log_sfhMapMin)] = 0
log_sfhMapMax[np.isneginf(log_sfhMapMax)] = 0



randSfhMap_maxerr = np.zeros((nAgeBins,nCells,nIterations))
print 'Randomizing SFH'
for cell in np.arange(nCells):
    #Max error method
    error_logM = np.maximum(log_sfhMapMax[:, cell] - log_sfhMap[:, cell], log_sfhMap[:, cell] - log_sfhMapMin[:, cell])
    cov_logM_cell = np.diag(error_logM**2.0)
    sfhBinCellRand = np.random.multivariate_normal(log_sfhMap[:, cell], cov_logM_cell, size = nIterations).T
#    sfhBinCellRand[sfhBinCellRand < 0.0] = 0.0
    randSfhMap_maxerr[:, cell, :] = sfhBinCellRand

if checkSFH == 'y':
    dtdplot.check_randSFH(log_sfhMap, log_sfhMapMax, log_sfhMapMin, randSfhMap_maxerr, binning='Unbinned', fileSuffix='max',  showPlot = True)
    print 'Plot saved in DTD_Plots/'


        
