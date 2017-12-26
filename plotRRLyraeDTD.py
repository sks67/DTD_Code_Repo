import numpy as np
import matplotlib.pyplot as plt
import corner
from astropy.io import fits
import dtdutils
reload(dtdutils)
import dtdplotutils as plot_dtd
reload(plot_dtd)
from scipy import interpolate
from scipy import stats
from astropy.io import ascii
from scipy import integrate
import sys

DTDpath = '/Users/sumits2k/Desktop/Research/SNResearch2/RadioSNRs/DTD/'

#~~~~~~~~ INPUT ~~~~~~~~~~~~~~~~~~~#
file_special_prefix = ''#'Fake_'
outfilePrefix = 'DTD_Plots/'
objName = 'RRLyrae'
obj_subtype_arr = ['All']
subfolderName = 'RRLyrae_linearSFH/'
filePrefix = DTDpath + file_special_prefix+'MCMC_DTD_fits/DTD_'+objName+'/' + subfolderName 
sfhBinType = 'Unbinned'
nIterations = 100
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if 1:
    masterchain_randSFH = []
    print 'Creating Master Chain...'
    for i in range(nIterations):
        if i%10 == 0.:
            print i
        chain_randSFH, header = fits.getdata(filePrefix + 'LMC_MCMC_DTD_' + objName + obj_subtype_arr[0] + '_' + sfhBinType + \
                                             '_Iter' + plot_dtd.string_filenum(i) + '.fits', 0, header=True)
        masterchain_randSFH.append(chain_randSFH)

    masterchain_randSFH = np.asarray(masterchain_randSFH)
    
nIterations = masterchain_randSFH.shape[0]

chain_nom_rrl = []
chain_low = []
chain_high = []
print 'Reading subtypes : ', obj_subtype_arr
for obj_subtype in obj_subtype_arr:
    x_nom, x_low, x_high = plot_dtd.read_chains(filePrefix=filePrefix, objName=objName, \
                                                    obj_subtype=obj_subtype, sfhBinType=sfhBinType)
    chain_nom_rrl.append(x_nom)
    chain_low.append(x_low)
    chain_high.append(x_high)

    
#Flatten the chains into a 2D array, with rows equal to #chain elements
#and columns = #parameters

print 'Flatten and average chains....'

ndim = masterchain_randSFH.shape[3]
chain_randSFH = masterchain_randSFH.reshape((-1, ndim))

#Choosing average chain from the superchain
nParams = int(chain_randSFH.shape[1])  #Number of DTD parameters
chainLength = int(chain_randSFH.shape[0]/nIterations*1.0)
print 'Length of chains : ', chainLength#Length of each DTD Markov Chain
chain_SFH = np.zeros((int(nParams), chainLength))
for par in range(nParams):
    chain_SFH[par] = np.random.choice(chain_randSFH[:,par], size=chainLength, replace=False)
    #You want the size of the randomly drawn chain to be the same size as the other chain, to 
    #normalize the statistical comparisons.
    
chain_SFH_rrl = chain_SFH.T

#plot_dtd.errorTable(chain_nom_rrl, chain_low, chain_high, [chain_SFH_rrl], objName, obj_subtype_arr, sfhBinType)

#trueDTD = np.concatenate(([0.]*14, [1.0e-5]*2))
trueDTD = np.array([])
if 1:
    plot_dtd.plot_dtds_SFH(chain_nom_rrl, [chain_SFH_rrl], trueDTD, objName, obj_subtype_arr, sfhBinType = [sfhBinType], mult_factor=1.0, \
                               colorScheme = 'bmh', stat_sad_crit = False, show_errors = True, show_uplims = True, savefigure = False, plot_title='OGLE III Result - Statistical', ylabel='RR Lyrae DTD', savefile = DTDpath + 'Writeup/'+objName+file_special_prefix+'_DTD_DetSigs_Statistical.pdf')

if 0:
    plot_dtd.plot_dtds_AllErrors(chain_nom_rrl, [chain_SFH_rrl], trueDTD, objName, obj_subtype_arr, sfhBinType = [sfhBinType], mult_factor=1.0, \
          colorScheme = 'bmh', stat_sad_crit = False, show_errors = True, show_uplims = True, savefigure = False, \
              savefile = DTDpath + 'Writeup/'+objName+file_special_prefix+'_DTD_noOutliers_Cepheids.pdf')
