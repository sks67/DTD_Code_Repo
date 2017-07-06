import numpy as np                        ##IMPORTS
import numpy.ma as ma
import os
import glob
import string
import emcee
import corner
from scipy import misc
from scipy import stats
from scipy import interpolate
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy import constants as const
from astropy.io import fits
from astropy.io import ascii
import sys

DTDpath = '/Users/sumits2k/Desktop/Research/SNResearch2/RadioSNRs/DTD/'
pathName = DTDpath + 'Output_SFH_Files/'

#MCMC Parameters
nWalkers = 80
chainLength = 3000
nBurninSteps = 4000
print 'Number of walkers, chain length, burn-in length: ', nWalkers, chainLength, nBurninSteps
nIterations = 50
nParamsPerPage = 6

#Read SFHs and object count
objClassName = 'Cepheids'
obj_subtype = 'All'
binningScheme = 'Unbinned'
refName = 'OGLE'
galaxy = 'LMC'

#Function to plot chains 
def plotchains(sampler,nParamsPerPage,filePrefix,fileSuffix):

    nParams =  (sampler.chain.shape)[2] 
    nWalkers =  (sampler.chain.shape)[0] 
    
    plt.clf()
    plt.figure(figsize=(20, 20))
    plt.rcParams['font.size'] = 14
    plt.rcParams['legend.fontsize'] = 14
    nPages = np.ceil(nParams/nParamsPerPage)
    thisPage = 0

    for parameter in range(nParams):

        parName = '$\Psi_{' + str(parameter) + '}$'
    
        indexPlot = 2*(parameter-thisPage*nParamsPerPage)+1
        plt.subplot(nParamsPerPage, 2, indexPlot)
        plt.xlabel('Step')
        plt.ylabel(parName)
        for i in range(nWalkers): plt.semilogy(sampler.chain[i,:,parameter], color='k')
        
        indexPlot += 1 
        plt.subplot(nParamsPerPage, 2, indexPlot)
        plt.xlabel(parName)
        plt.ylabel('N')
        #if parameter == 0 : plt.xlim(0,1e-9)
        #else : plt.xlim(0,3e-10)
        plt.hist(sampler.flatchain[:,parameter],bins = 20)

        if (np.mod(parameter,nParamsPerPage) == nParamsPerPage-1) or (parameter == nParams-1) :
            pageTitle = '_Params_'+str(parameter-nParamsPerPage+1)+'_'+str(parameter)
            plt.savefig(filePrefix+pageTitle+fileSuffix)
            ## plt.savefig('/Users/badenes/python/MC_MCMC/MCMC_Error/Error_Random/MCMC_'
            ##             +galaxy+'_'+objClassName+'_'+binningScheme+'_Iter'+str(iteration).zfill(3)+pageTitle+'_Burn-in.png')
            thisPage = thisPage + 1
            plt.clf()
            plt.figure(figsize=(20, 20))

    #Plot triangle
    ## labels = ['$\Psi_'+s+'$' for s in map(str,range(nAgeBins))]
    ## figure = triangle.corner(sampler.flatchain,labels=labels)
    ## figure.savefig('/Users/badenes/python/MC_MCMC/MCMC_Error/Error_Random/MCMC_'+galaxy+'_'+objClassName+'_'+binningScheme+'_Iter'+str(iteration).zfill(3)+'_Triangle.png')
    return


if galaxy == 'LMC' :
    if refName == 'Boyer' :
        nCells = 1227
    elif refName == 'OGLE' :
        nCells = 808
    elif refName == 'Reid&Parker' :
        nCells = 723
    else :
        nCells = 1376

if galaxy == 'SMC' :
    if refName == 'Boyer' :
        nCells = 378
    else :
        nCells = 428

if binningScheme == 'Coarse' :
    nAgeBins = 3
elif binningScheme == 'Massive' :
    nAgeBins = 4
elif binningScheme == 'Submedium' :
    nAgeBins = 5
elif (binningScheme == 'Medium') or (binningScheme == 'MediumB') or (binningScheme == 'MediumC'):
    nAgeBins = 6
elif binningScheme == 'Subfine' :
    nAgeBins = 7
elif binningScheme == 'Fine' :
    nAgeBins = 8
elif binningScheme == 'CrazyCeph':
    nAgeBins = 8
elif binningScheme == 'Ultrafine' :
    nAgeBins = 11
elif binningScheme == 'Unbinned' :
    if galaxy == 'LMC' :
        nAgeBins = 16    
    elif galaxy == 'SMC' :
        nAgeBins = 18

objMap = np.zeros(nCells, dtype = np.int)
sfhMap = np.zeros((nAgeBins,nCells))
sfhMapMin = np.zeros((nAgeBins,nCells))
sfhMapMax = np.zeros((nAgeBins,nCells))
sfhMapRange = np.zeros((nAgeBins,nCells))
with open(pathName+galaxy+'_SFH_Cells_'+objClassName+obj_subtype+'_'+binningScheme+'.dat', 'r') as f:
    for cell in np.arange(nCells) :
        cellLineWords = string.split(f.readline().strip())
        objMap[cell] = float(cellLineWords[1])
        sfhMap[:,cell] = np.asarray(cellLineWords[2:]).astype(np.float)
        cellLineWords = string.split(f.readline().strip())
        sfhMapMin[:,cell] = np.asarray(cellLineWords).astype(np.float)
        cellLineWords = string.split(f.readline().strip())
        sfhMapMax[:,cell] = np.asarray(cellLineWords).astype(np.float)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Calculating Covariance matrix with total stellar mass
# and fraction of masses in adjacent bins 
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#/Users/sumits2k/Desktop/Research/SNResearch2/RadioSNRs/DTD/DTD_Code_Repo/SFH_Error_Diagnosis/sfHMap_output
np.savetxt(DTDpath+'DTD_Code_Repo/sfHMap_output/sfhMap.txt', sfhMap)
np.savetxt(DTDpath+'DTD_Code_Repo/sfHMap_output/sfhMapMin.txt', sfhMapMin)
np.savetxt(DTDpath+'DTD_Code_Repo/sfHMap_output/sfhMapMax.txt', sfhMapMax)

#column vector of total stellar mass per cells. Tested with a single cell
tot_mass_per_cell = np.sum(sfhMap, axis=0)  

#column vector of relative stellar mass fractions. Tested with single cells
f_bini_binj = np.array([massvect[:-1]/massvect[1:] for massvect in sfhMap.T])
 
sfhMap_newcoord = np.column_stack((tot_mass_per_cell, f_bini_binj)).T 

#Some elements will have nans or infinities because many age bins don't have star-formation. We can mask these
#values before doing covariance calculations
ma_sfhMap_newcoord = ma.masked_invalid(sfhMap_newcoord)  

#Lets output to a file so that we can share with Jeff
f = open('covariance_outputs.txt', 'w')
sys.stdout = f
print '\n\nCOVARIANCE MATRIX\n\n'
print np.ma.cov(ma_sfhMap_newcoord)
print '\n\nCORRELATION MATRIX\n\n'
print np.ma.corrcoef(ma_sfhMap_newcoord)
f.close()


