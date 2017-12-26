import numpy as np                        ##IMPORTS
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
import dtdplotutils as dtdplot
reload(dtdplot)

DTDpath = '/Users/sumits2k/Desktop/Research/SNResearch2/RadioSNRs/DTD/'
pathName = DTDpath + 'Output_SFH_Files/'


#MCMC Parameters
nWalkers = 80
chainLength = 3000
nBurninSteps = 4000
print 'Number of walkers, chain length, burn-in length: ', nWalkers, chainLength, nBurninSteps
nIterations = 0
nParamsPerPage = 6

#Read SFHs and object count
objClassName = 'RRLyrae'
obj_subtype = 'All'
binningScheme = 'Unbinned'
refName = 'OGLE'
galaxy = 'LMC'
checkSFH = 'n'#raw_input('Save SFH comparisons? [y/n]: ')
subclassFolder = 'RRLyrae_linearSFH'
file_special_prefix = 'Fake_' #If, e.g., you're doing the fake maps. If not, then set galaxy_pref = ''

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
with open(pathName+file_special_prefix+galaxy+'_SFH_Cells_'+objClassName+obj_subtype+'_'+binningScheme+'.dat', 'r') as f:
    for cell in np.arange(nCells) :
        cellLineWords = string.split(f.readline().strip())
        objMap[cell] = float(cellLineWords[1])
        sfhMap[:,cell] = np.asarray(cellLineWords[2:]).astype(np.float)
        cellLineWords = string.split(f.readline().strip())
        sfhMapMin[:,cell] = np.asarray(cellLineWords).astype(np.float)
        cellLineWords = string.split(f.readline().strip())
        sfhMapMax[:,cell] = np.asarray(cellLineWords).astype(np.float)

print 'Objects found:', sum(objMap)

#Print total masses formed
## for timeBin in range(nAgeBins) :
##     print 'Bin %d: %10.2e %10.2e %10.2e' % (timeBin, np.sum(sfhMap[timeBin,:]), np.sum(sfhMapMin[timeBin,:]), np.sum(sfhMapMax[timeBin,:]))


#Calculate array of factorials for use in likelihoods
objMapFact = misc.factorial(objMap)
normThreshold = 25 #Threshold to use Normal instead of Poisson
if objMap.max() >= normThreshold :
    normApp = True 
    normIndexes = np.where(objMap >= normThreshold)
    poissIndexes = np.where(objMap < normThreshold)
    print 'Using Normal approximation for ', len(normIndexes[0]), ' cells'
else :
    normApp = False
    normIndexes = np.array([])
    poissIndexes = np.arange(nCells)
    print 'Using Poisson statistics'

#Define the likelihood function and the priors
def ln_like(psi, sfhMap, objMap, objMapFact, normApp, normIndexes, poissIndexes):
    l = np.dot(psi,sfhMap)
    if not normApp :
        return -1.0 * sum(l) + sum(np.log((l**objMap)/objMapFact)) 
    else :
        likeArr = np.zeros(nCells) 
        likeArr[normIndexes] = -0.5*np.log(2.0*np.pi*l) - ((objMap-l)**2.0)/(2.0*l)
        likeArr[poissIndexes] = -1.0 *l + np.log((l**objMap)/objMapFact)
        return sum(likeArr)    
 
def ln_prior(psi):
    if (psi < 0.0).any() : return -np.inf
    else : return 0.0

def ln_prob(psi, sfhMap, objMap, objMapFact, normApp, normIndexes, poissIndexes):
    lnprior = ln_prior(psi)
    if not np.isfinite(lnprior) : 
        return lnprior 
    return lnprior + ln_like(psi, sfhMap, objMap, objMapFact, normApp, normIndexes, poissIndexes)
    ## if np.isnan(ln_like(psi, sfhMap, objMap, objMapFact, normApp, normIndexes, poissIndexes)).any() : 
    ##     print 'Oh-oh'
    ##     pdb.set_trace()
    #return np.nansum(ln_prior(psi),ln_like(psi, sfhMap, objMap, objMapFact, normApp, normIndexes, poissIndexes))

#Exploration of SFH errors
randSfhMap = np.zeros((nAgeBins,nCells,nIterations)) 
print 'Randomizing SFH'
for cell in np.arange(nCells) :
    for ageBin in np.arange(nAgeBins) :
        if sfhMap[ageBin,cell] == 0.0 : #Best-fit is 0
            sfhBinCellRand = np.random.normal(0.0,sfhMapMax[ageBin,cell],nIterations)
            sfhBinCellRand[sfhBinCellRand < 0.0] = 0.0
            randSfhMap[ageBin,cell,:] = sfhBinCellRand
            
        ## elif sfhMapMin[ageBin,cell] == 0.0 : #Lower limit is 0
        ##     sfhBinCellRand = np.random.normal(sfhMap[ageBin,cell],sfhMapMax[ageBin,cell]-sfhMap[ageBin,cell],nIterations)
        ##     sfhBinCellRand[sfhBinCellRand < 0.0] = 0.0
        ##     randSfhMap[ageBin,cell,:] = sfhBinCellRand
            
        else :  #Best-fit is NOT 0
            sfhBinCellRand = np.random.normal(sfhMap[ageBin,cell],sfhMap[ageBin,cell]-sfhMapMin[ageBin,cell],nIterations)
            sfhBinCellRand[sfhBinCellRand < 0.0] = 0.0
            randSfhMap[ageBin,cell,:] = sfhBinCellRand
            ## low = sfhMapMin[ageBin,cell]
            ## mean = sfhMap[ageBin,cell]
            ## high = sfhMapMax[ageBin,cell]
            ## lowGauss = np.random.normal(mean,(mean-low),3*nIterations)
            ## lowGauss[lowGauss < 0.0] = 0.0
            ## lowGaussClip = lowGauss[lowGauss < mean]
            ## highGauss = np.random.normal(mean,high-mean,3*nIterations)
            ## highGaussClip = highGauss[highGauss > mean]
            ## frankenGauss = np.concatenate([lowGaussClip[0:nIterations],highGaussClip[0:nIterations]])
            ## np.random.shuffle(frankenGauss)
            ## randSfhMap[ageBin,cell,:] = frankenGauss[0:nIterations]

#Save randomized and nominal maps
#np.save('/Users/badenes/python/MC_MCMC/MCMC_Error/Error_Random/'+galaxy+'_'+refName+'_'+binningScheme+'_randomSFHMap',randSfhMap) 
#np.save('/Users/badenes/python/MC_MCMC/MCMC_Error/Error_Random/'+galaxy+'_'+refName+'_'+binningScheme+'_SFHMap',sfhMap) 
#np.save('/Users/badenes/python/MC_MCMC/MCMC_Error/Error_Random/'+galaxy+'_'+refName+'_'+binningScheme+'_SFHMapMin',sfhMapMin) 
#np.save('/Users/badenes/python/MC_MCMC/MCMC_Error/Error_Random/'+galaxy+'_'+refName+'_'+binningScheme+'_SFHMapMax',sfhMapMax) 

print 'Running ', nIterations, ' iterations.'
if checkSFH == 'y':
    dtdplot.check_randSFH(sfhMap, sfhMapMax, sfhMapMin, randSfhMap, binning='Unbinned', fileSuffix = subclassFolder,  showPlot = False)
    print 'Plot saved in DTD_Plots/'

#sys.exit()
print randSfhMap.shape
for iteration in np.arange(-3, nIterations) :

    if iteration == -3 :
        iterName = '_Nominal'
        sfhMapIter = sfhMap
        
    elif iteration == -2 :
        iterName = '_LowLim'
        sfhMapIter = sfhMapMin
        
    elif iteration == -1 :
        iterName = '_HighLim'
        sfhMapIter = sfhMapMax
        
    else :
      iterName = '_Iter'+str(iteration).zfill(3) 
      sfhMapIter = randSfhMap[:,:,iteration]
    
    print 'Iteration ', iterName
    filePrefix = DTDpath + 'Chains_Burnins_pngs/'+galaxy+'_'+objClassName+obj_subtype+'_'+binningScheme+iterName
    
    #Choose an initial set of DTDs for the walkers; normal in log around -10, sigma = 3
    log10PsiMean0 = -5.0
    log10PsiSigma0 = 0.1
    log10Psi0 = [log10PsiSigma0*np.random.randn(nAgeBins)+log10PsiMean0 for i in xrange(nWalkers)]
    psi0 = np.power(10.0,log10Psi0)

    ## #Choose randomized sfh map
    ## if iteration == 0 :
    ##     sfhMapIter = sfhMap 
    ## elif iteration == 1 :
    ##     sfhMapIter = sfhMapMin
    ## elif iteration == 2 :
    ##     sfhMapIter = sfhMapMax
    ## else :
   
    
    #Initialize the sampler with the chosen specs.
    print 'Initializing sampler'
    sampler = emcee.EnsembleSampler(nWalkers, nAgeBins, ln_prob, args=[sfhMapIter, objMap, objMapFact, normApp, normIndexes, poissIndexes])
    pos, prob, state = sampler.run_mcmc(psi0, nBurninSteps)

    #Plot burn-in 
#    plotchains(sampler,nParamsPerPage,filePrefix,'_Burn-in.png')         

    #Check that burn-in went well
    afb = np.mean(sampler.acceptance_fraction)
    afbThreshold = 0.01
    print 'Mean acceptance fraction during Burn-in:', afb
    if afb < afbThreshold :
        print 'Something went wrong'
        if step == nSteps-1 : sys.exit('End of loop.')
    else :
        print 'Success!'
        
    # Reset the chain to remove the burn-in samples.
    print 'Resetting chain'
    sampler.reset()

    # Starting from the final position in the burn-in chain, run sampler
    print 'Starting run'
#    stop
    sampler.run_mcmc(pos, chainLength, rstate0=state)

    # Print out the mean acceptance fraction.
    print 'Mean acceptance fraction:', np.mean(sampler.acceptance_fraction)

    #Store the chains to a FITS file: one table extension per parameter
    print 'Writing chains to FITS file'
    hdu = fits.PrimaryHDU(sampler.chain)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(DTDpath + file_special_prefix+'MCMC_DTD_fits/DTD_'+objClassName+'/'+ subclassFolder+'/'+galaxy+'_MCMC_DTD_'+objClassName+obj_subtype+'_'+binningScheme+iterName+'.fits', clobber = True)  
    hdulist.close()

    #Plot chains 
#    plotchains(sampler,nParamsPerPage,filePrefix+subclassFolder,'.png')
