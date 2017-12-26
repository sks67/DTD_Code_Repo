import numpy as np
import os
import string
import emcee
import corner
from scipy import misc
from matplotlib import pyplot as plt
import sys
import dtdutils

DTD_path = os.getenv('DTD')+'/'
sad_files_path = DTD_path + 'Output_SFH_Files/'


def nCells_of_Survey(galaxyName, surveyName):
    nCells_dict = {'LMC':{'Boyer':1227, 'OGLE':808, 'Reid&Parker':723, '':1376},                    'SMC':{'Boyer':378, '':428}}
    return nCells_dict[galaxyName][surveyName]

def nAgebins_scheme(binningScheme):
    bin_dict = {'Coarse':3, 'Massive':4, 'Medium':6, 'MediumB':6, 'MediumC':6, 'Subfine':7, 'Fine':8,               'CrazyCeph':8, 'Ultrafine':11, 'Unbinned': 16}
    return bin_dict[binningScheme]

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


# ### DTD Recovery Code

# In[13]:

object_Name = 'RRLyrae'
object_Subtype = 'All'
binning_type = 'Unbinned'
#MCMC Parameters
nWalkers = 80
chainLength = 3000
nBurninSteps = 4000
print 'Number of walkers, chain length, burn-in length: ', nWalkers, chainLength, nBurninSteps
nIterations = 0
normThreshold = 25 #Threshold to use Normal instead of Poisson
afbThreshold = 0.01

sad_file_name = sad_files_path + 'LMC_SFH_Cells_' + object_Name + object_Subtype + '_' + binning_type + '.dat'

#Open file containing stellar masses for a given binning strategy, and only use the 
#rows with best-fit stellar masses. 
sad_file_obj = open(sad_file_name, 'r')
sad_file = sad_file_obj.readlines()
sad_file_obj.close()

sad_cellnames = np.array([lines.split()[0] for lines in sad_file if len(lines.split())==18])
sads = np.array([map(float, lines.split()[2:]) if len(lines.split())==18 else map(float, lines.split())                  for lines in sad_file])

#Lots of things happening here. I'm reading each line of the file with the for statement, picking only lines that have
#best fit stellar mass (i.e. with 18 entries, 1st two being cell name and number of objects in cell). The split() splits
#up the full file string, and map(float,...) converts all the SAD entries into numbers!
dtd_old = np.concatenate(([0.]*14, [1.0e-5]*2))
    #Read SFHs and object count
objClassName = 'RRLyrae'
obj_subtype = 'All'
binningScheme = 'Unbinned'
refName = 'OGLE'
galaxy = 'LMC'
checkSFH = 'n'#raw_input('Save SFH comparisons? [y/n]: ')
subclassFolder = 'RRLyrae_linearSFH'
galaxy_pref = 'Fake_' #If, e.g., you're doing the fake maps. If not, then set galaxy_pref = ''

nCells = nCells_of_Survey(galaxy, refName)
nAgeBins = nAgebins_scheme(binningScheme)
sfhMap = sads[::3].T
sfhMapMin = sads[1::3].T
sfhMapMax = sads[2::3].T
ages, agebin_1, agebin_2 = dtdutils.sfh_ageBins(binning_type)

dtd_signal_array = []
dtd_mode_array = []
#plt.figure(figsize=(6,5))
nIterations = 100
for i in range(nIterations):
    print 'iteration ', i
    lambda_i = np.dot(sads[0::3], dtd_old)
    #print 100.*np.sum(lambda_i)/23459.  #Check if the fraction of old objects produced is consistent.
    objMap = np.random.poisson(lam=lambda_i)

    #Calculate array of factorials for use in likelihoods
    objMapFact = misc.factorial(objMap)
    if objMap.max() >= normThreshold :
        normApp = True 
        normIndexes = np.where(objMap >= normThreshold)
        poissIndexes = np.where(objMap < normThreshold)
     #   print 'Using Normal approximation for ', len(normIndexes[0]), ' cells'
    else :
        normApp = False
        normIndexes = np.array([])
        poissIndexes = np.arange(nCells)
    #    print 'Using Poisson statistics'

    #Choose an initial set of DTDs for the walkers; normal in log around -10, sigma = 3
    log10PsiMean0 = -5.0
    log10PsiSigma0 = 0.1
    log10Psi0 = [log10PsiSigma0*np.random.randn(nAgeBins)+log10PsiMean0 for i in xrange(nWalkers)]
    psi0 = np.power(10.0,log10Psi0)   

    #Initialize the sampler with the chosen specs.
#    print 'Initializing sampler'
    sampler = emcee.EnsembleSampler(nWalkers, nAgeBins, ln_prob, args=[sfhMap, objMap, objMapFact, normApp, normIndexes, poissIndexes])
    pos, prob, state = sampler.run_mcmc(psi0, nBurninSteps)         

    #Check that burn-in went well
    afb = np.mean(sampler.acceptance_fraction)
#    print 'Mean acceptance fraction during Burn-in:', afb
    if afb < afbThreshold :
        print 'Something went wrong'
        if step == nSteps-1 : sys.exit('End of loop.')
    #else :
        #continue
        #print 'Success!'

    # Reset the chain to remove the burn-in samples.
#    print 'Resetting chain'
    sampler.reset()

    # Starting from the final position in the burn-in chain, run sampler
#    print 'Starting run'
    #    stop
    sampler.run_mcmc(pos, chainLength, rstate0=state)

    # Print out the mean acceptance fraction.
#    print 'Mean acceptance fraction:', np.mean(sampler.acceptance_fraction)

    dtdchains = sampler.flatchain 

    dtd, dtd_error, isuplim = dtdutils.dtd_staterrors(dtdchains)
    dtd_no = np.array([signal if isuplim[i] else 0. for i, signal in enumerate(dtd)])
    dtd_yes = np.array([signal  if np.invert(isuplim[i]) else 0. for i,signal in enumerate(dtd)])
    dtd_signal_array.append(dtd_yes)
    #Appending all the DTD modes
    dtd_mode_array.append(np.array([dtdutils.dtd_mode(dtd_j) for dtd_j in dtdchains.T]))
    dtd_y = np.array([signal for i, signal in enumerate(dtd) if np.invert(isuplim[i])])
    ages_y = np.array([age for i, age in enumerate(ages) if np.invert(isuplim[i])])
#    dtd_yes_error = np.array([err for i, err in enumerate(dtd_error) if np.invert(isuplim[i])])
#    dtd_yes_error_low = dtd_yes_error[:,0]
#    dtd_yes_error_high = dtd_yes_error[:,1]
    
#    plt.errorbar(ages_y/1.0e6, dtd_y, yerr=(dtd_yes_error_low, dtd_yes_error_high), \
   #                  uplims=False, fmt='o', ls='', ecolor='k', label='')


# In[14]:

np.count_nonzero(dtd_signal_array, axis=0)/100.
np.savetxt('')


# In[15]:




# In[6]:

fig = corner.corner(dtdchains*1.0e6, quantiles=[0.16, 0.5, 0.84])


# In[ ]:



