#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Will check if the randomized star-formation histories have similar
# errors if generated :-
# a) from the covariance matrix.
# b) from treating bins as independent, and variance of each bin
#    being the sum of the rows of the covariance matrix
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

import numpy as np
import matplotlib.pyplot as plt
import sys
import dtdutils as dtd

path = '/Users/sumits2k/Desktop/Research/SNResearch2/RadioSNRs/DTD/DTD_Code_Repo/sfHMap_output/'
sfhMap = np.loadtxt(path+'sfhMap.txt')
sfhMapMin = np.loadtxt(path+'sfhMapMin.txt')
sfhMapMax = np.loadtxt(path+'sfhMapMax.txt')

M_fullMap = np.sum(sfhMap, axis=1) 
M_fullMapMin = np.sum(sfhMapMin, axis=1) 
M_fullMapMax = np.sum(sfhMapMax, axis=1) 

#Stellar mass per age-bin with errors for the full LMC-Cepheid survey area
logM_fullMap = np.log10(M_fullMap)
logM_fullMapMin = np.log10(M_fullMapMin)
logM_fullMapMax = np.log10(M_fullMapMax)

print 'Total stellar mass (x 10^9 Msun) = ', np.sum(M_fullMap)/1.0e9

#********************************************************#
#       COVARIANCE MATRIX FOR FULL LMC-CEPHEID AREA      #
#********************************************************#

#Even though log space has infinities, we can zero them. This means the bins with 0 Msun, have 1 Msun, which I believe                                  
#low enough as well on the astrophysical mass scale.                                                                                                       

log_sfhMap = np.log10(sfhMap)
log_sfhMapMin = np.log10(sfhMapMin)
log_sfhMapMax = np.log10(sfhMapMax)

log_sfhMap[np.isneginf(log_sfhMap)] = 0
log_sfhMapMin[np.isneginf(log_sfhMapMin)] = 0
log_sfhMapMax[np.isneginf(log_sfhMapMax)] = 0

cov_log_sfh = np.cov(log_sfhMap) + 0.3
corr_log_sfh = np.corrcoef(log_sfhMap)

print '\n\nCovariance matrix of full Star formation history map\n\n'
print np.around(cov_log_sfh, decimals=3)
print '\n\nCorrelation matrix of full SFH map\n\n'
print np.around(corr_log_sfh, decimals=3)

#********************************************************#
#       RANDOMIZED STAR FORMATION HISTORIES              #
#********************************************************#

massVect_covar = np.random.multivariate_normal(logM_fullMap, cov_log_sfh, size=500).T
meanM = np.mean(massVect_covar, axis=1)
stdM = np.std(massVect_covar, axis=1)
varM = np.var(massVect_covar, axis=1)

#Now assume each bin is indepedent, error in each bin is the variance, no covariance
cov_log_sfh_indep = np.diag(varM)

massVect_indep = np.random.multivariate_normal(logM_fullMap, cov_log_sfh_indep, size=500).T
meanM_indep = np.mean(massVect_indep, axis=1)
stdM_indep = np.std(massVect_indep, axis=1)

if 0: #Here I was summing up the rows of the covariance matrix. On second thought, I don't think thats what the paper means.
    tot_err_log_sfh = np.sum(cov_log_sfh, axis=1)
    cov_log_sfh_onlyvar = np.diag(tot_err_log_sfh)
    massVect_onlyvar = np.random.multivariate_normal(logM_fullMap, cov_log_sfh_onlyvar, size=500).T
    meanM_onlyvar = np.mean(massVect_onlyvar, axis=1)
    stdM_onlyvar = np.std(massVect_onlyvar, axis=1)
    
#********************************************************#
#       STAR FORMATION HISTORY PLOT                      #
#********************************************************#

ages, age_left, age_right = dtd.sfh_ageBins('Unbinned')
age_left_err = np.log10(ages) - np.log10(age_left)
age_right_err = np.log10(age_right) - np.log10(ages)
plt.figure(figsize=(7,6))
plt.title('Error simulation for Full SFH Map', fontsize=18)
#plt.errorbar(np.log10(ages), logM_fullMap, xerr=(age_left_err, age_right_err), yerr=(logM_fullMapMax - logM_fullMap, logM_fullMap - logM_fullMapMin), fmt='o', color='k', ecolor='k', elinewidth=3.0, label='Original')
plt.errorbar(np.log10(ages), meanM, xerr=(age_left_err, age_right_err),  yerr = stdM, ls='', fmt = '', ecolor = 'k', label='Covariant errors')
plt.errorbar(np.log10(ages), meanM_indep, ls='', yerr = stdM_indep, fmt = '', ecolor = 'r', elinewidth=3.0, alpha=0.3, label='Independent errors (var+covar)')
plt.xlabel('Log Age [yr]', fontsize=15)
plt.ylabel(r'Log M$_{*}$ [M$_{\odot}$]', fontsize=15)
plt.legend(loc=2)

plt.tight_layout()
plt.savefig(path + 'ErrorSimulation.png')
plt.show()

