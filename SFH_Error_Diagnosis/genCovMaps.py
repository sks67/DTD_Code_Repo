#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Generate randomized maps from covariance matrices and display 
# them
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import sys

cell = input('Cell number: ')
path = '/Users/sumits2k/Desktop/Research/SNResearch2/RadioSNRs/DTD/DTD_Code_Repo/SFH_Error_Diagnosis/sfHMap_output/'
sfhMap = np.loadtxt(path+'sfhMap.txt')
sfhMapMin = np.loadtxt(path+'sfhMapMin.txt')
sfhMapMax = np.loadtxt(path+'sfhMapMax.txt')

log_sfhMap = np.log10(sfhMap)
log_sfhMapMin = np.log10(sfhMapMin)
log_sfhMapMax = np.log10(sfhMapMax)

log_sfhMap[np.isneginf(log_sfhMap)] = 0
log_sfhMapMin[np.isneginf(log_sfhMapMin)] = 0
log_sfhMapMax[np.isneginf(log_sfhMapMax)] = 0
if 0:
    ma_log_sfhMap = ma.masked_invalid(log_sfhMap)
    ma_log_sfhMapMin = ma.masked_invalid(log_sfhMapMin)
    ma_log_sfhMapMax = ma.masked_invalid(log_sfhMapMax)

np.set_printoptions(precision=1)



cov_sfh = np.cov(sfhMap)
#cov_log_sfh = np.ma.cov(ma_log_sfhMap)
#corr_log_sfh = np.ma.corrcoef(ma_log_sfhMap)

cov_log_sfh = np.cov(log_sfhMap)
corr_log_sfh = np.corrcoef(log_sfhMap)

print '\n\nCovariance matrix of full Star formation history map\n\n'
print np.around(cov_log_sfh, decimals=3)
print '\n\nCorrelation matrix of full SFH map\n\n'
print np.around(corr_log_sfh, decimals=3)
#print np.diagonal(cov_sfh)
print '\n\nLog stellar mass [Msun] per age bin for cell', cell, ' (Harris&Zaritsky 2009) :-'

print '\n\nUpper limit:', log_sfhMapMax[:, cell]
print 'Median     :', log_sfhMap[:, cell]
print 'Lower limit:', log_sfhMapMin[:, cell]

mean_sigma_tot = np.maximum(log_sfhMapMax[:, cell] - log_sfhMap[:, cell], log_sfhMap[:, cell] - log_sfhMapMin[:, cell])
sum_offdiag_rows = np.sum(cov_log_sfh, axis=0) - np.diagonal(cov_log_sfh)
sigma_rand = mean_sigma_tot**2.0 - sum_offdiag_rows

print '\n\nTotal error (random + covariant) from above'
print '(Note: I assume error = max(upper limit - median, median - lower limit). So this is the max error per stellar mass\n'
print mean_sigma_tot

print '\n\nCovariant errors: Sum of off-diagonal elements per row in covariance matrix\n'
print sum_offdiag_rows

print '\n\nVariance = total error^2 - covariant errors\n'
print sigma_rand





