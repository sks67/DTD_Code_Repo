import numpy as np
import matplotlib.pyplot as plt
import os

DTDpath = os.getenv('DTD')
plotPath = DTDpath + '/Writeup/'

ages_edges = np.array([4.5, 8.9, 17.8, 35.5, 70.8, 125.9, 199.5, 316.2, 501.2, 794.3, 1300., 2000., 3200., 5000., 7900., 12600, 20000.])

log_age_cent = np.log10(ages_edges[:-1]) + 0.5*(np.log10(ages_edges[1:]) - np.log10(ages_edges[:-1]))
ages_cent = 10.0**log_age_cent

TRR_nondetect = np.array([0.25, 0.11, 0.19, 0.08, 0.06, 0.52, 0.28, 0.74, 0.58, 0.88, 0., 0., 0., 0., 0., 0.])*1.0e6
 
TRR_detect = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.72, 1.61, 1.84, 3.0, 1.72, 2.54])*1.0e6

TRR_error_low = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.12, 0.10, 0.25, 0.37, 0.39, 0.19])*1.0e6

TRR_error_up = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.15, 0.09, 0.25, 0.27, 0.3, 0.18])*1.0e6

show_errors = False
show_uplims = False

#BMH color scheme
clr = [u'#348ABD', u'#A60628',u'#467821',  u'#D55E00',u'#CC79A7',\
                            u'#56B4E9', u'#009E73', u'#F0E442', u'#0072B2']
with plt.style.context(('bmh')):
            #This temporarily sets the plotting style. Do not use plt.style.use()                                                                              
    x_axis = np.vstack((ages_edges[:-1], ages_edges[1:])).T.flatten()
    y_axis = np.repeat(TRR_detect,2)
    y_axis_2 = np.repeat(TRR_nondetect, 2)
    print x_axis.shape
    print y_axis.shape
   
    plt.plot(x_axis, y_axis, color = clr[0], ls = '-', lw=2.0)
    plt.errorbar(ages_cent, TRR_detect, yerr=(TRR_error_low, TRR_error_up), fmt='', ls='', ecolor='k')
    plt.plot(x_axis, y_axis_2, color = clr[0], ls = '--', lw=2.0)
    plt.gca().fill_between(x_axis, 0., y_axis, facecolor=clr, alpha=0.8)

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Delay Time [Myrs]', fontsize=15)
plt.xlim(5.0, 2.0e4)
plt.ylim(5.0e4,1.0e7)                                                                                                                                   
plt.ylabel(r'RR Lyrae Lifetimes, $T_{RR}$ [yr]', fontsize=15)
#plt.text(700, 2.0e-4, 'RR-Lyrae DTD in LMC\n(OGLE III + \nZaritsky \'09 SFH maps)')
plt.tick_params(axis='y', labelsize=14)
plt.tick_params(axis='x', labelsize=14)
plt.legend(loc=2, numpoints=1)
plt.savefig(plotPath + 'RRLyraeLifetimes.pdf', dpi=100)
plt.show()
