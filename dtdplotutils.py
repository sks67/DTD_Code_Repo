from astropy.io import fits
import dtdutils
reload(dtdutils)
import numpy as np
import matplotlib.pyplot as plt



#All the data are in units of PNe/1.0e6 solar masses.
def read_chains(filePrefix='', objName='', obj_subtype='All', sfhBinType=''):
    #Returns the nominal, high, low chains
    data_nominal, header = fits.getdata(filePrefix+'LMC_MCMC_DTD_'+objName+obj_subtype+'_'+sfhBinType+'_Nominal.fits', 0, header=True)
    data_low, header = fits.getdata(filePrefix+'LMC_MCMC_DTD_'+objName+obj_subtype+'_'+sfhBinType+'_Lowlim.fits', 0, header=True)
    data_high, header = fits.getdata(filePrefix+'LMC_MCMC_DTD_'+objName+obj_subtype+'_'+sfhBinType+'_Highlim.fits', 0, header=True)
    ndim = data_nominal.shape[2]
    data_nominal = data_nominal.reshape((-1, ndim))
    data_low = data_low.reshape((-1, ndim))
    data_high = data_high.reshape((-1, ndim))
    return data_nominal, data_low, data_high

def string_filenum(ind):
    #Returns the file number of randomized SFH DTDs in 00x format
    return ''.join(['0' for i in range(3 - len(str(ind)))])+str(ind)

    
def plot_dtds(data_nom_arr, data_low_arr, data_high_arr, objName, obj_subtype, sfhBinType, \
              colorScheme='bmh', show_errors = False):
    
    color_options = {'bmh':[u'#348ABD', u'#A60628', u'#7A68A6', u'#467821', u'#D55E00',u'#CC79A7', u'#56B4E9', u'#009E73', u'#F0E442', u'#0072B2'],
                    'seaborn-colorblind':[u'#0072B2', u'#009E73', u'#D55E00', u'#CC79A7', u'#F0E442', u'#56B4E9']}
    colors = iter(color_options[colorScheme])
    labels = iter(obj_subtype)
    
    plt.rc('font', family='serif')
    plt.figure(figsize=(6,5))
    for data_nominal, data_low, data_high in zip(data_nom_arr, data_low_arr, data_high_arr):
        
        ages, agebin_1, agebin_2 = dtdutils.sfh_ageBins(sfhBinType)

        #Extracting the signals (or lack thereof) with statistical errors
        dtd, dtd_error, isuplim = dtdutils.dtd_staterrors(data_nominal*1.0e6)

        #Extracting the systematic errors due to different SFHs
        dtd_sad_error_low, dtd_sad_error_high = dtdutils.dtd_saderrors(data_low*1.0e6, 
                                                                      data_nominal*1.0e6, data_high*1.0e6)


        #Separating bins with upper limits
        dtd_no      = np.array([np.array([signal, int(i)]) for i, signal in enumerate(dtd) if isuplim[i]])
        ages_no     = ages[dtd_no[:,1].astype(int)]
        agebin_1_no = agebin_1[dtd_no[:,1].astype(int)]
        agebin_2_no = agebin_2[dtd_no[:,1].astype(int)]

        #Separating bins with signals
        dtd_yes = np.array([np.array([signal,i]) for i, signal in enumerate(dtd) if np.invert(isuplim[i])])
        dtd_yes_error = np.array([err for i, err in enumerate(dtd_error) if np.invert(isuplim[i])])
        dtd_yes_error_low = dtd_yes_error[:,0]
        dtd_yes_error_high = dtd_yes_error[:,1]
        ages_yes     = ages[dtd_yes[:,1].astype(int)]
        agebin_1_yes = agebin_1[dtd_yes[:,1].astype(int)]
        agebin_2_yes = agebin_2[dtd_yes[:,1].astype(int)]

        #Total error
        dtd_yes_sad_error_low = dtd_sad_error_low[dtd_yes[:,1].astype(int)]
        dtd_yes_sad_error_high = dtd_sad_error_high[dtd_yes[:,1].astype(int)]
        dtd_yes_toterr_low = np.sqrt(dtd_yes_error_low**2 + dtd_yes_sad_error_low**2)
        dtd_yes_toterr_high = np.sqrt(dtd_yes_error_high**2 + dtd_yes_sad_error_high**2)

        with plt.style.context((colorScheme)):
            #This temporarily sets the plotting style. Do not use plt.style.use()
            clr = colors.next()
            plt.plot(np.vstack((agebin_1_yes/1.0e6, agebin_2_yes/1.0e6)).T.flatten(), np.repeat(dtd_yes[:,0],2),\
                     color=clr, ls='-',lw=2.0, label=labels.next())
            plt.plot(np.vstack((agebin_1_no/1.0e6, agebin_2_no/1.0e6)).T.flatten(), np.repeat(dtd_no[:,0],2),\
                     color=clr, ls='--', lw=2.0)
            if show_errors:
                plt.errorbar(ages_yes/1.0e6, dtd_yes[:,0], yerr=(dtd_yes_error_low, dtd_yes_error_high), \
                             uplims=False, fmt='', ls='', ecolor='k', label='$\sigma$ = stat')
                plt.errorbar(ages_yes/1.0e6, dtd_yes[:,0], yerr=(dtd_yes_toterr_low, dtd_yes_toterr_high), \
                             uplims=False, fmt='', ls='', ecolor='k', label='$\sigma$ = stat+SFH')
            #plt.axvspan(np.repeat(dtd,2), )

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Time [Myrs]', fontsize=15)
    plt.xlim(5.0, 2.0e4)
    plt.ylim(1.0e-1,5.0e2)
    plt.title('OGLE {0} - {1}'.format(objName, sfhBinType))
    plt.ylabel(r'$\rm{\Psi\ T_{vis}}$ $\rm{\left(Number/10^6\ M_{\odot}\right)}$', fontsize=15)
    plt.tick_params(axis='y', labelsize=14)
    plt.tick_params(axis='x', labelsize=14)
    plt.legend(numpoints=1, loc=3)
    #plt.savefig(outfilePrefix+'DTD_'+objName+'_'+sfhBinType+'.pdf', dpi=150)
    plt.show()
