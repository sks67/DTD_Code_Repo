import numpy as np
import matplotlib.pyplot as plt
import corner
from astropy.io import fits
import dtdutils
reload(dtdutils)
from scipy import interpolate
from scipy import stats
from astropy.io import ascii
from scipy import integrate
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table, Column

DTDpath = '/Users/sumits2k/Desktop/Research/SNResearch2/RadioSNRs/DTD/'

def errorTable(object_name, age_1, age_2, dtd, dtd_error, isdetect, det_sig):
    """
    Create a Table of best fit DTDs, the errors and detection significances. Can write both ASCII
    and AASTEX formatted tables. 

    """
    
    dtdTable = Table()
    dtdTable['Age (left edge)'] = age_1
    dtdTable['Age (right edge)'] = age_2
    dtdTable['DTD (best fit)'] = np.around(dtd, decimals=7)
    dtdTable['DTD (error)'] = np.around(dtd_error, decimals=7)
    dtdTable['Detection?'] = ~isdetect
    dtdTable['Det. Significance'] = np.around(det_sig, decimals = 2)
    dtdTable.write(object_name + '_' + 'DTD_Error_Table.txt', overwrite=True, format='ascii')

def averagePDF(param, bins):
    p_x_array = np.zeros((nIterations, bins.size-1))
    for i in range(nIterations):
        p_x_array[i], binedges = np.histogram(masterchain_randSFH[i,:,:,param].flatten()*1.0e6, density=False, bins=bins)

    return np.mean(p_x_array.T, axis=1)


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
              colorScheme='bmh', show_errors = False, show_uplims=False, show_title=True):
    
    color_options = {'bmh':[u'#348ABD', u'#A60628',u'#467821',  u'#D55E00',u'#CC79A7',\
                            u'#56B4E9', u'#009E73', u'#F0E442', u'#0072B2'],
                    'seaborn-colorblind':[u'#0072B2', u'#009E73', u'#D55E00', u'#CC79A7', u'#F0E442',\
                                          u'#56B4E9'],
                     'rgb':['r', 'g', 'b', 'k', 'm'],
                    'seq-blue':reversed(['#fff7fb','#ece7f2','#d0d1e6','#a6bddb','#74a9cf','#3690c0','#0570b0',\
                                 '#045a8d','#023858'])}
    colors = iter(color_options[colorScheme])
    labels = iter(obj_subtype)
    
    plt.rc('font', family='serif')
    plt.figure(figsize=(6,5))
    for data_nominal, data_low, data_high in zip(data_nom_arr, data_low_arr, data_high_arr):
        
        ages, agebin_1, agebin_2 = dtdutils.sfh_ageBins(sfhBinType)

        #Extracting the signals (or lack thereof) with statistical errors
        dtd, dtd_error, isuplim = dtdutils.dtd_staterrors(data_nominal)
        
        #dtd_SFH, dtd_error_SFH, isuplim_SFH = dtdutils.dtd_staterrors(data_randSFH*1.0e6)
        
        
        #Extracting the systematic errors due to different SFHs
        dtd_sad_error_low, dtd_sad_error_high = dtdutils.dtd_saderrors(data_low, data_nominal, data_high)

        dtd_no = np.array([signal if isuplim[i] else 0. for i, signal in enumerate(dtd)])
        dtd_yes = np.array([signal if np.invert(isuplim[i]) else 0. for i,signal in enumerate(dtd)])
        
        if show_errors:
            dtd_y = np.array([signal for i, signal in enumerate(dtd) if np.invert(isuplim[i])])
            ages_y = np.array([age for i, age in enumerate(ages) if np.invert(isuplim[i])])
            
            dtd_yes_error = np.array([err for i, err in enumerate(dtd_error) if np.invert(isuplim[i])])          
            dtd_yes_error_low = dtd_yes_error[:,0]
            dtd_yes_error_high = dtd_yes_error[:,1]
            
            dtd_yes_sad_error_low = np.array([err for i, err in enumerate(dtd_sad_error_low) if np.invert(isuplim[i])])
            dtd_yes_sad_error_high = np.array([err for i, err in enumerate(dtd_sad_error_high) if np.invert(isuplim[i])])
   
            dtd_yes_toterr_low = np.sqrt(dtd_yes_error_low**2 + dtd_yes_sad_error_low**2)
            dtd_yes_toterr_high = np.sqrt(dtd_yes_error_high**2 + dtd_yes_sad_error_high**2)

            
        #Separating bins with upper limits
        if False:
            agebin_1_no = agebin_1[dtd_no[:,1].astype(int)]
            agebin_2_no = agebin_2[dtd_no[:,1].astype(int)]
    
        
        #Separating bins with signals
        if False:
            agebin_1_yes = agebin_1[dtd_yes[:,1].astype(int)]
            agebin_2_yes = agebin_2[dtd_yes[:,1].astype(int)]
            dtd_yes_sad_error_low = dtd_sad_error_low[dtd_yes[:,1].astype(int)]
            dtd_yes_sad_error_high = dtd_sad_error_high[dtd_yes[:,1].astype(int)]
       
            
        #Total error

        with plt.style.context(('bmh')):
            #This temporarily sets the plotting style. Do not use plt.style.use()
            clr = colors.next()
            if False:
                plt.plot(np.vstack((agebin_1_yes/1.0e6, agebin_2_yes/1.0e6)).T.flatten(), np.repeat(dtd_yes[:,0],2),\
                         color=clr, ls='-',lw=2.0, label=labels.next())
                plt.plot(np.vstack((agebin_1_no/1.0e6, agebin_2_no/1.0e6)).T.flatten(), np.repeat(dtd_no[:,0],2),\
                         color=clr, ls='--', lw=2.0)
                plt.errorbar(ages_y/1.0e6, dtd_y[:,0], yerr=(dtd_yes_error_low, dtd_yes_error_high), \
                             uplims=False, fmt='', ls='', ecolor='k', label='$\sigma$ = stat')
                plt.errorbar(ages_yes/1.0e6, dtd_yes[:,0], yerr=(dtd_yes_toterr_low, dtd_yes_toterr_high), \
                             uplims=False, fmt='', ls='', ecolor='k', label='$\sigma$ = stat+SFH')
            #plt.axvspan(np.repeat(dtd,2), )
            
            x_axis = np.vstack((agebin_1/1.0e6, agebin_2/1.0e6)).T.flatten()
            y_axis = np.repeat(dtd_yes,2)
            plt.plot(x_axis, y_axis, color = clr, ls = '-', lw=2.0)
            plt.gca().fill_between(x_axis, 0., y_axis, facecolor=clr, alpha=0.8, label=labels.next())
            if show_errors:
                plt.errorbar(ages_y/1.0e6, dtd_y, yerr=(dtd_yes_toterr_low, dtd_yes_toterr_high), \
                         uplims=False, fmt='', ls='', ecolor='k')#, label='$\sigma$ = stat+SFH')
            if show_uplims:
                plt.plot(x_axis, np.repeat(dtd_no,2), color=clr, ls='--', lw=2.0)


    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Time [Myrs]', fontsize=15)
    plt.xlim(5.0, 2.0e4)
    plt.ylim(1.0e-9,1.0e-3)
    if show_title:
        plt.title('OGLE {0} - {1}'.format(objName, sfhBinType))
    plt.ylabel(r'$\rm{\Psi\ T_{vis}}$ $\rm{\left(Number/10^{6}\ M_{\odot}\right)}$', fontsize=15)
    plt.tick_params(axis='y', labelsize=14)
    plt.tick_params(axis='x', labelsize=14)
    plt.legend(numpoints=1)
#    plt.savefig(outfilePrefix+'DTD_'+objName+'_'+sfhBinType+'.pdf', dpi=150)
    plt.show()

def color_options(key):
    colors = {'bmh':[u'#348ABD', u'#A60628',u'#467821',  u'#D55E00',u'#CC79A7',\
                            u'#56B4E9', u'#009E73', u'#F0E442', u'#0072B2'],
                    'seaborn-colorblind':[u'#0072B2', u'#009E73', u'#D55E00', u'#CC79A7', u'#F0E442',\
                                          u'#56B4E9'],
                     'rgb':['r', 'g', 'b', 'k', 'm'],
                    'seq-blue':reversed(['#fff7fb','#ece7f2','#d0d1e6','#a6bddb','#74a9cf','#3690c0','#0570b0',\
                                 '#045a8d','#023858'])}

    return colors[key]

def plot_dtds_SFH(data_nom_arr, data_randSFH_arr, trueDTD, objName, obj_subtype, sfhBinType, \
              colorScheme='bmh', mult_factor=1., stat_sad_crit = True, show_errors = False, show_uplims=False, plot_title='',\
                  ylabel='', savefigure=True, savefile=''):

    colors = iter(color_options(colorScheme))
    labels = iter(obj_subtype)
    sfhBinType = sfhBinType*len(obj_subtype) if len(sfhBinType)==1 else sfhBinType
    bintype = iter(sfhBinType)
    
    plt.figure(figsize=(7,6))
    
    #Upper and lower 1sigma errors from combined randomized SFH and statistical errors.
    #dtd_stat, dtd_sfh, dtd_1sig_low, dtd_1sig_up, dtd_2sig = dtd_stat_randomSFH_errors(data_randSFH, data_nom_arr)
    
    
    for data_nominal, data_randSFH in zip(data_nom_arr, data_randSFH_arr):
        
        ages, agebin_1, agebin_2 = dtdutils.sfh_ageBins(bintype.next())
        print ages

        
        #DTD errors
        if stat_sad_crit:
            #Extracting dtd per bin based on combined statistical and randomized SFH errors
#            dtd, dtd_error, actual_dtd, actual_dtd_err, isuplim, det_sig = dtdutils.dtd_stat_randomsad_errors(data_nominal, data_randSFH, mult_factor=mult_factor)
            dtd, dtd_error, isuplim, det_sig = dtdutils.dtd_staterrors(data_nominal)            
        else:
            #Extracting dtd per bin based only on systematic errors
            dtd, dtd_error, actual_dtd, actual_dtd_err, isuplim, det_sig = dtdutils.dtd_randomsad_errors(data_randSFH, mult_factor=mult_factor)

        print '\n\n DTDs : ', dtd
        print '\n\n DTD Errors: ', dtd_error
        print '\n\nDetections? : ', isuplim

#        errorTable(objName, agebin_1, agebin_2, dtd, dtd_error, isuplim, det_sig)
        
        dtd_no = np.array([signal if isuplim[i] else 0. for i, signal in enumerate(dtd)])
        dtd_yes = np.array([signal  if np.invert(isuplim[i]) else 0. for i,signal in enumerate(dtd)])
        
        if show_errors:
            dtd_y = np.array([signal for i, signal in enumerate(dtd) if np.invert(isuplim[i])])
            ages_y = np.array([age for i, age in enumerate(ages) if np.invert(isuplim[i])])
            
            dtd_yes_error = np.array([err for i, err in enumerate(dtd_error) if np.invert(isuplim[i])])  
            print dtd_yes_error.shape
            dtd_yes_error_low = dtd_yes_error[:,0]
            dtd_yes_error_high = dtd_yes_error[:,1]

    

#        with plt.style.context(('bmh')):
            #This temporarily sets the plotting style. Do not use plt.style.use()
        clr = colors.next()
        
        x_axis = np.vstack((agebin_1/1.0e6, agebin_2/1.0e6)).T.flatten()
        y_axis = np.repeat(dtd_yes,2)
        
        plt.plot(x_axis, y_axis, color = clr, ls = '-', lw=2.0, label='Recovered DTD')
        plt.gca().fill_between(x_axis, 0., y_axis, facecolor=clr, alpha=0.8)#, label=labels.next())
        
        if show_errors:
            #                plt.errorbar(ages/1.0e6, actual_dtd, yerr = (actual_dtd_err[:, 0], actual_dtd_err[:, 1]), ls='', fmt='', ecolor='k')
            plt.errorbar(ages_y/1.0e6, dtd_y, yerr=(dtd_yes_error_low, dtd_yes_error_high), \
                         uplims=False, fmt='', ls='', ecolor='k', label='')#label='$\sigma$ = stat+randomSFH')
#                plt.errorbar(ages_y/1.0e6, dtd_y, yerr=(dtd_randsfh_err_low[10:], dtd_randsfh_err_up[10:]), \
 #                        uplims=False, fmt='', ls='', ecolor='y', label='$\sigma$ = stat+SFH (Randomized)')
        if show_uplims:
           # plt.plot(x_axis, np.repeat(dtd_no,2), color=clr, ls='--', lw=2.0)
            #switch to color=clr when needed
            plt.errorbar(ages/1.0e6, dtd_no, yerr=dtd_no*0.25, xerr=((ages-agebin_1)/1.0e6, (agebin_2 - ages)/1.0e6),  fmt='', ls='', ecolor='k', color='k', uplims=True)

        if trueDTD.size>0:
            y_trueDTD = np.repeat(trueDTD, 2)
            plt.plot(x_axis, y_trueDTD, color='k', ls='--', lw=2.0, label='True DTD')

    [plt.text(agebin_1[i]/1.0e6, dtd[i]*(10**0.2), str(np.around(det_sig[i], decimals=1))+r'$\sigma$', fontsize=8) for i in range(ages.size)]
    plt.title(plot_title, fontsize=15)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Delay Time [Myrs]', fontsize=15)
    plt.xlim(4.0, 2.0e4)
#    plt.ylim(6.0e-7,2.0e-4)
    plt.ylim(9.0e-8, 2.0e-4)
#    plt.ylabel(r'$\mathrm{\Psi T_{RR}}$ $\rm{\left(Number/M_{\odot}\right)}$', fontsize=15)
    plt.ylabel(ylabel+r' $\rm{\left(Number/M_{\odot}\right)}$', fontsize=15)
    plt.tick_params(axis='y', labelsize=14)
    plt.tick_params(axis='x', labelsize=14)
    plt.legend(loc=1, numpoints=1)
    if savefigure:
        plt.savefig(savefile, dpi=100)
    plt.show()

def plot_dtds_AllErrors(data_nom_arr, data_randSFH_arr, trueDTD, objName, obj_subtype, sfhBinType, \
              colorScheme='bmh', mult_factor=1., stat_sad_crit = True, show_errors = False, show_uplims=False, show_title=False,\
                  savefigure=True, savefile=''):
    """
    Same as plot_dtds_SFH, but this plots all the DTDs as histograms, with modes and errors bars
    from statistical and SAD uncertainties. This is to check the ratio of the errors in the different age bins.
    
    CREATED: 11/24/2017
    """
 

    colors = iter(color_options(colorScheme))
    labels = iter(obj_subtype)
    sfhBinType = sfhBinType*len(obj_subtype) if len(sfhBinType)==1 else sfhBinType
    bintype = iter(sfhBinType)
    
    plt.figure(figsize=(7,6))
    
    #Upper and lower 1sigma errors from combined randomized SFH and statistical errors.
    #dtd_stat, dtd_sfh, dtd_1sig_low, dtd_1sig_up, dtd_2sig = dtd_stat_randomSFH_errors(data_randSFH, data_nom_arr)
    
    
    for data_nominal, data_randSFH in zip(data_nom_arr, data_randSFH_arr):
        
        ages, agebin_1, agebin_2 = dtdutils.sfh_ageBins(bintype.next())
        dtd_nom, dtd_sfh, dtd_stat_up, dtd_stat_down, dtd_sfh_up, dtd_sfh_down, det_sig_stat, det_sig_sfh = dtdutils.dtd_all_errors(data_nominal, data_randSFH)
        dtd = (dtd_nom + dtd_sfh)/2.
        dtd_n, dtd_error, actual_dtd, actual_dtd_err, isuplim, det_sig = dtdutils.dtd_randomsad_errors(data_randSFH, mult_factor=mult_factor) 
        dtd_no = np.array([signal if isuplim[i] else 0. for i, signal in enumerate(dtd)])
        dtd_yes = np.array([signal  if np.invert(isuplim[i]) else 0. for i,signal in enumerate(dtd)])
        
        if show_errors:
            dtd_y = np.array([signal for i, signal in enumerate(dtd) if np.invert(isuplim[i])])
            ages_y = np.array([age for i, age in enumerate(ages) if np.invert(isuplim[i])])
            
            dtd_yes_error = np.array([err for i, err in enumerate(dtd_error) if np.invert(isuplim[i])])  
            print dtd_yes_error.shape
            dtd_yes_error_low = dtd_yes_error[:,0]
            dtd_yes_error_high = dtd_yes_error[:,1]

    
        print 'Detection Sigma (Stat): ', det_sig_stat
        print 'Detection Sigma (SFH): ', det_sig_sfh

#        with plt.style.context(('bmh')):
            #This temporarily sets the plotting style. Do not use plt.style.use()
        clr = colors.next()
        
        x_axis = np.vstack((agebin_1/1.0e6, agebin_2/1.0e6)).T.flatten()
        y_axis = np.repeat(dtd_yes,2)

        plt.plot(x_axis, y_axis, color = clr, ls = '-', lw=2.0, label='Recovered DTD')
        plt.gca().fill_between(x_axis, 0., y_axis, facecolor=clr, alpha=0.5)#, label=labels.next())

#        plt.errorbar(ages_y/1.0e6, dtd_y, yerr=(dtd_yes_error_low, dtd_yes_error_high), \
#                         uplims=False, fmt='', ls='', ecolor='k', label='')#label='$\sigma$ = stat+randomSFH')        

        plt.errorbar(ages/1.0e6, dtd_nom, yerr=(dtd_stat_down, dtd_stat_up),  xerr=((ages-agebin_1)/1.0e6, (agebin_2 - ages)/1.0e6), fmt='o', ls='', color='#d95f02', ecolor='#d95f02', label='Stat error')
        plt.errorbar(ages/1.0e6, dtd_sfh, yerr=(dtd_sfh_down, dtd_sfh_up),  xerr=((ages-agebin_1)/1.0e6, (agebin_2 - ages)/1.0e6), fmt='o', ls='', color='#66a61e', ecolor='#66a61e', lw=2.0, alpha=0.5, label='SAD error')
        if trueDTD.size>0:
            y_trueDTD = np.repeat(trueDTD, 2)
            plt.plot(x_axis, y_trueDTD, color='k', ls='--', lw=2.0, label='True DTD')
        

    [plt.text((agebin_1[i]*10**(0.05))/1.0e6, dtd[i]*(10**(0.25)), str(np.around(det_sig_sfh[i], decimals=1))+r'$\sigma$', fontsize=8) for i in range(ages.size)] 
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Delay Time [Myrs]', fontsize=15)
    plt.xlim(5.0, 2.0e4)
#    plt.ylim(6.0e-7,2.0e-4)
    plt.ylim(1.0e-9, 1.0e-3)
    if show_title:
        plt.title('OGLE {0}'.format(objName))
#    plt.ylabel(r'$\mathrm{\Psi T_{RR}}$ $\rm{\left(Number/M_{\odot}\right)}$', fontsize=15)
    plt.ylabel(r'RR Lyrae DTD, $\rm{\left(Number/M_{\odot}\right)}$', fontsize=15)
    plt.tick_params(axis='y', labelsize=14)
    plt.tick_params(axis='x', labelsize=14)
    plt.legend(loc=1, numpoints=1)
    if savefigure:
        plt.savefig(savefile, dpi=100)
    plt.show()


def show_MESA():
    plt.figure(figsize=(6.5,8))    
    prefix = DTDpath + 'MESA_sims/'
    mass = np.arange(0.8, 2.1, 0.1)
    age_CHeB_z00 = np.loadtxt(prefix + 'CoreHeBAges_solarZ.txt')
    age_CHeB_z01 = np.loadtxt(prefix + 'CoreHeBAges_025solarZ.txt')
    age_CHeB_z001 = np.loadtxt(prefix + 'CoreHeBAges_01solarZ.txt')
    
    #Plot luminosity vs age
    plt.plot(age_CHeB_z00/1.0e6, mass, 'o', color='#bd0026', label=r'$Z=Z_{\odot}$')
    plt.plot(age_CHeB_z01/1.0e6, mass, 'o', color='#fd8d3c', label=r'$Z=0.25 Z_{\odot}$')
    plt.plot(age_CHeB_z001/1.0e6, mass, 'o', color='#fed976',label=r'$Z=0.1 Z_{\odot}$')
    plt.hlines(y=mass, xmin=age_CHeB_z001/1.0e6, xmax=age_CHeB_z00/1.0e6, linestyle='solid', color='k')
    plt.text(7,1.4, 'MIST models\n(single, non-rotating,\nsolar-scaled stars)', fontsize=11)
    plt.title('Age when Core He-burn starts [Myrs]', fontsize=14)
    plt.xscale('log')
    #Plot marker for the diferent evolutionary phases
    plt.xticks(fontsize=0)
    plt.tick_params(axis='y', labelsize=14)
    plt.ylabel(r'ZAMS mass [$M_{\odot}$]', fontsize=15)
    plt.xlim(5.0, 2.0e4)
    plt.ylim(0.7,2.2)
    plt.legend(numpoints=1, loc=2)
    plt.show()
    
def dtd_stat_randomSFH_errors(chain_SFH, chain_nom):

    dtd_errlow_sfh = np.zeros(chain_SFH.shape[1])
    dtd_errup_sfh = np.zeros(chain_SFH.shape[1])
    dtd_95_sfh = np.zeros_like(dtd_errup_sfh)
    dtd_errlow_stat = np.zeros_like(dtd_errup_sfh)
    dtd_errup_stat = np.zeros_like(dtd_errup_sfh)
    dtd_95_stat = np.zeros_like(dtd_errup_sfh)
    mode_sfh = np.zeros_like(dtd_errup_sfh)
    mode_stat = np.zeros_like(dtd_errup_sfh)

    #Error estimates fromthe SFH chain
    for i, data in enumerate(chain_SFH.T):
        mode_sfh[i] = dtdutils.dtd_mode(data)
        dtd_errup_sfh[i] = dtdutils.hpd(data, 0.68)[1] - mode_sfh[i]
        dtd_errlow_sfh[i] = mode_sfh[i] - dtdutils.hpd(data, 0.68)[0]
        dtd_95_sfh[i] = dtdutils.hpd(data, 0.95)[1] - mode_sfh[i]

    #Error estimate from the nominal chain
    chain_stat = np.asarray(chain_nom[0])
    for i, data in enumerate(chain_stat.T):
        mode_stat[i] = dtdutils.dtd_mode(data)
        dtd_errup_stat[i] = dtdutils.hpd(data, 0.68)[1] - mode_stat[i] 
        dtd_errlow_stat[i] = mode_stat[i] - dtdutils.hpd(data, 0.68)[0]
        dtd_95_stat[i] = dtdutils.hpd(data, 0.95)[1] - mode_stat[i]
        
    dtd_2sigma = np.sqrt(dtd_95_sfh**2 + dtd_95_stat**2)
    dtd_1sigma_up = np.sqrt(dtd_errup_sfh**2 + dtd_errup_stat**2)
    dtd_1sigma_low = np.sqrt(dtd_errlow_sfh**2 + dtd_errlow_stat**2)
    
    return mode_stat, mode_sfh, dtd_1sigma_low, dtd_1sigma_up, dtd_2sigma


def check_randSFH(logM, logM_up, logM_low, randSFH, binning = 'Unbinned', fileSuffix = 'trunc', showPlot = False):
    """                                                                                                                                                       
    Check if random maps are consistent with the median and errors                                                                                            
    """
    #Pick 5 random cells whose errors we will compare                                                                                                      
    cells = np.random.randint(low = 0, high = logM.shape[-1], size = 10)
    ages, age_left, age_right = dtdutils.sfh_ageBins( binning )
    age_left_err = np.log10(ages) - np.log10(age_left)
    age_right_err = np.log10(age_right) - np.log10(ages)
    with PdfPages(DTDpath+'DTD_Plots/check_randSFH_'+fileSuffix+'.pdf') as pdf:
        for cell in cells:
            logM_cell = logM[:, cell]
            logM_up_cell = logM_up[:, cell] - logM_cell
            logM_low_cell = logM_cell - logM_low[:, cell]
            meanM_indep = np.mean(randSFH[:, cell, :], axis=1)
            stdM_indep = np.std(randSFH[:, cell, :], axis=1)

            plt.figure(figsize=(7,6))
            plt.title('Error Comparison for cell '+str(cell), fontsize=18)
            plt.errorbar(np.log10(ages), logM[:, cell], xerr=(age_left_err, age_right_err),  yerr = (logM_low_cell, logM_up_cell), ls='', fmt = 'o', ecolor =\
 'k', color='k', label='HZ09 error')
            plt.errorbar(np.log10(ages), meanM_indep, ls='', yerr = stdM_indep, fmt = 'o', ecolor = 'r', color='r', elinewidth=3.0, alpha=0.3, label='SFH map\
 error')
            plt.xlabel('Log Age [yr]', fontsize=15)
            plt.ylabel(r'Log M$_{*}$ [M$_{\odot}$]', fontsize=15)
            plt.legend(loc=2)
#            plt.ylim(-0.2, 10)
            plt.tight_layout()
            pdf.savefig()

