import numpy as np                        ##IMPORTS                                                                                                         
import os
import glob
import aplpy
import string
import sys
from scipy import misc
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy import constants as const
from astropy import units as u
from astropy import wcs
from astropy.io import fits
from astropy.io import ascii
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy import coordinates as coord
from astropy.time import Time

def sfh_ageBins(scheme):
    """
    Returns the age bins for the SFH Binning scheme used.

    Parameters:
    ----------
    """
    #SFH bins and binning schemes from Zaritsky (2009)
    logAgeLimsArr = np.asarray([6.65,6.95,7.25,7.55,7.85,8.1,8.3,8.5,8.7,8.9,9.1,9.3,9.5,9.7,9.9,10.1,10.3])
    ageLimsArr = 10.0**logAgeLimsArr
    ageIntervalsArr = ageLimsArr[1:]-ageLimsArr[:-1]

    binningSchemeNames = ['Coarse','Massive','Medium','MediumB','MediumC','Subfine','Fine','CrazyCeph','Unbinned']

    binningSchemes = [[[0,1,2,3],[4,5,6,7,8],[9,10,11,12,13,14,15]]
                  ,[[0],[1],[2],[3,4,5,6,7,8,9,10,11,12,13,14,15]]
                  ,[[0,1,2],[3,4,5,6,7,8],[9,10],[11],[12,13],[14,15]]
                  ,[[0,1,2],[3,4,5,6,7,8],[9,10],[11,12],[13],[14,15]]
                  ,[[0,1,2],[3,4,5,6],[7,8,9,10],[11,12],[13],[14,15]]
                  ,[[0,1,2],[3,4,5],[6,7,8],[9,10],[11,12],[13],[14,15]]
                  ,[[0,1,2],[3,4,5],[6,7,8],[9,10],[11],[12],[13],[14,15]]
                  ,[[0],[1],[2,3],[4],[5,6,7],[8,9],[10,11,12],[13,14,15]]    
                  ,[[0],[1],[2],[3],[4],[5],[6],[7],[8],[9],[10],[11],[12],[13],[14],[15]]]

    #Creating a dictionary of binning schemes and their bins
    binDict = {sch:bins for sch, bins in zip(binningSchemeNames, binningSchemes) }

    #Selecting the age bins based on 'scheme'
    ageBins = [ageLimsArr[sch[0]] for sch in binDict[scheme]]
    ageBins.append(ageLimsArr[-1])
    logAgeCentroids = np.log10(ageBins)[:-1] + (np.log10(ageBins[1:])-np.log10(ageBins[:-1]))/2.0
    ages = (10**logAgeCentroids)
    
    return (ages, np.array(ageBins[:-1]), np.array(ageBins[1:]))

def hpd(trace, mass_frac) :
    """
    Returns highest probability density region given by
    a set of samples.

    Parameters
    ----------
    trace : array
        1D array of MCMC samples for a single variable
    mass_frac : float with 0 < mass_frac <= 1
        The fraction of the probability to be included in
        the HPD.  For example, `massfrac` = 0.95 gives a
        95% HPD.
        
    Returns
    -------
    output : array, shape (2,)
        The bounds of the HPD

    Courtesy: http://bebi103.caltech.edu.s3-website-us-east-1.amazonaws.com/2015/tutorials/l06_credible_regions.html
    """
    # Get sorted list
    d = np.sort(np.copy(trace))

    # Number of total samples taken
    n = len(trace)
    
    # Get number of samples that should be included in HPD
    n_samples = np.floor(mass_frac * n).astype(int)
    
    # Get width (in units of data) of all intervals with n_samples samples
    int_width = d[n_samples:] - d[:n-n_samples]
    
    # Pick out minimal interval
    min_int = np.argmin(int_width)
    
    # Return interval
    return np.array([d[min_int], d[min_int+n_samples]])

def dtd_mode(data):
    n, bins = np.histogram(data, bins=100)
    x = bins[:-1] + (bins[1:] - bins[:-1]) / 2.0
    return x[np.where(n == n.max())][0]
        
def dtd_staterrors(data_nominal):
    """
    Returns the statistical errors on the dtd bins based on highest posterior criterion

    Parameter
    ---------

    data_nominal: ndarray
                MCMC chain contain dtd for the nominal SFH

    Return
    ------
    med: ndarray, float
           Median value of the DTD in the given bin or 2-sigma upper limit

    staterror: ndarray, float
           Statistical errors (or just median - mode if upper limit)

    isuplim: ndarray, boolean
           Tells whether each bin has a detection or not
    """
    med = np.zeros(data_nominal.shape[1])
    staterror = []
    isuplim = []
    for i, data in enumerate(data_nominal.T):

        hpd_2sigma = hpd(data, 0.95) 
        hpd_1sigma = hpd(data, 0.68) 
    
        n, bins = np.histogram(data, bins=100)
        x = bins[:-1] + (bins[1:] - bins[:-1]) / 2.0
        mode = x[np.where(n == n.max())][0]   #Extra because sometimes two modes might correspond to n.max()


        if (2.0*mode - hpd_2sigma[1]) < 0:
            isuplim.append(True)
            med[i]      = hpd_2sigma[1] 
            staterror.append(med[i]-mode)
  
    
        else:
            isuplim.append(False)
            med[i]      = mode
            staterror.append(np.array([ med[i] - hpd_1sigma[0], hpd_1sigma[1] - med[i]]))

    return (med, staterror, isuplim)

def dtd_stat_randomsad_errors(data_nominal, data_SFH):
    """
    Returns the DTD per bin with combined statistical and sad errors 

    Parameter
    ---------

    data_nominal: ndarray
                MCMC chain contain dtd for the nominal SFH

    Return
    ------
    med: ndarray, float
           Median value of the DTD in the given bin or 2-sigma upper limit

    staterror: ndarray, float
           Statistical errors (or just median - mode if upper limit)

    isuplim: ndarray, boolean
           Tells whether each bin has a detection or not
    """
    med = []
    staterror = []
    isuplim = []
    for (data_nom, data_sfh) in zip(data_nominal.T, data_SFH.T):

        hpd_nom_2sigma = hpd(data_nom, 0.95) 
        hpd_nom_1sigma = hpd(data_nom, 0.68)
        mode_nom = dtd_mode(data_nom)

        hpd_sfh_2sigma = hpd(data_sfh, 0.95)
        hpd_sfh_1sigma = hpd(data_sfh, 0.68)
        mode_sfh = dtd_mode(data_sfh)

        mode = (mode_nom + mode_sfh)/2.0
        dtd_2sigma = np.sqrt((hpd_sfh_2sigma[1] - mode_sfh)**2 + (hpd_nom_2sigma[1] - mode_nom)**2)
        dtd_1sigma_low = np.sqrt((mode_sfh - hpd_sfh_1sigma[0])**2 + (mode_nom - hpd_nom_1sigma[0])**2)
        dtd_1sigma_up = np.sqrt((hpd_sfh_1sigma[1] - mode_sfh)**2 + (hpd_nom_1sigma[1] - mode_sfh)**2)

        if (mode - dtd_2sigma) < 0:
            isuplim.append(True)
            med.append(mode + dtd_2sigma) 
            staterror.append(dtd_2sigma)
  
    
        else:
            isuplim.append(False)
            med.append(mode)
            staterror.append(np.array([dtd_1sigma_low , dtd_1sigma_up]))

    return (med, staterror, isuplim)

def dtd_randomsad_errors(data_SFH):
    """
    Returns the DTD per bin with sad errors only

    Parameter
    ---------

    data_nominal: ndarray
                MCMC chain contain dtd for the nominal SFH

    Return
    ------
    med: ndarray, float
           Median value of the DTD in the given bin or 2-sigma upper limit

    staterror: ndarray, float
           Statistical errors (or just median - mode if upper limit)

    isuplim: ndarray, boolean
           Tells whether each bin has a detection or not
    """
    med = []
    error = []
    isuplim = []
    for data_sfh in data_SFH.T:


        hpd_sfh_2sigma = hpd(data_sfh, 0.95)
        hpd_sfh_1sigma = hpd(data_sfh, 0.68)
        mode_sfh = dtd_mode(data_sfh)

        dtd_2sigma = hpd_sfh_2sigma[1] - mode_sfh
        dtd_1sigma_low = mode_sfh - hpd_sfh_1sigma[0]
        dtd_1sigma_up = hpd_sfh_1sigma[1] - mode_sfh

        if (mode_sfh - dtd_2sigma) < 0:
            isuplim.append(True)
            med.append(mode_sfh + dtd_2sigma) 
            error.append(dtd_2sigma)
  
    
        else:
            isuplim.append(False)
            med.append(mode_sfh)
            error.append(np.array([dtd_1sigma_low , dtd_1sigma_up]))

    return (med, error, isuplim)


def dtd_saderrors(data_low, data_nominal, data_high):
    """
    Returns the errors due to star-formation history on the dtd bins

    Parameter
    ---------
    data_nominal: ndarray
    data_low    : ndarray
    data_high   : ndarray
    Return
    ------
    The +,- SAD errors
    
    """

    med_nominal = np.array([dtd_mode(data) for data in data_nominal.T])
    med_low     = np.array([dtd_mode(data) for data in data_low.T])
    med_high    = np.array([dtd_mode(data) for data in data_high.T])

    return (med_nominal - med_low, med_high - med_nominal)

def dtd_stat_sad_errors(data_nominal, data_low, data_high):
    """
    Returns the statistical errors on the dtd bins based on highest posterior criterion

    Parameter
    ---------

    data_nominal: ndarray
                MCMC chain contain dtd for the nominal SFH

    Return
    ------
    med: ndarray, float
           Median value of the DTD in the given bin or 2-sigma upper limit

    staterror: ndarray, float
           Statistical errors (or just median - mode if upper limit)

    isuplim: ndarray, boolean
           Tells whether each bin has a detection or not
    """
    med = np.zeros(data_nominal.shape[1])
    staterror = []
    isuplim = []
    sad_err_low, sad_err_up = dtd_saderrors(data_low, data_nominal, data_high)

    for i, (data_n, data_l, data_h) in enumerate(zip(data_nominal.T, data_low.T, data_high.T)):

        mode = dtd_mode(data_n)


        if (mode - 2.0*sad_err_up[i] < 0):
            isuplim.append(True)
            med[i]      = mode + 2.0*sad_err_up[i]
            staterror.append(2.0*sad_err_up)
  
    
        else:
            isuplim.append(False)
            med[i]      = mode
            staterror.append(np.array([ sad_err_low[i], sad_err_up[i]]))

    return (med, staterror, isuplim)

def is_notdetect(dtd, dtd_sigma):
    """
    Returns an array of truths if signal is below 2 sigma in a delay-time bin.

    """   
    return [True if dtd[i]-2.0*dtd_sigma[i] > 0. else False for i in range(dtd.size)]

def dtd_with_errors(data_low, data_nominal, data_high, stat_errortype='percentileA', only_staterror=False):
    """
    Returns the DTD from the chains with error bars.
    
    Arguments
    ---------
    stat_errortype: string
             Pertains to calculation of statistical errors. The types are :-
             
            'percentileA' - asymmetric 1sigma regions, difference between
                            84th and median, and median and 16th percentile.
            'percentileB' - symmetric 1sigma regions, difference between 
                            84th and 16th percentile, divided by 2.
            'hpd'         - highest posterior density method.
            
    only_staterror: boolean
            If true, returns only the statistical errors. Else the statistical
            and systematic errors (e.g. from SAD) added in quadrature.
            
    Returns
    -------
    
    psi_Tpn_median_nominal: float
        The median dtd value
        
    ..and the errors.
        
    """
    psi_Tpn_median_low = np.array([np.median(data_low[:,binnum]) for binnum in range(6)])
    psi_Tpn_median_high = np.array([np.median(data_high[:,binnum]) for binnum in range(6)])

    psi_Tpn_median_nominal = np.array([np.median(data_nominal[:,binnum]) for binnum in range(6)])
    
    if stat_errortype == 'percentileA':
        psi_Tpn_1sigma_nominal_up = np.array([(np.percentile(data_nominal[:,binnum], 84) - psi_Tpn_median_nominal[binnum]) for binnum \
                                      in range(6)])
        psi_Tpn_1sigma_nominal_down = np.array([(psi_Tpn_median_nominal[binnum] - np.percentile(data_nominal[:,binnum], 16)) for binnum \
                                      in range(6)])
        
    if stat_errortype == 'percentileB':
        psi_Tpn_1sigma_nominal_up = (np.array([(np.percentile(data_nominal[:,binnum], 84) - np.percentile(data_nominal[:,binnum], 16)) for binnum \
                                      in range(6)]))*0.5
        psi_Tpn_1sigma_nominal_down = psi_Tpn_1sigma_nominal_up
        
    if only_staterror == True:
        return (psi_Tpn_median_nominal, psi_Tpn_1sigma_nominal_down, psi_Tpn_1sigma_nominal_up)
    
    else:
        psi_Tpn_1sigma_sad_up = psi_Tpn_median_high - psi_Tpn_median_nominal
        psi_Tpn_1sigma_sad_down = psi_Tpn_median_nominal - psi_Tpn_median_low
        psi_Tpn_error_up = np.sqrt(psi_Tpn_1sigma_nominal_up**2 + psi_Tpn_1sigma_sad_up**2)
        psi_Tpn_error_down = np.sqrt(psi_Tpn_1sigma_nominal_down**2 + psi_Tpn_1sigma_sad_down**2)
        return (psi_Tpn_median_nominal, psi_Tpn_error_down, psi_Tpn_error_up)

def object_ra_decs(fileName, objClassName, obj_subtype='All'):

    print objClassName+' in the LMC.'
    print 'Reading object catalog from file '+fileName

    #RRLyrae from OGLE is in ASCII format, so it can easily read with this codeblock instead of running for loops..
    data = ascii.read(fileName)
    data = data if obj_subtype=='All' else data[data['Type'] == obj_subtype]
    radec_string = [r+' '+d for r,d in zip(data['RA'], data['Dec'])]
    c = SkyCoord(radec_string, unit='deg')
    return (c.ra.deg, c.dec.deg)        


    
def objectRADecs(fileName, firstLine, lastLine, objClassName):
    """
    (Under construction)
    Returns the RA and Declinations from a catalog file. 

    Parameters
    ----------------------
    fileName: string
            Name of file.

    firstLine: int
            Line number where object list starts.

    lastLine: int
            Line number where object list ends.
            
    objClassName: string
            Type of object


    Returns
    ----------------------
    objListRA: list
            The list of RAs for the objects

    objListDec: list
            The list of Decs for the objects
    """

    print objClassName+' in the LMC.'
    print 'Reading object catalog from file '+fileName
    f = open(fileName,'r')
    lineList = f.readlines()
    f.close()
    objListRA = []
    objListDec = []

    #RRLyrae from OGLE is in ASCII format, so it can easily read with this codeblock instead of running for loops..
    if objClassName == 'RRLyrae':
        data = ascii.read(fileName)
        radec_string = [r+' '+d for r,d in zip(data['RA'], data['Dec'])]
        c = SkyCoord(radec_string, unit='deg')
        return (c.ra.deg, c.dec.deg)        

    for line in range(firstLine, lastLine) :
        words = lineList[line].split()
        
        #    objRA =  float(words[3]) + float(words[4])/60.0 + float(words[5])/3600.0                                                                         
    #if float(words[6]) < 0 :                                                                                                                                 
    #    objDec =  float(words[6]) - float(words[7])/60.0 - float(words[8])/3600.0                                                                            
    #else :                                                                                                                                                   
    #    objDec =  float(words[6]) + float(words[7])/60.0 + float(words[8])/3600.0                                                                

        #print line, objRA, objDec                                                                                                                          

        if objClassName == 'RCrBor' :
            raHMSWords = words[3].split(':')
            objRA = float(raHMSWords[0]) + float(raHMSWords[1])/60.0 + float(raHMSWords[2])/3600.0
            decDMSWords = words[4].split(':')
            objDec = float(decDMSWords[0]) - float(decDMSWords[1])/60.0 - float(decDMSWords[2])/3600.0
        elif objClassName == 'PNe' :
            raHMSWords = words[1].split(':')
            objRA = float(raHMSWords[0]) + float(raHMSWords[1])/60.0 + float(raHMSWords[2])/3600.0
            decDMSWords = words[2].split(':')
            objDec = float(decDMSWords[0]) - float(decDMSWords[1])/60.0 - float(decDMSWords[2])/3600.0
        elif objClassName == 'PNe0' :
            raHMSWords = words[1].split(':')
            objRA = float(raHMSWords[0]) + float(raHMSWords[1])/60.0 + float(raHMSWords[2])/3600.0
            decDMSWords = words[2].split(':')
            objDec = float(decDMSWords[0]) - float(decDMSWords[1])/60.0 - float(decDMSWords[2])/3600.0
        elif objClassName == 'Novae' :
            objRA = float(words[2]) + float(words[3])/60.0 + float(words[4])/3600.0
            objDec = float(words[5]) - float(words[6])/60.0 - float(words[7])/3600.0
        elif objClassName == 'SNRs' :
            objRA = float(words[3]) + float(words[4])/60.0 + float(words[5])/3600.0
            objDec = float(words[6]) - float(words[7])/60.0 - float(words[8])/3600.0
        elif objClassName[0:2] == 'WR' :
            objRA = float(words[0])
            objDec = float(words[1])
        elif objClassName[0:2] == 'Os' :
            objRA = float(words[0])
            objDec = float(words[1])
        elif objClassName == 'RSG' :
            objRA = float(words[0])
            objDec = float(words[1])
        else :
            objRA = float(words[4])
            objDec = float(words[5])

        if objClassName == 'PNe' :
            if words[6] == 'T' :
                objListRA.append(objRA)
                objListDec.append(objDec)
        else :
            objListRA.append(objRA)
            objListDec.append(objDec)


    print 'Objects found in file: ', len(objListRA)
    return (objListRA, objListDec)

def cellsinSFHSurvey(cellNames, cellCentersRA, cellCentersDec):
    """
    (Under Construction)
    Returns the locations of SFH cells that fall inside the survey area of the object

    Parameters
    ----------

    cellNames: 

    cellCentersRA: 

    cellCentersDec:

    Returns
    -------

    cellNames_inSurvey:

    """

    fieldCenter = SkyCoord('05h22m00s', '-69d00m00s', frame=coord.FK5(equinox=Time('B1950', scale='utc')))
    fieldCenter.transform_to(coord.FK5(equinox=Time(2000, format='jyear', scale='utc')))
    cellNames = np.asarray(cellNames)
    surveyExtent = 5.25/2.0 #degrees                                                                                                                         

    print 'Field Center (RA, Dec) : ', fieldCenter.ra.deg, fieldCenter.dec.deg

    is_RAinSurvey = (cellCentersRA >= (fieldCenter.ra.deg - surveyExtent)) & (cellCentersRA <= (fieldCenter.ra.deg + surveyExtent))
    cellCentersRA = cellCentersRA[is_RAinSurvey]
    cellCentersDec = cellCentersDec[is_RAinSurvey]
    cellNames = cellNames[is_RAinSurvey]
    is_DecinSurvey = (cellCentersDec >= (fieldCenter.dec.deg - surveyExtent)) & (cellCentersDec <= (fieldCenter.dec.deg + surveyExtent))
    cellCentersRA_inSurvey = cellCentersRA[is_DecinSurvey]
    cellCentersDec_inSurvey = cellCentersDec[is_DecinSurvey]
    cellNames_inSurvey = cellNames[is_DecinSurvey]
    return cellNames_inSurvey
