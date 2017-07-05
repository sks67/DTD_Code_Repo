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

#Read Generic catalogs

## fileName = '/users/badenes/observations/LMC/OGLE_LMC_Cepheids.txt'
## firstLine = 7
## lastLine = 3382
## objClassName = 'Cepheids'
## refName = 'OGLE'

## fileName = '/users/badenes/observations/LMC/OGLE_LMC_Cepheids_f.txt'
## firstLine = 7
## lastLine = 1856
## objClassName = 'Cepheids_Fund'
## refName = 'OGLE'

#fileName = '/users/badenes/observations/LMC/OGLE_LMC_Cepheids_1O.txt'
#firstLine = 7
#lastLine = 1245
#objClassName = 'Cepheids_1O'
#refName = 'OGLE'

## fileName = '/users/badenes/observations/LMC/OGLE_LMC_DeltaScuti.txt'
## firstLine = 7
## lastLine = 2781
## objClassName = 'DeltaScuti'
## refName = 'OGLE'

## fileName = '/users/badenes/observations/LMC/OGLE_LMC_RRLyrae.txt'
## firstLine = 7
## lastLine = 24913
## objClassName = 'RRLyrae'
## refName = 'OGLE'

## fileName = '/users/badenes/observations/LMC/OGLE_LMC_RCrBor.txt'
## firstLine = 7
## lastLine = 30
## objClassName = 'RCrBor'
## refName = 'OGLE'

## fileName = '/users/badenes/observations/LMC/LMC_table1.dat'
## firstLine = 0
## lastLine = 44
## objClassName = 'Novae'
## refName = 'Pietsch_etal'

## fileName = '/users/badenes/observations/LMC/lmcsnrs_Allcatalogs.txt'
## firstLine = 4
## lastLine = 57
## objClassName = 'SNRs'
## refName = 'Badenes_etal'

fileName = 'lmc_dtd.dat'
firstLine = 0
lastLine = 573
objClassName = 'PNe'
refName = 'Reid&Parker'

## fileName = '/users/badenes/observations/LMC/lmcRSG.txt'
## firstLine = 0
## lastLine = 138
## objClassName = 'RSG'
## refName = 'Smith'

#fileName = '/users/badenes/observations/LMC/lmc_dtd_0.0.dat'
#firstLine = 0
#lastLine = 161
#objClassName = 'PNe0'
#refName = 'Reid&Parker'

print objClassName+' in the LMC.'
print 'Reading object catalog from file '+fileName
f = open(fileName,'r')
lineList = f.readlines()
f.close()
objListRA = []
objListDec = []
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

nObjs = len(objListRA)
print 'Objects found in file: ', nObjs

#Read SFH map
sfhFileName = 'lmc_sfh.dat'
print 'Reading SFH from file '+sfhFileName
nCells = 0
nAgeBins = 16
cellNames = []
cellCentersRA = []
cellCentersDec = []
sfhMap = []
with open(sfhFileName, 'r') as f:
    for i in range(17): f.readline()
    while True:
        line = f.readline() #Hyphens
        if not line : 
            print 'End of file reached after cell ', nCells
            break #Check EOF
        
        words = f.readline().split()
        cellNames.append(words[1])  #Cell name    
        #print 'Reading cell ' + cellNames[-1]
        words = f.readline().split() #Cell center
        cellCentersRA.append(float(words[1][0:2]) + float(words[2][0:2])/60.0)
        cellCentersDec.append(float(words[3][0:3]) - float(words[4][0:2])/60.0)
        line = f.readline() #Hyphens
        sfhCell = np.zeros([3,nAgeBins])
        for bin in range(nAgeBins) :
            floats = [float(x) for x in f.readline().split()]
            floatArr = np.asarray(floats)
            sfhCell[0,nAgeBins-1-bin] = floatArr[[1,4,7,10]].sum() #Best Fit
            sfhCell[1,nAgeBins-1-bin] = floatArr[[2,5,8,11]].sum() #Lower limit
            sfhCell[2,nAgeBins-1-bin] = floatArr[[3,6,9,12]].sum() #Upper limit
            f.readline() #Blank line

        sfhMap.append(sfhCell)
        nCells += 1
        
sfhMapArr = np.asarray(sfhMap)
sfhMapArr = sfhMapArr/1e6 #Appears to be SFH for 'x' number of cells, 3 columns (median, upper, lower)limit and in different agebins        
cellCentersRA = np.asarray(cellCentersRA)
cellCentersDec = np.asarray(cellCentersDec)

#Modified cell centers: read from file
cellCentersFileName = 'LMC_SFH_Cell_Centers_Corrected.txt'
undivCells = [item for item in range(nCells) if len(cellNames[item]) == 2]
divCells = [item for item in range(nCells) if len(cellNames[item]) == 5]
nColumns = 24
nRows = 19
with open(cellCentersFileName, 'r') as f:
    for cell in range(nCells) :
        words = f.readline().split()
        cellIndex = cellNames.index(words[0])
        cellCentersRA[cellIndex] = float(words[1])
        cellCentersDec[cellIndex] = float(words[2])
        ## if len(words[0]) == 2 :
        ##     verticesRA = cellCentersRA[cellIndex] + 0.5*raInc*np.array([-1,-1,1,1]) 
        ##     verticesDec = cellCentersDec[cellIndex] + 0.5*decInc*np.array([-1,1,1,-1]) 
            #fig.show_rectangles(cellCentersRA[cellIndex],cellCentersDec[cellIndex],raInc,decInc,color='blue')
            # else :
            # fig.show_rectangles(cellCentersRA[cellIndex],cellCentersDec[cellIndex],0.5*raInc,0.5*decInc,color='blue')
        # print cellNames[cellIndex]


#~~~~~~~~~~~ Not sure what this part of the code is doing~~~~~~~~~~~~~~~~~~~#
raInc = cellCentersRA[3]-cellCentersRA[0] #BA-AA
decInc = cellCentersDec[1]-cellCentersDec[0] #AB-AA
raOrigin = 5.5-13.0*raInc
decOrigin = -72-0.75*decInc
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#Relevant variables - cellNames, cellCentersRA, cellCentersDec - these are all in degrees.
fieldCenter = SkyCoord('05h22m00s', '-69d00m00s', frame=coord.FK5(equinox=Time('B1950', scale='utc')))
fieldCenter.transform_to(coord.FK5(equinox=Time(2000, format='jyear', scale='utc')))
cellNames = np.asarray(cellNames)
surveyExtent = 5.25/2.0 #degrees
print fieldCenter.ra.deg, fieldCenter.dec.deg

is_RAinSurvey = (cellCentersRA >= (fieldCenter.ra.deg - surveyExtent)) & (cellCentersRA <= (fieldCenter.ra.deg + surveyExtent))
cellCentersRA = cellCentersRA[is_RAinSurvey]  
cellCentersDec = cellCentersDec[is_RAinSurvey]
cellNames = cellNames[is_RAinSurvey]
is_DecinSurvey = (cellCentersDec >= (fieldCenter.dec.deg - surveyExtent)) & (cellCentersDec <= (fieldCenter.dec.deg + surveyExtent))
cellCentersRA_inSurvey = cellCentersRA[is_DecinSurvey]  
cellCentersDec_inSurvey = cellCentersDec[is_DecinSurvey]
cellNames_inSurvey = cellNames[is_DecinSurvey]

print cellNames_inSurvey
