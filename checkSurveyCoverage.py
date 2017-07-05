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

def cellsinSFHSurvey(cellNames, cellCentersRA, cellCentersDec):
    return cellNames
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

fileName = 'OGLE_LMC_RRLyrae.txt'
firstLine = 7
lastLine = 24913
objClassName = 'RRLyrae'
refName = 'OGLE'

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

##fileName = 'InputFiles/lmc_dtd.dat'
##firstLine = 0
##lastLine = 573
##objClassName = 'PNe'
##refName = 'Reid&Parker'

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


#Read SFH map
sfhFileName = 'InputFiles/lmc_sfh.dat'
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
    
#Write file with number of objects in each cell that falls within the survey area
if refName == 'OGLE' : 
    surveyCells = [line.strip() for line in open('LMC_SFH_Cells_OGLE.txt')]
elif refName == 'Reid&Parker' :
    surveyCells = [line.strip() for line in open('LMC_SFH_Cells_Reid&Parker.txt')]
#    surveyCells = cellsinSFHSurvey(cellNames, cellCentersRA, cellCentersDec)
else : 
    surveyCells = cellNames

nSurveyCells = len(surveyCells)
cellInSurvey = np.zeros(nCells,dtype='bool')
nObjAcc = 0
objMap = np.zeros(nCells)
for cell in range(nCells) :
    #Determine if cell is within survey area                                                                                                                  
    if cellNames[cell] in surveyCells :
        cellInSurvey[cell] = True

print 'HZ cells in survey area: ', cellInSurvey.sum()

#If requested, read LMC image and produce a map of objects
if True :
    fig = aplpy.FITSFigure('LMC60.M0NHC.FITS')
    fig.set_theme('publication')
    fig.recenter(5.35*15.0,-68.75,radius=4.5)
    fig.set_axis_labels_font(size=16)
    fig.set_tick_labels_font(size=16)
    fig.add_grid()
    fig.grid.set_linestyle('dotted')
    #fig.grid.set_yspacing(0.2)
    fig.grid.set_color('black')
    fig.show_grayscale()
    fig.add_label(15.0*5.75,-65.5, objClassName+' from '+refName, color = 'black',size = 16)
    #fig.add_label(15.0*5.6,-65.5, objClassName+' from Reid & Parker (2010)', color = 'black',size = 16)
    fig.tick_labels.set_yformat('dd')
    fig.tick_labels.set_xformat('hh:mm')
    fig.show_markers(cellCentersRA,cellCentersDec,marker='o',s=20.0,edgecolor='blue')
    fig.show_markers(cellCentersRA[cellInSurvey],cellCentersDec[cellInSurvey],marker='o',facecolor='blue')
    #Overlay Reid&Parker coverage for PNe 
    if refName == 'Reid&Parker' :  #Made some changes to this codeblock. Using a combination of SkyCoord and transform_to instead
        # of coord.FK5 and precess_to()
        fieldCenter = SkyCoord('05h22m00s', '-69d00m00s', frame=coord.FK5(equinox=Time('B1950', scale='utc')))
        fieldCenter.transform_to(coord.FK5(equinox=Time(2000, format='jyear', scale='utc')))
        print fieldCenter.ra.deg, fieldCenter.dec.deg
        fig.show_rectangles(fieldCenter.ra.deg, fieldCenter.dec.deg, 5.25, 5.25, color = 'orange')

    fig.save('LMC_'+objClassName+'_Map.pdf')

sys.exit([0])
