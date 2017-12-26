import numpy as np                        ##IMPORTS
import os
import glob
import aplpy
import string
import sys
import dtdutils
reload(dtdutils)
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy import units as u
from astropy import wcs 
from astropy.io import fits
from astropy.io import ascii 
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy import coordinates as coordf
from astropy.time import Time



#Read Generic catalogs
DTDpath = '/Users/sumits2k/Desktop/Research/SNResearch2/RadioSNRs/DTD/'

if True:
    fileName = DTDpath + 'InputFiles/OGLE_LMC_RRLyrae_noOutliers.txt'
    firstLine = 7
    lastLine = 24913
    objClassName = 'RRLyrae'
    refName = 'OGLE'

#NOTE - if you're downloading data from OGLE, remove the '#' in the .txt file in the 
#headings row where 'ID', 'RA', 'Dec' etc is.

obj_subtype = 'All'
objListRA, objListDec = dtdutils.object_ra_decs(fileName, objClassName, obj_subtype=obj_subtype)

#Read SFH map

sfhFileName = DTDpath + 'MC_SFH_Maps/lmc_sfh.dat'
outPathName = DTDpath + 'Output_SFH_Files/MixingTest/'

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
cellCentersFileName = DTDpath + 'MC_SFH_Maps/LMC_SFH_Cell_Centers_Corrected.txt'
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
    surveyCells = [line.strip() for line in open(DTDpath + 'MC_SFH_Maps/LMC_SFH_Cells_OGLE.txt')]
elif refName == 'Reid&Parker' :
    surveyCells = [line.strip() for line in open(DTDpath + 'MC_SFH_Maps/LMC_SFH_Cells_Reid&Parker.txt')]
#    surveyCells = cellsinSFHSurvey(cellNames, cellCentersRA, cellCentersDec)
else : 
    surveyCells = cellNames
    
nSurveyCells = len(surveyCells)
objRAArr = 15.0*np.asarray(objListRA)
objDecArr = np.asarray(objListDec)
cellInSurvey = np.zeros(nCells,dtype='bool')
nObjAcc = 0
objMap = np.zeros(nCells)
for cell in range(nCells) :
    #Determine if cell is within survey area
    if cellNames[cell] in surveyCells :
        cellInSurvey[cell] = True
        if len(cellNames[cell]) == 2 :
            cellHalfSizeRA = raInc/2.0
            cellHalfSizeDec = decInc/2.0
        else :
            cellHalfSizeRA = raInc/4.0
            cellHalfSizeDec = decInc/4.0
            
        objInCell = np.logical_and(abs(objRAArr - cellCentersRA[cell]) <= cellHalfSizeRA, 
                                   abs(objDecArr - cellCentersDec[cell]) <= cellHalfSizeDec)
        
        nObjInCell = objInCell.sum()
        nObjAcc += nObjInCell
        objMap[cell] = nObjInCell
            
print 'HZ cells in survey area: ', cellInSurvey.sum()
print 'Total objects: ', nObjAcc, objMap.sum()


#Bin SFH map and write final file
logAgeArr = np.asarray([6.8, 7.1, 7.4, 7.7, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0, 10.2])
ageArr = 10.0**logAgeArr
logAgeLimsArr = np.asarray([6.65,6.95,7.25,7.55,7.85,8.1,8.3,8.5,8.7,8.9,9.1,9.3,9.5,9.7,9.9,10.1,10.3])
ageLimsArr = 10.0**logAgeLimsArr
ageIntervalsArr = ageLimsArr[1:]-ageLimsArr[:-1]
#ageIntervalsArr = np.asarray([4101224.2, 9404646.0, 18764724., 37440568., 54185272., 75594208.,
#                              1.1980914e+08, 1.8988461e+08, 3.0094621e+08, 4.7696742e+08, 7.5594214e+08, 
#                              1.1980915e+09, 1.8988457e+09, 3.0094623e+09, 4.7696753e+09, 5.9120108e+09])

#binningSchemeNames = ['Coarse','Massive','Medium','MediumB','MediumC','Subfine','Fine','CrazyCeph','RRLyraeLB','Unbinned']
binningSchemeNames = ['Unbinned']
nSchemes = len(binningSchemeNames)
#binsList = [[0, 1, 2], [3, 4, 5, 6, 7, 8], [9, 10], [11, 12], [13], [14, 15]] 
binningSchemes = [range(16)]
if 0:
    binningSchemes = [[[0,1,2,3],[4,5,6,7,8],[9,10,11,12,13,14,15]]
                  ,[[0],[1],[2],[3,4,5,6,7,8,9,10,11,12,13,14,15]]
                  ,[[0,1,2],[3,4,5,6,7,8],[9,10],[11],[12,13],[14,15]]
                  ,[[0,1,2],[3,4,5,6,7,8],[9,10],[11,12],[13],[14,15]]
                  ,[[0,1,2],[3,4,5,6],[7,8,9,10],[11,12],[13],[14,15]]
                  ,[[0,1,2],[3,4,5],[6,7,8],[9,10],[11,12],[13],[14,15]]
                  ,[[0,1,2],[3,4,5],[6,7,8],[9,10],[11],[12],[13],[14,15]]
                  ,[[0],[1],[2,3],[4],[5,6,7],[8,9],[10,11,12], [13,14,15]]
                  ,[[0,1,2,3,4,5,6,7,8,9],[10],[11],[12],[13,14],[15]]
                  ,range(16)]


corrList = []
corrListNorm = []
sfhMapBinnedList = []
sfhMapTotal = np.reshape(np.dot(sfhMapArr[:,0,:],ageIntervalsArr),nCells)

print '\n\nCreating SFH files ending in .dat\n'
for scheme in range(nSchemes):
    print binningSchemeNames[scheme]
    nBinsScheme = len(binningSchemes[scheme])
    sfhMapBinned = np.zeros((nBinsScheme,3,nCells))
    sfhMapBinned_scaled = np.zeros((nBinsScheme, nCells))
    for bin in range(nBinsScheme) :
        sfhMapBinned[bin,0,:] =  np.reshape(np.dot(sfhMapArr[:,0,[(binningSchemes[scheme])[bin]]],ageIntervalsArr[(binningSchemes[scheme])[bin]]),nCells)
        sfhMapBinned[bin,1,:] =  np.reshape(np.dot(sfhMapArr[:,1,[(binningSchemes[scheme])[bin]]],ageIntervalsArr[(binningSchemes[scheme])[bin]]),nCells)
        sfhMapBinned[bin,2,:] =  np.reshape(np.dot(sfhMapArr[:,2,[(binningSchemes[scheme])[bin]]],ageIntervalsArr[(binningSchemes[scheme])[bin]]),nCells)
        #sfhMapBinned[bin,2,:] =  np.reshape(np.dot(sfhMapArr[:,2,[binsList[bin]]],ageIntervalsArr[[binsList[bin]]]),nCells)

    #Using scaled sfhMapBinned[:, 0, :]
    f = 0.7
    M_j = [np.sum(sfhMapBinned[j, 0, :]) for j in range(nBinsScheme)]
    M_i = [np.sum(sfhMapBinned[:, 0, i]) for i in range(nCells)]
    M_tot = np.sum(sfhMapBinned[:, 0, :])
    print 'M_j = ', M_j
    print 'Total Mass (M_tot) = ', M_tot
    print 'Total Mass (Sum of M_j) = ', np.sum(M_j)
    print 'Total Mass (Sum of M_i) = ', np.sum(M_i)

    for j in range(nBinsScheme):
        print 'j = ', j
        for i in range(nCells):
            sfhMapBinned_scaled[j, i] = (1-f)*sfhMapBinned[j, 0, i] + f*((M_i[i]*M_j[j])/M_tot)
            
    sfhMapBinned[:, 0, :] = sfhMapBinned_scaled

    sfhMapBinnedList.append(sfhMapBinned)

    #Write SFH_Cells files
    #corrListScheme = np.zeros(nBinsScheme)
    #corrListSchemeNorm = np.zeros(nBinsScheme)
    with open(outPathName+'LMC_SFH_Cells_'+objClassName+obj_subtype+'_'+binningSchemeNames[scheme]+'_f_'+str(f)+'.dat', 'w') as f:
        for cell in range(nCells) :
            if cellNames[cell] in surveyCells :
                f.write('%s  %i ' % (cellNames[cell], objMap[cell]))
                f.write((nBinsScheme*'%0.3e  ') % tuple(sfhMapBinned[:,0,cell]))
                f.write('\n')
                f.write('         ')
                f.write((nBinsScheme*'%0.3e  ') % tuple(sfhMapBinned[:,1,cell]))
                f.write('\n')
                f.write('         ')
                f.write((nBinsScheme*'%0.3e  ') % tuple(sfhMapBinned[:,2,cell]))
                f.write('\n')
                #corrListScheme += objMap[cell]*sfhMapBinned[:,0,cell]
                #corrListSchemeNorm += sfhMapBinned[:,0,cell]
                #stellarMassFormed[:,0] += sfhMapBinned[:,0,cell]

                #corrList.append(corrListScheme/corrListSchemeNorm)

print 'Done!'
