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
if False:
   fileName = DTDpath+'InputFiles/OGLE_LMC_ClassCeph.txt'
   objClassName = 'Cepheids'
   refName = 'OGLE'

if False:
   fileName = DTDpath+'InputFiles/OGLE_LMC_Type2Ceph.txt'
   objClassName = 'Type2_Cepheids'
   refName = 'OGLE'

if False:
   fileName = DTDpath+'InputFiles/OGLE_LMC_AnomCeph.txt'
   objClassName = 'Anom_Cepheids'
   refName = 'OGLE'

if False:
   fileName = DTDpath+'InputFiles/OGLE_LMC_DeltaScuti.txt'
   objClassName = 'DeltaScuti'
   refName = 'OGLE'

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

if True:
    fileName = DTDpath + 'InputFiles/OGLE_LMC_RRLyrae_noOutliers.txt'
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

if False:
    fileName = DTDpath + 'InputFiles/lmc_dtd.dat'
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


#NOTE - if you're downloading data from OGLE, remove the '#' in the .txt file in the 
#headings row where 'ID', 'RA', 'Dec' etc is.

obj_subtype = 'All'
objListRA, objListDec = dtdutils.object_ra_decs(fileName, objClassName, obj_subtype=obj_subtype)

#Read SFH map

sfhFileName = DTDpath + 'MC_SFH_Maps/lmc_sfh.dat'
outPathName = DTDpath + 'Output_SFH_Files/'
plotPathName = DTDpath + 'Writeup/'
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

print cellCentersRA[cellInSurvey].shape
print cellCentersDec[cellInSurvey].shape


line_list = [np.array([[cellCentersRA[0], cellCentersRA[-1]], [cellCentersDec[0], cellCentersDec[0]]])]
#If requested, read LMC image and produce a map of objects
print '\n\nMaking object Map...'
if True :
    fig = aplpy.FITSFigure(DTDpath + 'InputFiles/LMC60.M0NHC.FITS')
    fig.set_theme('publication')
    fig.recenter(5.35*15.0,-68.75,radius=4.5)
    fig.set_axis_labels_font(size=16)
    fig.set_tick_labels_font(size=16)
    fig.add_grid()
    fig.grid.set_linestyle('dotted')
    #fig.grid.set_yspacing(0.2)
    fig.grid.set_color('black')
    fig.show_grayscale()
    fig.show_markers(15.0*np.asarray(objListRA), np.asarray(objListDec), edgecolor='red', marker='.', s=1)
    fig.add_label(15.0*5.75,-65.5, objClassName+' from '+refName, color = 'black',size = 16)
    #fig.add_label(15.0*5.6,-65.5, objClassName+' from Reid & Parker (2010)', color = 'black',size = 16)
    fig.tick_labels.set_yformat('dd')
    fig.tick_labels.set_xformat('hh:mm')
#    fig.show_markers(cellCentersRA,cellCentersDec,marker='o',edgecolor='blue')
#    fig.show_markers(cellCentersRA[cellInSurvey],cellCentersDec[cellInSurvey],marker='o',edgecolor='green')
#    fig.show_rectangles(xw=cellCentersRA, yw=cellCentersDec, width=0.2, height=0.2, color='k')
#Overlay Reid&Parker coverage for PNe 
    if refName == 'Reid&Parker' :  #Made some changes to this codeblock. Using a combination of SkyCoord and transform_to instead
        # of coord.FK5 and precess_to()
        fieldCenter = SkyCoord('05h22m00s', '-69d00m00s', frame=coord.FK5(equinox=Time('B1950', scale='utc')))
        fieldCenter.transform_to(coord.FK5(equinox=Time(2000, format='jyear', scale='utc')))
        #print fieldCenter.ra.deg, fieldCenter.dec.deg
        fig.show_rectangles(fieldCenter.ra.deg, fieldCenter.dec.deg, 5.25, 5.25, color = 'orange')

    fig.save(outPathName+'LMC_'+objClassName+'_'+obj_subtype+'_Map.pdf')


#Bin SFH map and write final file
logAgeArr = np.asarray([6.8, 7.1, 7.4, 7.7, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0, 10.2])
ageArr = 10.0**logAgeArr
logAgeLimsArr = np.asarray([6.65,6.95,7.25,7.55,7.85,8.1,8.3,8.5,8.7,8.9,9.1,9.3,9.5,9.7,9.9,10.1,10.3])
ageLimsArr = 10.0**logAgeLimsArr
ageIntervalsArr = ageLimsArr[1:]-ageLimsArr[:-1]
#ageIntervalsArr = np.asarray([4101224.2, 9404646.0, 18764724., 37440568., 54185272., 75594208.,
#                              1.1980914e+08, 1.8988461e+08, 3.0094621e+08, 4.7696742e+08, 7.5594214e+08, 
#                              1.1980915e+09, 1.8988457e+09, 3.0094623e+09, 4.7696753e+09, 5.9120108e+09])

binningSchemeNames = ['Coarse','Massive','Medium','MediumB','MediumC','Subfine','Fine','CrazyCeph','RRLyraeLB','Unbinned']
nSchemes = len(binningSchemeNames)
#binsList = [[0, 1, 2], [3, 4, 5, 6, 7, 8], [9, 10], [11, 12], [13], [14, 15]] 
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
for scheme in range(nSchemes) :
    nBinsScheme = len(binningSchemes[scheme])
    sfhMapBinned = np.zeros((nBinsScheme,3,nCells))
    for bin in range(nBinsScheme) :
        sfhMapBinned[bin,0,:] =  np.reshape(np.dot(sfhMapArr[:,0,[(binningSchemes[scheme])[bin]]],ageIntervalsArr[(binningSchemes[scheme])[bin]]),nCells)
        sfhMapBinned[bin,1,:] =  np.reshape(np.dot(sfhMapArr[:,1,[(binningSchemes[scheme])[bin]]],ageIntervalsArr[(binningSchemes[scheme])[bin]]),nCells)
        sfhMapBinned[bin,2,:] =  np.reshape(np.dot(sfhMapArr[:,2,[(binningSchemes[scheme])[bin]]],ageIntervalsArr[(binningSchemes[scheme])[bin]]),nCells)
        #sfhMapBinned[bin,2,:] =  np.reshape(np.dot(sfhMapArr[:,2,[binsList[bin]]],ageIntervalsArr[[binsList[bin]]]),nCells)

    sfhMapBinnedList.append(sfhMapBinned)

#Plot maps
createMaps = True
if createMaps :                    
    #Create map grid
    gridKey = np.zeros([2*nColumns,2*nRows],dtype=int)
    gridKey[:,:] = -1
    for cell in range(nCells) :
        colLetter = cellNames[cell][0]
        rowLetter = cellNames[cell][1]
        if len(cellNames[cell]) == 5 :
            colNum = int(cellNames[cell][3])
            rowNum = int(cellNames[cell][4])
            gridKey[2*(string.uppercase.find(colLetter))+colNum,2*(string.uppercase.find(rowLetter))+rowNum] = cell
            #print cellNames[cell], 2*(string.uppercase.find(colLetter))+colNum, 2*(string.uppercase.find(rowLetter))+rowNum
        else :
            gridKey[2*(string.uppercase.find(colLetter)),2*(string.uppercase.find(rowLetter))] = cell
            gridKey[2*(string.uppercase.find(colLetter))+1,2*(string.uppercase.find(rowLetter))] = cell
            gridKey[2*(string.uppercase.find(colLetter)),2*(string.uppercase.find(rowLetter))+1] = cell
            gridKey[2*(string.uppercase.find(colLetter))+1,2*(string.uppercase.find(rowLetter))+1] = cell
            #print cellNames[cell], 2*(string.uppercase.find(colLetter)), 2*(string.uppercase.find(rowLetter))   

   
    ## #Map the sfh
    sfhMapImages = np.zeros([nAgeBins,2*nColumns,2*nRows])
    sfhMapMassiveImages = np.zeros([4,2*nColumns,2*nRows])
    sfhMapRRLyraeLBImages = np.zeros([6,2*nColumns, 2*nRows])
    for row in range(2*nRows) :
        for column in range(2*nColumns) :
            if gridKey[column,row] >= 0 :
                sfhMapImages[:,column,row] = (sfhMap[gridKey[column,row]])[0]
                sfhMapMassiveImages[:,column,row] = (sfhMapBinnedList[1])[:,0,gridKey[column,row]]
                sfhMapRRLyraeLBImages[:,column,row] = (sfhMapBinnedList[-2])[:,0,gridKey[column,row]]

    plt.rcdefaults()
    plt.rcParams.update({'figure.autolayout':'True'})
    plt.rcParams.update({'font.size': 8})
    plt.rcParams.update({'mathtext.default':'regular'})
    plt.rcParams.update({'mathtext.fontset':'stixsans'})
    plt.rcParams.update({'axes.linewidth': 0.5})
    plt.rcParams.update({'xtick.major.size': 2})
    plt.rcParams.update({'xtick.major.width': 0.25 })
    plt.rcParams.update({'xtick.minor.size': 1})
    plt.rcParams.update({'xtick.minor.width': .25 })
    plt.rcParams.update({'ytick.major.size': 2})
    plt.rcParams.update({'ytick.major.width': .25 })
    plt.rcParams.update({'ytick.minor.size': 1})
    plt.rcParams.update({'ytick.minor.width': .25 })

#Plot unbinned, no objects
print 'Making plots :- \n\n'
print 'Unbinned Map...'
if False:
    plotFileName = outPathName+'LMC_Unbinned_Maps.pdf'
    plt.figure(1,figsize = [11.0, 8.5])
    plt.clf()
    f, axArr = plt.subplots(4, 4, sharex = True, sharey = True, squeeze = True)   
    for ageBin in range(nAgeBins) :
        row = int(np.floor(ageBin/4))
        column = ageBin-4*row
        axArr[row,column].set_xlim([2*nColumns,0])
        axArr[row,column].set_ylim([0,2*nRows])
        axArr[row,column].axis('off')
        axArr[row,column].imshow(np.transpose(sfhMapImages[ageBin,:,:]),origin='lower',interpolation='none')
        axArr[row,column].text(47,33,'Bin '+str(ageBin)+' log(t)='+str(logAgeArr[ageBin]),color='w')
    plt.savefig(plotFileName)

heatmap, xedges, yedges = np.histogram2d((15.0*np.asarray(objListRA)-cellCentersRA[0])/(0.5*raInc),(np.asarray(objListDec)-cellCentersDec[0])/(0.5*decInc), bins=30) 
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

colorscheme = ['#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026']

 #Plot unbinned
print 'Unbinned Map with objects...'
if True:
    plotFileName = plotPathName+'LMC_'+objClassName+'_'+obj_subtype+'_'+refName+'_Unbinned_Maps.pdf'
    plt.figure(1,figsize = [11.0, 8.5])
    plt.clf()
    f, axArr = plt.subplots(4, 4, sharex = True, sharey = True, squeeze = True) 
    f.subplots_adjust(hspace=-0.1, wspace=-0.1)
    for ageBin in range(nAgeBins) :
        row = int(np.floor(ageBin/4))
        column = ageBin-4*row
        row, column = 3-row, 3-column
        axArr[row,column].set_xlim([2.*nColumns,0])
        axArr[row,column].set_ylim([0,2.*nRows])
        axArr[row,column].axis('on')
        axArr[row,column].imshow(np.transpose(sfhMapImages[ageBin,:,:]), origin='lower', interpolation='none', cmap = 'Greys')
        fact1, unit1 = (1.0e6, ' Myrs') if logAgeLimsArr[ageBin]<9. else (1.0e9, ' Gyrs')
        age1 = np.around(10**logAgeLimsArr[ageBin]/fact1, decimals=1)
        fact2, unit2 = (1.0e6, ' Myrs') if logAgeLimsArr[ageBin+1]<9. else (1.0e9, ' Gyrs')
        age2 = np.around(10**logAgeLimsArr[ageBin+1]/fact2, decimals=1)

        axArr[row,column].text(0.1, 0.9, str(age1)+' '+unit1+' - '+str(age2)+' '+unit2, fontsize=6, transform=axArr[row,column].transAxes)
#        axArr[row,column].plot((15.0*np.asarray(objListRA)-cellCentersRA[0])/(0.5*raInc),(np.asarray(objListDec)-cellCentersDec[0])/(0.5*decInc),'r.',ms=0.5, alpha=0.1)
        axArr[row, column].contour(heatmap.T, extent=extent, levels=[20, 40, 70, 100, 150, 200, 250], linewidths=1.0, colors=colorscheme, origin='lower')
    plt.savefig(plotFileName)

##RRLyraeLB Map
print 'RRLyraeLB Map with objects...'
if True:
    plotFileName = outPathName+'LMC_'+objClassName+'_'+obj_subtype+'_'+refName+'_RRLyraeLB_Maps.pdf'
    plt.figure(1,figsize = [11.0, 8.5])
    plt.clf()
    scheme = -2
    print binningSchemes[scheme][0][-1]
    nBinsScheme = len(binningSchemes[scheme])
    f, axArr = plt.subplots(2, 3, sharex = True, sharey = True, squeeze = True) 
    f.tight_layout()  
    for ageBin in range(nBinsScheme) :
        row = int(np.floor(ageBin/3))
        column = ageBin-3*row
        axArr[row,column].set_xlim([2.*nColumns,0])
        axArr[row,column].set_ylim([0,2.*nRows])
        axArr[row,column].axis('off')
        axArr[row,column].imshow(np.transpose(sfhMapRRLyraeLBImages[ageBin,:,:]), origin='lower', interpolation='none', cmap = 'Greys')
        ageBin_correct = binningSchemes[scheme][ageBin][-1]
        fact, unit = (1.0e6, ' Myrs') if logAgeArr[ageBin_correct]<9. else (1.0e9, ' Gyrs')
        age = np.around(10**logAgeArr[ageBin_correct]/fact, decimals=1)

        axArr[row,column].text(0.1, 0.9, 'Age = '+str(age)+unit, fontsize=6, transform=axArr[row,column].transAxes)
        #axArr[row,column].plot((15.0*np.asarray(objListRA)-cellCentersRA[0])/(0.5*raInc),(np.asarray(objListDec)-cellCentersDec[0])/(0.5*decInc),'r.',ms=0.5)
        axArr[row, column].contour(heatmap.T, extent=extent, levels=[20, 60, 100, 200], linewidths=1.0, colors=colorscheme, origin='lower')
    plt.savefig(plotFileName)

#Plothistograms
print 'Unbinned Map histograms...'
if True:
    plotFileName = outPathName+'LMC_'+objClassName+'_'+obj_subtype+'_'+refName+'_Unbinned_Histograms.pdf'
    plt.figure(1,figsize = [11.0, 8.5])
    plt.clf()
    f, axArr = plt.subplots(4, 4, sharex = True, sharey = True, squeeze = True) 
    f.tight_layout()  
    for ageBin in range(nAgeBins) :
        row = int(np.floor(ageBin/4))
        column = ageBin-4*row
        bins = np.logspace(0, 4.2, 20)
        axArr[row,column].hist(np.transpose(sfhMapImages[ageBin,:,:]).flatten(), bins=bins, histtype='step', color='b')
        axArr[row,column].text(0.1, 0.9, 'Bin '+str(ageBin)+' log(t)='+str(logAgeArr[ageBin]), fontsize=4, transform=axArr[row,column].transAxes)
        axArr[row,column].set_xscale('log')
        axArr[row,column].set_xlim(0., 2.0e4)
        axArr[row,column].set_ylim(0., 500)
        #axArr[row,column].plot((15.0*np.asarray(objListRA)-cellCentersRA[0])/(0.5*raInc),(np.asarray(objListDec)-cellCentersDec[0])/(0.5*decInc),'r.',ms=0.5)
    plt.savefig(plotFileName)


 #Plot Massive
print 'Massive Maps...'
if False:
    plotFileName = outPathName+'LMC_'+objClassName+'_'+obj_subtype+'_'+refName+'_Massive_Maps.pdf'
    scheme = 1
    nBinsScheme = len(binningSchemes[scheme])
    plt.figure(1,figsize = [11.0, 8.5])
    plt.clf()
    f, axArr = plt.subplots(2, 2, sharex = True, sharey = True, squeeze = True)   
    for ageBin in range(nBinsScheme) :
         row = int(np.floor(ageBin/2))
         column = ageBin-2*row
         #print row, column
         axArr[row,column].set_xlim([2*nColumns,0])
         axArr[row,column].set_ylim([0,2*nRows])
         axArr[row,column].axis('off')
         axArr[row,column].imshow(np.transpose(sfhMapMassiveImages[ageBin,:,:]),cmap='Greys',origin='lower',interpolation='none')
         axArr[row,column].set_title('Bin '+str(ageBin)+' log(t)='+str(logAgeArr[ageBin]), fontsize=6)
         axArr[row, column].contour(heatmap.T, extent=extent, levels=[20, 60, 100, 200], colors=colorscheme, origin='lower')
#axArr[row,column].plot((15.0*np.asarray(objListRA)-cellCentersRA[0])/(0.5*raInc),(np.asarray(objListDec)-cellCentersDec[0])/(0.5*decInc),'r.',ms=1.0)
    
    plt.savefig(plotFileName)


#Plot Correlations
if False:
    plotFileName = outPathName+'LMC_'+objClassName+'_'+obj_subtype+'_'+refName+'_Corr.pdf'
    scheme = 2
    print binningSchemeNames[scheme]
    nBinsScheme = len(binningSchemes[scheme])
    plt.figure(1,figsize = [11.0, 8.5])
    plt.clf()
    f, axArr = plt.subplots(2, 3, squeeze = True)   
    for scheme in range(nSchemes) :
        row = int(np.floor(scheme/3))
        column = scheme-3*row
        print row, column
       # axArr[row,column].text(0.5,0.95*(corrList[scheme]).max(),binningSchemeNames[scheme])
       # axArr[row,column].bar(range(len(corrList[scheme])),corrList[scheme])  
        axArr[row,column].text(0.5,0.95,binningSchemeNames[scheme])
       
    
    plt.savefig(plotFileName)

print 'Done!'
