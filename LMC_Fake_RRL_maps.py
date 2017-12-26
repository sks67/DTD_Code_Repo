#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 
# DESCRIPTION: Generate Poisson-based fake RRL maps, and 
# redo the DTD method to assess if we are recovering the 
# correct signal. In this case, I am inserting 
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



import os
import numpy as np
DTD_path = os.getenv('DTD')+'/'
sad_files_path = DTD_path + 'Output_SFH_Files/'

#INPUT PARAMETERS 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
object_Name = 'RRLyrae'
object_Subtype = 'All'
binning_type = 'Unbinned'
dtd_old = np.concatenate(([0.]*14, [1.0e-5], [1.0e-5]))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

sad_file_name = sad_files_path + 'LMC_SFH_Cells_' + object_Name + object_Subtype + '_' + binning_type + '.dat'

#Open file containing stellar masses for a given binning strategy, and only use the 
#rows with best-fit stellar masses. 
sad_file_obj = open(sad_file_name, 'r')
sad_file = sad_file_obj.readlines()
sad_file_obj.close()

sad_cellnames = np.array([lines.split()[0] for lines in sad_file if len(lines.split())==18])
sads = np.array([map(float, lines.split()[2:]) if len(lines.split())==18 else map(float, lines.split()) \
                 for lines in sad_file])

#Lots of things happening here. I'm reading each line of the file with the for statement, picking only lines that have
#best fit stellar mass (i.e. with 18 entries, 1st two being cell name and number of objects in cell). The split() splits
#up the full file string, and map(float,...) converts all the SAD entries into numbers!

lambda_i = np.dot(sads[0::3], dtd_old)
print 100.*np.sum(lambda_i)/23459.  #Check if the fraction of old objects produced is consistent.
k_i = np.random.poisson(lam=lambda_i)

nCells = lambda_i.size
nBinsScheme = sads.shape[-1]
if 1:
    with open(sad_files_path+'Fake_LMC_SFH_Cells_'+object_Name+object_Subtype+'_'+binning_type+'.dat', 'w') as f:
        for cell in range(nCells) :
            f.write('%s  %i ' % (sad_cellnames[cell], k_i[cell]))
            f.write((nBinsScheme*'%0.3e  ') % tuple(sads[3*cell]))
            f.write('\n')
            f.write('         ')
            f.write((nBinsScheme*'%0.3e  ') % tuple(sads[3*cell + 1]))
            f.write('\n')
            f.write('         ')
            f.write((nBinsScheme*'%0.3e  ') % tuple(sads[3*cell + 2]))
            f.write('\n')
