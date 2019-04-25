import os
import warning
import numpy as np
from scipy.stats import  linregress

def EMCCD_CMOS_calib(folder = './',csv = True, ROIsize=9):

    #Extracting images and csv files
    files_tif = []
    files_csv = []
    # ~ filenames = []
    for file in os.listdir(folder):
        if file.endswith('.tif'):
            files_tif.append(file)
            # ~ filenames.append(file.split('_')[0])
        elif file.endswith('.csv'):
            files_csv.append(file)
            
    
    # ~ filenames = list(set(filenames))
    # ~ filenames.sort()
    files_tif.sort()
    files_csv.sort()
    #files = {name: [] for name in filenames}
        
    # Let us subdivide files with name of the cell
    # ~ for file in files_tif:
        # ~ filename = file.split('_')[0]
        # ~ files[filename].append(file)
    
    files = files_tif
    filesCMOS = [file for file in files if file.find('CMOS')>=0]
    if len(filesCMOS) == 0:
        raise 'Sorry, there are no files with CMOS in their name.'
    
    
    filesEMCCD = [file for file in files if file.find('EMCCD')>=0]
    if len(filesEMCCD)==0:
        raise 'We stop, no files with EMCCD in their name, for '+name+'.'
        
    
    ptsC = []
    ptsE = []
    for i,names in enumerate(zip(filesCMOS,filesEMCCD)): #real traces
        fCMOS, fEMCCD = names
    
        # Change these lines to appropiate reading of tif files:
        CMOS = imread(folder+fCMOS)
        EMCCD= imread(folder+fEMCCD)
    
        CMOS = CMOS[:,:CMOS.shape[1]//2]
    
        if not csv:
            # Produce csv - pass it to subroutine    
            pass
        else:
            fCMOScsv, fEMCCDcsv =  (filesCMOScsv[i],filesEMCCDcsv[i])
            # ReadROIS!? Show me csv files
            CMOS_pts = ReadROIs(loadtxt(fCMOS))
            EMCCD_pts  = ReadROIs(loadtxt(fEMCCD))
        
        #CenterROIs function
        CMOS_pts_center=CenterROIs(CMOS,CMOS_pts,ROIsize);
        EMCCD_pts_center=CenterROIs(EMCCD,EMCCD_pts,ROIsize);
        ptsC.append(CMOS_pts_center)
        ptsE.append(EMCCD_pts_center)

    ptsC = np.array(ptsC)
    ptsE = np.array(ptsE)
    Xc=ptsC[:,0]
    Yc=ptsC[:,1]
    Xe=ptsE[:,0]
    Ye=ptsE[:,1]

    fitobject = linregress(Xc,Xe)
    #slope, intercept, r_value, p_value, std_err
    ax, bx, _, _, _ = fitobject

    fitobject = linregress(Yc,Ye)
    #slope, intercept, r_value, p_value, std_err
    ax, bx, _, _, _ = fitobject

    return(ax,ay,bx,by)
