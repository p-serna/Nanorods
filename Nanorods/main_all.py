import sys
import os

#-------------------------------------------------------------
#
# In this main script we will call to the different subroutines in order for each file,
# if preferred you can run each step for all the files at the same time by running runall_*
# scripts
#
#   1. Extracting ROIs
#   2. Producing a movie with selected ROIs
#   3. Analysis of blinking NRs
#     3.1. Fourier transform and Spectrogram
#     3.2. Selecting blinking periods - statistics ON state 
#     3.3. ??
#   4. Diffusion analysis
#     4.1. Fits to 2d gaussians
#     4.2. MSD fits
#     4.3. Laterality//anisotropy
#     4.4. ??


# First we need to set the working directories.


wdirs = ["18_11_29_pd23_11_div6_25Hzsqwave","18_12_10_pd3_12_div7_WIS_NR-BeRST","18_12_12_pd7_12_div5_WIS_NR-BeRST",
"18_12_12_pd7_12_div5_WIS_NR-BeRST_DM590","19_01_30_pd25_01_div5_NR_BeRST",
"19_02_05_pd1_02_div4_NR_BeRST"]


#os.system('python3 runall_peakext.py ') 
#os.system('python3 runall_ROIsfits.py ') 
#os.system('python3 runall_inMask.py ') 
#os.system('python3 sptrack/runall_correct_drift.py ') 
#os.system('python3 sptrack/runall_msdestimates.py ') 
homed = "/mnt/data/Anastasia/"

dfiles = []
for dirt in wdirs:
    basedir = homed+dirt
    files = os.listdir(basedir)
    if dirt[-1] != '/':
        dirt = dirt+'/'
    for f in files:
        if f[-4:]=='.tif': 
            try:
                i = int(f[-5])
                dfiles.append(homed+dirt+f)
            except:
                pass
                
for cfile in dfiles:
    os.system('CUDA_VISIBLE_DEVICES="" python3 extractingPeaks_general.py '+cfile) 
    
for cfile in dfiles:
    wdir = ''
    cf2 = cfile.split(".")[0].split("/")
    for fs in cf2[:-1]:
        wdir = wdir+fs+'/'
    wdir = wdir+cf2[-1]+'output/'
    os.system('python3 ROIS_analysis_all_speeding.py '+wdir) 
    
for cfile in dfiles:
    os.system('python3 onlymovieAll.py '+cfile) 
    
for cfile in dfiles:
    os.system('python3 inMask.py '+cfile) 
    

for cfile in dfiles:
    wdir = ''
    cf2 = cfile.split(".")[0].split("/")
    for fs in cf2[:-1]:
        wdir = wdir+fs+'/'
    wdir = wdir+cf2[-1]+'output/sptrack/'
    os.system('python3 sptrack/Correct_drift.py '+wdir) 
    os.system('python3 sptrack/msdestimate_longer.py '+wdir+' 1') 
    
