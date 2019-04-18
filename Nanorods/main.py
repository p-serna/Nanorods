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

wdirs = ["/mnt/data/Anastasia/18_11_29_pd23_11_div6_25Hzsqwave/",
"/mnt/data/Anastasia/QDs_prelabsem1811/","/mnt/data/Anastasia/Initial/","/mnt/data/Anastasia/Glass/",
"/mnt/data/Anastasia/test_18_07_11_pd29_06_div12_WISSNR","/mnt/data/Anastasia/19_02_20_pd15_02_div5_NR_BeRST_MOVILITYONLY"]

dfiles = []
for dirt in wdirs:
    basedir = dirt
    files = os.listdir(basedir)
    if dirt[-1] != '/':
        dirt = dirt+'/'
    for f in files:
        if f[-4:]=='.tif': 
            try:
                i = int(f[-5])
                dfiles.append(dirt+f)
            except:
                pass

for cfile in dfiles:
    #   1. Extracting ROIs
    os.system('CUDA_VISIBLE_DEVICES="" python3 extractingPeaks_general.py '+cfile) 
    #   1.1 Selecting ROIS inside the mask
    os.system('python3 inMask.py '+cfile) 

    #   2. Producing a movie with selected ROIs
    os.system('python3 onlymovieAll.py '+cfile)
    wdir = ''
    cf2 = cfile.split(".")[0].split("/")
    for fs in cf2[:-1]:
        wdir = wdir+fs+'/'
    wdir = wdir+cf2[-1]+'output/'

    #   2.1 Producing ROIS gifs!
    os.system('python3 rois_GIFs.py '+wdir) 

    #     3.1. Fourier transform and Spectrogram    
    fs = 100
    bt = 25
    bw = 4
    os.system('python3 signal_analysis/Spectsearch.py '+str(fs)+' '+str(bt)+' '+str(bw)+' '+wdir)
    # ~ os.system('CUDA_VISIBLE_DEVICES="" python3 extractingPeaks_general.py '+cfile) 
    #     3.2. Selecting blinking periods - statistics ON state 
    # ~ os.system('CUDA_VISIBLE_DEVICES="" python3 extractingPeaks_general.py '+cfile) 
    #     4.1. Fits to 2d gaussians
    os.system('python3 ROIS_analysis_all_speeding.py '+wdir) 
    #     4.2. Drift estimate
    wdir = ''
    cf2 = cfile.split(".")[0].split("/")
    for fs in cf2[:-1]:
        wdir = wdir+fs+'/'
    wdir = wdir+cf2[-1]+'output/sptrack/'
    os.system('python3 sptrack/Correct_drift.py '+wdir) 
    #     4.3. MSD estimate
    # ~ wdir = ''
    # ~ cf2 = cfile.split(".")[0].split("/")
    # ~ for fs in cf2[:-1]:
        # ~ wdir = wdir+fs+'/'
    # ~ wdir = wdir+cf2[-1]+'output/sptrack/'
    os.system('python3 sptrack/msdestimate_longer.py '+wdir) 
    #     4.4. MSD analysis
    #     4.5. Laterality
    
