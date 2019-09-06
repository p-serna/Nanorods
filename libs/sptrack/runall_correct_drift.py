import sys
import os

wdirs = ["/mnt/data/Anastasia/18_11_29_pd23_11_div6_25Hzsqwave/",
"/mnt/data/Anastasia/QDs_prelabsem1811/","/mnt/data/Anastasia/Initial/","/mnt/data/Anastasia/Glass/",
"/mnt/data/Anastasia/test_18_07_11_pd29_06_div12_WISSNR","/mnt/data/Anastasia/Glass"]

wdirs = ["/mnt/data/Anastasia/18_12_10_pd3_12_div7_WIS_NR-BeRST","/mnt/data/Anastasia/18_12_12_pd7_12_div5_WIS_NR-BeRST",
"/mnt/data/Anastasia/18_12_12_pd7_12_div5_WIS_NR-BeRST_DM590","/mnt/data/Anastasia/19_01_30_pd25_01_div5_NR_BeRST",
"/mnt/data/Anastasia/19_02_05_pd1_02_div4_NR_BeRST"]

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
    wdir = ''
    cf2 = cfile.split(".")[0].split("/")
    for fs in cf2[:-1]:
        wdir = wdir+fs+'/'
    wdir = wdir+cf2[-1]+'output/sptrack/'
    os.system('python3 Correct_drift.py '+wdir) 
    
