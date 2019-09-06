import sys
import os

wdirs = ["/mnt/data/Anastasia/18_11_29_pd23_11_div6_25Hzsqwave/",
"/mnt/data/Anastasia/18_12_10_pd3_12_div7_WIS_NR-BeRST/",
"/mnt/data/Anastasia/18_12_12_pd7_12_div5_WIS_NR-BeRST/",
"/mnt/data/Anastasia/18_12_12_pd7_12_div5_WIS_NR-BeRST_DM590/",
"/mnt/data/Anastasia/19_01_30_pd25_01_div5_NR_BeRST/",
"/mnt/data/Anastasia/19_02_05_pd1_02_div4_NR_BeRST/"]


dfiles = []
for dirt in wdirs:
    basedir = dirt
    files = os.listdir(basedir)
    if dirt[-1] != '/':
        dirt = dirt+'/'
    files.sort()
    for f in files:
        if f[-4:]=='.tif': 
            try:
                i = int(f[-5])
                dfiles.append(dirt+f)
            except:
                pass

sys.path.append("/users/bssn/serna/GitIBENS/Nanorods/")

from numpy import *
from matplotlib.pylab import *
import sys
import os
import cv2
from tools.extracting_ROIfrommovie import *
import pickle

based = "/mnt/data/Anastasia/Experimental_data/random_traces/"
files = os.listdir(based)

Xs, Ys, dXs, dYs = [], [], [], []
files.sort()
for f in files:
    if f[:5] == "DataX":
        print("X:",f)
        Xs.append(load(f))
    elif f[:5] == "DataY":
        print("Y:",f)
        Ys.append(load(f))
    elif f[:6] == "Datadx":
        print("dx:",f)
        dXs.append(load(f))
    elif f[:6] == "Datady":
        print("dy:",f)
        dYs.append(load(f))


rois = []


ki = 0
for i in range(6):
    wdir = wdirs[i]
    files = os.listdir(wdir)
    files.sort()
    dfiles = []
    for f in files:
        if f[-4:]=='.tif': 
            try:
                idd = int(f[-5])
                dfiles.append(wdir+f)
            except:
                pass

    x1,y1 = 1.0*Xs[i],1.0*Ys[i]
    x2,y2 = 1.0*dXs[i],1.0*dYs[i]
    x2 = x1+x2
    y2 = y1+y2
    print(x1.shape,y1.shape)
    if len(dfiles)!= x1.shape[0]:
        print("Wrong number of files!")
        
    for k in range(len(dfiles)):
        cname = dfiles[k]
        posr = column_stack((x1[k,:500],y1[k,:500]))
        posb = column_stack((x2[k,:500],y2[k,:500]))
        pos = array(row_stack((posr,posb)),dtype=int)
        roit = extractROIS(cname,pos)
        #rois.append(roit)
        with open("rois_random"+str(ki).zfill(2)+".npy","wb") as f:
            pickle.dump(roit,f)
        ki += 1



