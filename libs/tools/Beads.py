import sys
import os
from numpy import *
from matplotlib.pylab import *
import sys
import os
import cv2
from sub.subs import *
from sub.cor2img import transfpar,coincidROI
#sys.path.append("../")


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

cfile = dfiles[0]

visuallog = False

    
wdir = ''
for fs in cfile.split("/")[:-1]:
    wdir = wdir+fs+'/'

outputdir = cfile.split(".")[0]+'output/'
if not os.path.isdir(outputdir):
    try:
        os.system("mkdir "+outputdir)
    except ValueError:
        print("I cannot create the folder for the output!")
        raise SystemExit(0)
        
movie = readtifFile(cfile)
tduration, mheight,mwidth = movie.shape

im1 = mean(movie,axis=0)

imA = im1[:,:512]
imB = im1[:,512:]

posA = extractpeaks2d(imA,kernel=5)

posB = extractpeaks2d(imB,kernel=5)

height,width = imA.shape
posA = stripborder(posA,height,width)
posB = stripborder(posB,height,width)

savetxt(outputdir+"posA.dat",posA)
savetxt(outputdir+"posB.dat",posA)
if visuallog:
    ion() 
    visualization(imA,posA[:200])
    
    
fposA = posA[:200,:]
fposB = posB[:200,:]

for i,cl in enumerate(fposA):
    roi = selROI(movie,cl[0],cl[1],t=0,size=2,xlim=(0,mwidth//2),ylim=(0,mheight))
    roi = roi.reshape(tduration,5*5)
    
    save(outputdir+"roi_sA"+str(i).zfill(4)+".npy",roi)
for i,cl in enumerate(fposB):
    roi = selROI(movie[:,:,512:],cl[0],cl[1],t=0,size=2,xlim=(0,mwidth//2),ylim=(0,mheight))
    roi = roi.reshape(tduration,5*5)
    
    save(outputdir+"roi_sB"+str(i).zfill(4)+".npy",roi)
