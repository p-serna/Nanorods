
# Inicialo como
# CUDA_VISIBLE_DEVICES="" ipython3

from numpy import *
from matplotlib.pylab import *
import sys
import os
sys.path.append("/users/bssn/serna/GitIBENS/Nanorods")
#from keras.models import load_model
import cv2
from sub.subs import *
from sub.cor2img import transfpar,coincidROI


visuallog = False

if len(sys.argv)>1:
    cfile = sys.argv[1] 
else:
    cfile = 0
    wdir = "./"
    print("No argument, we do it for file 0")


print("We will do analysis for file "+cfile)

    
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



# ----------------------------------------------------------------
# We read the file.     

movie = readBigtifFile(cfile)
tduration, mheight,mwidth = movie.shape

# ----------------------------------------------------------------
# We could correct the counts per frame with times     

# ~ times = gettimes(f)


#-----------------------------------------------------------------
# This part is to identify the local maxima

im1 = mean(movie,axis=0)

imA = im1[:,:512]
imB = im1[:,512:]

posA = extractpeaks2d(imA,kernel=5)

posB = extractpeaks2d(imB,kernel=5)

height,width = imA.shape
posA = stripborder(posA,height,width)
posB = stripborder(posB,height,width)

#savetxt(outputdir+"posA.dat",posA)
#savetxt(outputdir+"posB.dat",posA)

fposA = posA[:100,:]
fposB = posB[:100,:]

fsel = coincidROI(fposA,fposB,err=10)
err0 = 10
while fsel.shape[0]<50:
    err0 +=1
    fsel = coincidROI(fposA,fposB,err=err0)

fposAt2 = fposA[fsel[:,0],:]
fposBt2 = fposB[fsel[:,1],:]
v0 = fposBt2.mean(axis=0)-fposAt2.mean(axis=0)
v0 = v0+randn(2)

v = transfpar(fposAt2,fposBt2,transform = -1,x0 = v0)
v0 = v[:2]

fposAt2 = array(fposAt2+v,dtype=int)

newsh = fsel.shape[0]
err0 = 1
fsel = coincidROI(fposAt2,fposBt2,err=1)
while fsel.shape[0]<newsh//2:
    err0 +=.5
    fsel = coincidROI(fposAt2,fposBt2,err=err0)

fposAt3 = fposAt2[fsel[:,0],:]
fposBt3 = fposBt2[fsel[:,1],:]
v0n = fposBt3.mean(axis=0)-fposAt3.mean(axis=0)
v = transfpar(fposAt3,fposBt3,transform = -1,x0 = v0)
v1 = v[:2]
fposAt3 = array(fposAt3+v1,dtype=int)

v0 = v0+v1



visuallog = True
if visuallog:
    visualization(imB,fposBt3,figname="Red channel")
    visualization(imB,fposAt3,figname="Red channel",color='green')
    outputdir = "/mnt/data/Anastasia/Experimental_data/fieldsred/"

    files = os.listdir(outputdir)
    ifilen = str(len(files)).zfill(2)
    savefig(outputdir+"translationinredchannel"+ifilen+".png")
    
print("v translation:",v0)

sel = permutation(min(posA.shape[0],2000))[:500]
fposA = posA[sel,:]
#fposAinB =  fposA*v1+v0+dot(v2,fposA.transpose()).transpose()
fposAinB =  fposA+v0
fposB = array(fposAinB,dtype=int)

rois = zeros((tduration,25,1000))

for i,cl in enumerate(fposA):
    roi = selROI(movie,cl[0],cl[1],t=0,size=2,xlim=(0,mwidth//2),ylim=(0,mheight))
    rois[:,:,i] = roi.reshape(tduration,5*5)
    
for i,cl in enumerate(fposB):
    roi = selROI(movie[:,:,512:],cl[0],cl[1],t=0,size=2,xlim=(0,mwidth//2),ylim=(0,mheight))
    rois[:,:,i+500] = roi.reshape(tduration,5*5)
    
outputdir = "/mnt/data/Anastasia/Experimental_data/randompeaks2/"
files = os.listdir(outputdir)
save(outputdir+"rois_random"+ifilen+".npy",rois)

outputdir = "/mnt/data/Anastasia/Experimental_data/randompeaks2/posA/"
savetxt(outputdir+"posA"+ifilen+".dat",fposA)
outputdir = "/mnt/data/Anastasia/Experimental_data/randompeaks2/posB/"
savetxt(outputdir+"posB"+ifilen+".dat",fposA)
