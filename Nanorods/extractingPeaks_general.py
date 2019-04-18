
# Inicialo como
# CUDA_VISIBLE_DEVICES="" ipython3

from numpy import *
from matplotlib.pylab import *
import sys
import os
from keras.models import load_model
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

savetxt(outputdir+"posA.dat",posA)
savetxt(outputdir+"posB.dat",posA)

#--------------------- Visualization   
if visuallog:
    ion() 
    visualization(imA,posA[:200])

# ----------------------------------------------------------------
# Now that we have peaks, we can classify them
model = load_model("/users/bssn/serna/anastasia/ROIS_raw/FullMovies/classifier20181121_17.h5")
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"   # see issue #152
os.environ["CUDA_VISIBLE_DEVICES"] = ""
from tensorflow.python.client import device_lib
print(device_lib.list_local_devices())
from sub.sub_Classifier import extract_data_Class20181121 as extract_data
from sub.sub_Classifier import formatClass20181121 as formatXt

 
ths = zeros(len(posA))
for i,cl in enumerate(posA):
    Xtest = []
    roi = selROI(movie,cl[0],cl[1],t=0,size=1,xlim=(0,mwidth//2),ylim=(0,mheight))
    roi = roi.reshape(tduration,3*3)

    for j in range(2):      
        xtt = extract_data(roi)
        Xtest.append(xtt)
    Xtest = formatXt(Xtest)
    Ytestp = model.predict(Xtest)
    Yres = Ytestp.flatten()
    Yres = mean(Yres)
    ths[i] = Yres

fposA = posA[ths>0.5,:]

if visuallog:
    visualization(imA,fposA[:300],figname="Blue channel")
    
savetxt(outputdir+"FposA.dat",fposA)

ths = zeros(len(posB))
for i,cl in enumerate(posB):
    Xtest = []
    roi = selROI(movie[:,:,512:],cl[0],cl[1],t=0,size=1,xlim=(0,mwidth//2),ylim=(0,mheight))
    roi = roi.reshape(tduration,3*3)
    for j in range(2):      
        xtt = extract_data(roi)
        Xtest.append(xtt)
    Xtest = formatXt(Xtest)
    Ytestp = model.predict(Xtest)
    Yres = Ytestp.flatten()
    Yres = mean(Yres)
    ths[i] = Yres

savetxt(outputdir+"FposB.dat",posB[ths>0.5,:])

fposB = posB[ths>0.5,:]

for i,cl in enumerate(fposA):
    roi = selROI(movie,cl[0],cl[1],t=0,size=2,xlim=(0,mwidth//2),ylim=(0,mheight))
    roi = roi.reshape(tduration,5*5)
    
    save(outputdir+"roi_sA"+str(i).zfill(4)+".npy",roi)
for i,cl in enumerate(fposB):
    roi = selROI(movie[:,:,512:],cl[0],cl[1],t=0,size=2,xlim=(0,mwidth//2),ylim=(0,mheight))
    roi = roi.reshape(tduration,5*5)
    
    save(outputdir+"roi_sB"+str(i).zfill(4)+".npy",roi)



fposAt = posA[:100,:]
fposBt = posB[:100,:]

fsel = coincidROI(fposAt,fposBt,err=10)
err0 = 10
while fsel.shape[0]<50:
    err0 +=1
    fsel = coincidROI(fposAt,fposBt,err=err0)

fposAt2 = fposAt[fsel[:,0],:]
fposBt2 = fposBt[fsel[:,1],:]
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

fposAinB =  array(fposA+v0,dtype=int)
#################### CRITERION FOR SELECTION OF ROIS!!!

selection = coincidROI(fposAinB,posB,err = 2)
savetxt(outputdir+"roi_coincidentROIS.dat",selection)

if visuallog:
    visualization(imB,fposB[:300],figname="Red channel")
    visualization(imB,fposAinB[:300],figname="Red channel",color='green')



if visuallog:
    show()
