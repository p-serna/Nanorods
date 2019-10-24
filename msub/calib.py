import os
import numpy as np
from Nanorods.sub.cor2img import transfpar,coincidROI
from matlab.tifmethods import readtifImage

def extractpeaks2d0(img,ROIsize = 5, kernel = 5):
    km = np.ones((kernel,kernel),np.float32)
    J = np.zeros(img.shape)
    rs2 = kernel//2
    mx,my = img.shape

    for x in range(0,mx-kernel):
        for y in range(0,my-kernel):
            temp = km*img[x:(x+kernel),y:(y+kernel)]
            J[x+rs2,y+rs2] = temp.max()

    its = []
    pts = [] 
    for x in range(2*ROIsize+1,mx-2*ROIsize):
        for y in range(2*ROIsize,my-2*ROIsize):
            if J[x, y]==img[x, y]:
                pts.append([y, x])
                its.append(img[x, y])
    
    sel = (-np.array(its)).argsort()
    #print(array(its)[sel][:4])
    return(np.array(pts)[sel])


def afintransf(posA,posB,err0 = 10,nbeads = 20):

    fsel = coincidROI(posA,posB,err=err0)
    while fsel.shape[0]<nbeads:
        err0 +=1
        fsel = coincidROI(posA,posB,err=err0)

    posAt1 = posA[fsel[:,0],:]
    posBt1 = posB[fsel[:,1],:]
    v0 = posBt1.mean(axis=0)-posAt1.mean(axis=0)
    # We add noise to refit!?
    v0 = v0+np.random.randn(2)
    v = transfpar(posAt1,posBt1,transform = -1,x0 = v0)
    v0 = v[:2]

    posAt1 = np.array(posAt1+v,dtype=int)

    newsh = fsel.shape[0]
    # By pixel
    err1 = 1
    fsel = coincidROI(posAt1,posBt1,err=err1)
    while fsel.shape[0]<newsh//2:
        err1 +=.5
        fsel = coincidROI(posAt1,posBt1,err=err1)

    posAt2 = posAt1[fsel[:,0],:]
    posBt2 = posBt1[fsel[:,1],:]
    v0n = posBt2.mean(axis=0)-posBt2.mean(axis=0)
    v = transfpar(posAt2,posBt2,transform = -1,x0 = v0)
    v1 = v[:2]
    #print(fsel.shape)
    return(v0+v1)


def CMOS_2fields_calib(folder='./',nbeads = 20):
    files_tif = []
    for file in os.listdir(folder):
        if file.endswith('.tif'):
            files_tif.append(file)
            
    files_tif.sort()
    
    files = files_tif
    filesCMOS = [file for file in files if file.find('CMOS')>=0]
    if len(filesCMOS) == 0:
        raise ValueError('Sorry, there are no files with CMOS in their name.')
    #print(filesCMOS)
     
    vs = []
    for i,fCMOS in enumerate(filesCMOS): 
        cmosim = readtifImage(folder+fCMOS)
        sY, sX = cmosim.shape
        sX2 = sX//2
        try:        
            ptsB = extractpeaks2d0(cmosim[:,:sX2])
            ptsR = extractpeaks2d0(cmosim[:,sX2:])

            v = afintransf(ptsB[:200,:],ptsR[:200,:],nbeads = 100)
            vs.append(v)
        except Exception as e: 
            print('Problem in file '+fCMOS+ ': ',e)
            
    vs = np.array(vs)
    return(vs.mean(axis=0))



