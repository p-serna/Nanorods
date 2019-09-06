from numpy import *
from matplotlib.pylab import *
import sys
import os
import pytiff
from sub.subs import *
from sub.cor2img import transfpar,coincidROI

visuallog = False

if len(sys.argv)>1:
    cfile = sys.argv[1] 
else:
    ifile = 0
    wdir = "./"
    print("No argument, we do it for file number 0")


print("We will do analysis for file "+cfile)

    
wdir = ''
for fs in cfile.split("/")[:-1]:
    wdir = wdir+fs+'/'

outputdir = cfile.split(".")[0]+'output/'
        

fposB = loadtxt(outputdir+"FposB.dat")
fposA = loadtxt(outputdir+"FposA.dat")

cellname = cfile.split("/")[-1].split("_")[0]
mask = pytiff.Tiff(wdir+cellname+"_Mask_CMOS.tif")[:,:]

v = transfpar(fposB,fposA,transform = 2)
v0 = v[:2]
v1 = v[2:4]
v2 = array([[0,v[4]],[v[5],0]])
fposBinA =  fposB*v1+v0+dot(v2,fposB.transpose()).transpose()



Av = zeros((fposA.shape[0]+fposB.shape[0]),dtype=int)
ipos = array(fposA,dtype=int)

lA = len(ipos)
for i in range(len(ipos)):
    Av[i] =  mask[ipos[i,1],ipos[i,0]]//255
    

ipos = array(fposBinA,dtype=int)
sel = (ipos[:,0]<mask.shape[1])*(ipos[:,1]<mask.shape[0])*(ipos[:,0]>=0)*(ipos[:,1]>=0)
#ipos = ipos[sel,:]
for i in range(len(ipos)):
    if sel[i]:
        Av[i+lA] =  mask[ipos[i,1],ipos[i,0]]//255
    else:
        Av[i+lA] = 0
        
        
savetxt(outputdir+"ROIs_Inmask.dat",Av)
