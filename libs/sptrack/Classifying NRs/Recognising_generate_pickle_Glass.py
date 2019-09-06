# coding: utf-8
# # Recognising NRs

import sys
import os
sys.path.append("/export/home1/users/bssn/serna/GitIBENS/Nanorods")
sys.path.append("/export/home1/users/bssn/serna/GitIBENS/Nanorods/sptrack/Classifying NRs")
from numpy import *
from matplotlib.pylab import *
import scipy.optimize as opt
import pickle
from scipy.optimize import minimize,least_squares
from scipy.stats import linregress
from signal_analysis.core import dblgausfit,dblgaussd
from msd import extractmsd,extractDe

# Small definition
def bimod(x):
    m1 = mean(x)
    m2 = std(x)
    m3 = mean((x-m1)**3)/m2**3
    m4 = mean((x-m1)**4)/m2**4
    return((m3**2+1)/m4)
def running_mean(x, N):
    cumsumt = cumsum(insert(x, 0, 0)) 
    return (cumsumt[N:] - cumsumt[:-N]) / float(N)
# Folder dirs where we are doing the analysis
basef = "/mnt/data/Anastasia/"
wdirs = ["/mnt/data/Anastasia/Glass"]

      
# We are sorting them by the presence of a tif file... if not, this have to be changed
# to include those files without a tif present


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

dfiles.pop(1)


# This part is to set which NRs/ROIS are inside a mask. To produce this file the scripts
# in inMask.py has had to be run.


for i,cfile in enumerate(dfiles):
    wdir = ''
    cf2 = cfile.split(".")[0].split("/")
    for fs in cf2[:-1]:
        wdir = wdir+fs+'/'
    mfile = wdir+cf2[-1]+'output/ROIs_Inmask.dat'
    try:
        mt = loadtxt(mfile)
        print(mt.shape)
        if i == 0:
            imsk = mt.flatten()
        else:
            imsk = concatenate((imsk,mt.flatten()))
    except:
        pass

#
# Extracting MSD of each field stored as tif.


dataL = []
nwdirs = []
for i,cfile in enumerate(dfiles):
    wdir = ''
    cf2 = cfile.split(".")[0].split("/")
    for fs in cf2[:-1]:
        wdir = wdir+fs+'/'
    wdir = wdir+cf2[-1]+'output/sptrack/'
    nwdirs.append(wdir)
    
    with open(wdir+'data_msdB.pickle', 'rb') as handle:
        datat = pickle.load(handle)
    dataL = [[d for d in datat if type(d[1])==ndarray]]


    xdst,xdet,xdent = extractmsd(dataL,verbose=True)


    Dest, Dept = extractDe(dataL,xdet,xdst)
    print("Part 1, NR:",i)
    xdst[:,-1] = i
    xdet[:,-1] = i
    xdent[:,-1] = i
    Dest[:,-1] = i
    Dept[:,-1] = i
    if i == 0:
        xds = xdst
        xde = xdet
        xden = xdent
        Des = Dest
        Dep = Dept
    else:
        xds = row_stack((xds,xdst))
        xde = row_stack((xde,xdet))
        xden = row_stack((xden,xdent))
        Des = row_stack((Des,Dest))
        Dep = row_stack((Dep,Dept))
        
    #dataL.append(datat)

# We store things here temporarily

with open('/mnt/data/Anastasia/glass_statsB_temp.pkl','wb') as f:
    pickle.dump([xds,xde,xden,Des,Dep],f)

# Second part: fit to 2 gaussians and extract distance and other properties of the ROI

driftcorrected = True
distmks = [0.125**2,0.25**2,0.5**2,0.75**2,1.0**2-1e-9] 

imsk = []
for i,dirt in enumerate(nwdirs):
    basedir = dirt
    files = os.listdir(basedir)
    dfiles = []
    if driftcorrected:
        for f in files:
            if f[-4:]=='.npy'and f[:4]=='posh' and f[-6:-4]=='DC': dfiles.append(basedir+f)
    else:
        for f in files:
            if f[-4:]=='.npy'and f[:4]=='posh' and f[-6:-4]!='DC': dfiles.append(basedir+f)
    dfiles.sort()
    
    dgpt = []
    idxt = []
    amps = []
    tpas=[]
    for ni,name in enumerate(dfiles):
        # This line is only for glass
        imsk.append(0)
        
        popts = load(name)
        posx = popts[:,5]
        posy = popts[:,6]
        amp = popts[:,0]
        tsr = arange(popts.shape[0])
        th = 0
        th2 = 13.0
        sel = (popts[:,-1]==0)*(posx>0.1)*(posy>0.1)*(posx<4.9)*(posy<4.9)*(amp>exp(th))*(amp<exp(th2))
        Tf = len(posx)
        if sel.sum() > 500:
            amps.append([mean(amp[sel]),mean(popts[sel,4]),bimod(amp[sel])])
            idxt.append(ni)

            m1,m2 = (mean(amp[sel]),std(amp[sel]))
            amb = (amp-m1)/m2
            am = (amp[sel]-m1)/m2
            h = histogram(am,arange(min(am),max(am),0.2))
            hd = h[0]/sum(h[0])/0.2
            hx = (h[1][1:]+h[1][:-1])/2.0
            dgfit = dblgausfit(hx,hd,par0=array([0.5,-1,.5,1,.5]))
            #dgpt.append(concatenate(([m1,m2],dgfit.x)))
            dgt = concatenate(([m1,m2],dgfit.x))

            dgpt.append(dgt)
            
            s1 = max(dgt[6],dgt[4])
            yt = abs((dgt[3]-dgt[5])/s1)
            
            sel2 = sel*(amb>min(dgt[5]+2*dgt[6],dgt[3]+2*dgt[4]))
            if sel2.sum()<500:
                sel2 = sel*(amb>max(dgt[5]-2*dgt[6],dgt[3]-2*dgt[4]))
            if sel2.sum()>500:
                sel = sel2

                
            tsr = tsr[sel]
            x = posx[sel]
            y = posy[sel]
            cm = [mean(x),mean(y)]
            x = x-cm[0]
            y = y-cm[1]
            r = x*x+y*y
            maxr = max(r)
            xm = running_mean(x,10)
            ym = running_mean(y,10)
            maxr2 = max(xm*xm+ym*ym)
            xm = running_mean(x,20)
            ym = running_mean(y,20)
            tm = running_mean(tsr,20)
            maxr3 = max(xm*xm+ym*ym)
            rm = xm*xm+ym*ym
            
            marks = zeros(8)
            marks[:3] = (maxr,maxr2,maxr3)
            
            sel = r>=distmks[0]*maxr
            marks[4] = tsr[sel][0]
            sel = rm>=distmks[0]*maxr3
            marks[5] = tm[sel][0]
            sel = r>=distmks[1]*maxr
            marks[6] = tsr[sel][0]
            sel = rm>=distmks[1]*maxr3
            marks[7] = tm[sel][0]
            
            tpas.append(marks)
            
    print("Part 2, NR:",i)
    #print(i,idx0.shape,len(dfiles))
    if i == 0:
        dgps = column_stack((array(dgpt),zeros(len(dgpt))+i))
        ampst = column_stack((array(amps),zeros(len(dgpt))+i))
        idx0 = array(idxt)
        dtpas = column_stack((array(tpas),zeros(len(tpas))+i))

    else:
        dgps = row_stack((dgps,column_stack((array(dgpt),zeros(len(dgpt))+i))))
        ampst = row_stack((ampst,column_stack((array(amps),zeros(len(dgpt))+i))))
        idx0 = concatenate((idx0,array(idxt)+idx0[-1]+1))
        dtpas = row_stack((dtpas,column_stack((array(tpas),zeros(len(tpas))+i))))


imsk = array(imsk)

# Storing last things in hard drive

with open('/mnt/data/Anastasia/glass_statsB.pkl','wb') as f:
    #pickle.dump([idx0,dgps,ampst,xds,xde,xde,nDes,Dep],f)
    pickle.dump([imsk,nwdirs,dgps,ampst,idx0,xds,xde,xden,Des,Dep,dtpas],f)

# And goodbye
print("Yay! Done!")
