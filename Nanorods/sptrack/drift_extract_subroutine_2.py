import sys
import os
import pickle
from numpy import *
from matplotlib.pylab import *
import os
import pickle

def extractdrift(dfiles,ulim = 10.8, llim = 0 ,weighted=False):
    popts = load(dfiles[0])
    t = arange(popts.shape[0])
    drift = zeros((popts.shape[0],5))
    driftw = zeros((popts.shape[0],4))
    hdrift = zeros(5)
    for ni,name in enumerate(dfiles):
        popts = load(name)
        posx = popts[:,5]
        posy = popts[:,6]
        amp = popts[:,0]
        th = llim
        th2 = ulim # 50e3 
        sel = (popts[:,-1]==0)*(posx>0.1)*(posy>0.1)*(posx<4.9)*(posy<4.9)*(amp>exp(th))*(amp<exp(th2))
        if sel.sum()>1000 and max(amp)<exp(th2):
            posx = posx[sel]-posx[0]
            posy = posy[sel]-posy[0] 
            drift[sel,0] += posx 
            drift[sel,1] += posy
            drift[sel,2] += posx**2
            drift[sel,3] += posy**2
            drift[sel,-1] += 1
            
            driftw[sel,0] += posx/popts[sel,1]**2
            driftw[sel,1] += posy/popts[sel,2]**2
            driftw[sel,-2] += 1/popts[sel,1]**2
            driftw[sel,-1] += 1/popts[sel,2]**2
            
            hdrift[0] += mean(posx[1:]-posx[:-1]) 
            hdrift[1] += mean(posy[1:]-posy[:-1])
            hdrift[2] += mean(posx[1:]-posx[:-1])**2 
            hdrift[3] += mean(posy[1:]-posy[:-1])**2
            hdrift[-1] += 1


    drift[:,0] = drift[:,0]/drift[:,-1]
    drift[:,1] = drift[:,1]/drift[:,-1]
    drift[:,2] = sqrt((drift[:,2]/drift[:,-1]-drift[:,0]**2)/drift[:,-1])
    drift[:,3] = sqrt((drift[:,3]/drift[:,-1]-drift[:,1]**2)/drift[:,-1])
    
    drift[0,:] = 0

    driftw[:,0] = driftw[:,0]/driftw[:,-2]
    driftw[:,1] = driftw[:,1]/driftw[:,-1]

    driftw[0,:] = 0

    hdrift[0] = hdrift[0]/hdrift[-1]
    hdrift[1] = hdrift[1]/hdrift[-1]
    hdrift[2] = sqrt(hdrift[2]/hdrift[-1]-hdrift[0]**2)/sqrt(hdrift[-1])
    hdrift[3] = sqrt(hdrift[3]/hdrift[-1]-hdrift[1]**2)/sqrt(hdrift[-1])
    
    if weighted:
        return(driftw,hdrift)
        
    return(drift,hdrift)
    



wdirs = ["/mnt/data/Anastasia/18_11_29_pd23_11_div6_25Hzsqwave/"]

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
wdir = ''
cf2 = cfile.split(".")[0].split("/")
for fs in cf2[:-1]:
    wdir = wdir+fs+'/'
wdir = wdir+cf2[-1]+'output/sptrack/'

wdir+'data_msd.pickle'


if len(sys.argv)>1:
    wdir = sys.argv[1]
    if wdir[-1] != '/':
        wdir = wdir+'/'
    if len(sys.argv)>2:
        interactive = sys.argv[2]=='1'
        plot = True
    else:
        interactive = False
        plot = False

else:
    wdir = "./"
    interactive = False
    print("No argument, we do it for current folder")


dirt = wdir
basedir = dirt
files = os.listdir(basedir)

dfiles = []
for f in files:
    if f[-4:]=='.npy'and f[:4]=='posh' and f[-6:-4]!='DC': dfiles.append(basedir+f)

dfiles.sort()

drift, hdrift = extractdrift(dfiles,ulim = 10.8, llim = 0 ,weighted=False)

nframes = drift.shape[0]
time = 60.0
dd = hdrift[:4]*325*nframes/time

print(wdir.split("/")[-3],' in ',wdir.split("/")[-4],'has a drift of (%.3f +- %.3f,%.3f +- %.3f) nm/s \n' % (dd[0],dd[2],dd[1],dd[3]))

if plot:
    plot(0+t*hdrift[0]-drift[:,0],'.-',alpha=0.5) 
    # ~ plot(0+t*hdrift[0]-driftw[:,0],'.-',alpha=0.5)

    plot(0+t*hdrift[1]-drift[:,1],'.-',alpha=0.5) 
    # ~ plot(0+t*hdrift[1]-driftw[:,1],'.-',alpha=0.5)

    figure()
    plot(drift[range(0,6000,1),0],drift[range(0,6000,1),1],'C1.-',alpha=0.5)
    plot(drift[range(0,6000,100),0],drift[range(0,6000,100),1],'C0.-',alpha=0.5)
    #plot(driftw[range(0,6000,100),0],driftw[range(0,6000,100),1],'.-',alpha=0.5)

    plot(0+t*hdrift[0],0+t*hdrift[1],'k--',alpha=0.5)

for ni,name in enumerate(dfiles):
    popts = load(name)
    t = arange(popts.shape[0])
    popts[:,5] = popts[:,5]-t*hdrift[0]
    popts[:,6] = popts[:,6]-t*hdrift[1]

    nname = name.split(".")
    nname = nname[0]+"_DC."+nname[1]
    
    save(nname,popts)
