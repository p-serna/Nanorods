import sys
import os
import pickle
from numpy import *
from matplotlib.pylab import *
import os


def Sarlebimodal(x):
    m1 = mean(x)
    m2 = var(x)
    m3 = mean((x-m1)**3)/m2**(3.0/2.0)
    m4 = mean((x-m1)**4)/m2**(2.0)-3
    n = len(x)
    fincor = 3*(n-1)**2/((n-2)*(n-3))
    return((m3**2+1)/(m4+fincor))
    
    
def extractdrift(dfiles,ulim = 14, llim = 0.0,bimodalth = 0.5 ,weighted=False):
    popts = load(dfiles[0])
    t = arange(popts.shape[0])
    drift = zeros((popts.shape[0],5))
    driftw = zeros((popts.shape[0],4))
    hdrift = zeros(5)
    if weighted:
        hdrift = zeros(8)
    for ni,name in enumerate(dfiles):
        popts = load(name)
        posx = popts[:,5]
        posy = popts[:,6]
        amp = popts[:,0]
        th = 0
        th2 = ulim # 50e3 
        sel = (popts[:,-1]==0)*(posx>0.1)*(posy>0.1)*(posx<4.9)*(posy<4.9)*(amp>exp(th))*(amp<exp(th2))
        if sel.sum()>1000 and max(amp)<exp(th2):
            y = popts[:,0]/sqrt(popts[:,4])
            if Sarlebimodal(y[sel])>bimodalth:
                y = (y-mean(y[sel]))/std(y[sel])
                sel2 = y>llim
                sel = sel*sel2
            else:
                sel2 = y>20.0
                sel = sel*sel2
            if sum(sel)>1000:
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
                
                if weighted:
                    wx = std(posx[1:]-posx[:-1])/sqrt(len(posx)-1) 
                    wy = std(posy[1:]-posy[:-1])/sqrt(len(posy)-1)
                    hdrift[0] += mean(posx[1:]-posx[:-1])/wx**2 
                    hdrift[1] += mean(posy[1:]-posy[:-1])/wy**2
                    hdrift[2] += 1.0/wx**2
                    hdrift[3] += 1.0/wy**2
                    hdrift[4] +=  mean((posx[1:]-posx[:-1]))**2/wx**2
                    hdrift[5] +=  mean((posy[1:]-posy[:-1]))**2/wy**2
                    hdrift[6] +=  1.0/wx**4
                    hdrift[7] +=  1.0/wy**4
                else:
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

    if weighted:
        hdrift[0] = hdrift[0]/hdrift[2]
        hdrift[1] = hdrift[1]/hdrift[3]
        hdrift[6] = 1.0-hdrift[6]/hdrift[2]**2
        hdrift[7] = 1.0-hdrift[7]/hdrift[3]**2
        hdrift[4] = (hdrift[4]/hdrift[2]-hdrift[0]**2)/hdrift[6]
        hdrift[5] = (hdrift[5]/hdrift[3]-hdrift[1]**2)/hdrift[7]        
        hdrift[2] = sqrt(hdrift[4])
        hdrift[3] = sqrt(hdrift[5])
    else:
        hdrift[0] = hdrift[0]/hdrift[-1]
        hdrift[1] = hdrift[1]/hdrift[-1]
        hdrift[2] = sqrt(hdrift[2]/hdrift[-1]-hdrift[0]**2)/sqrt(hdrift[-1])
        hdrift[3] = sqrt(hdrift[3]/hdrift[-1]-hdrift[1]**2)/sqrt(hdrift[-1])
    
    if weighted:
        return(drift,hdrift[:4])
        
    return(drift,hdrift)
    

def driftSingle(popts,th =0.0,th2 = 14.0,minlength=1000):
    posx = popts[:,5]
    posy = popts[:,6]
    amp = popts[:,0]
    ampbg = popts[:,0]/popts[:,4]
     

    sel = (popts[:,-1]==0)*(posx>0.1)*(posy>0.1)*(posx<4.9)*(posy<4.9)*(amp>exp(th))*(amp<exp(th2))
    if sel.sum()>minlength:
        dx = posx[sel][1:]-posx[sel][:-1]
        dy = posy[sel][1:]-posy[sel][:-1]
        return([dx.mean(),dy.mean(),amp[sel].mean(),dx.std(),dy.std(),ampbg[sel].mean()])
    else:
        return([nan,nan,nan,nan,nan,nan])
        
def extractdrift_clean(dfiles,ulim = 14, llim = 0.0, pxth = 0.1 ,weighted=False):
    popts = load(dfiles[0])
    t = arange(popts.shape[0])
    drift = zeros((popts.shape[0],5))
    driftw = zeros((popts.shape[0],4))
    hdrift = zeros(5)
    if weighted:
        hdrift = zeros(8)
    driftv = zeros((len(dfiles),6))
    for ni,name in enumerate(dfiles):
        popts = load(name)
        driftv[ni,:] = driftSingle(popts,llim,ulim)
        
    xt = driftv*1.0
    xt = xt[isfinite(xt[:,1]),:]
    sx,sy = xt[:,:2].std(axis=0)*2
    mx = xt.mean(axis=0)
    xt = xt[(abs(xt[:,0]-mx[0])<sx)*(abs(xt[:,1]-mx[1])<sy),:]
    xt = xt[(xt[:,3]<pxth)*(xt[:,4]<pxth),:]
    hdrift = concatenate((xt.mean(axis=0)[:2],xt.std(axis=0)[:2]/sqrt(xt.shape[0])))
    return(driftv,hdrift)


def main_estimate(cfile):
    wdir = cfile.split('.')[0]+'output/sptrack/'
    dirt = wdir
    basedir = dirt
    files = os.listdir(basedir)

    dfiles = []
    for f in files:
        if f[-4:]=='.npy'and f[:4]=='posh' and f[-6:-4]!='DC': dfiles.append(basedir+f)

    dfiles.sort()
    nfiles = dfiles 
    drift, hdrift = extractdrift_clean(nfiles,ulim = 16, llim = 0.0 ,weighted=True)
                
    time = 60.0
    dd = hdrift[:4]*325*6000/time

    print(wdir.split("/")[-3],' in ',wdir.split("/")[-4],'has a drift of (%.3f +- %.3f,%.3f +- %.3f) nm/s \n' % (dd[0],dd[2],dd[1],dd[3]))

    return(hdrift)


def main():
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

    #drift, hdrift = extractdrift(dfiles,ulim = 14, llim = 0.0 ,weighted=True)

    # ~ ys = []
    # ~ nfiles = []
    # ~ for fl in dfiles:
        # ~ popts = load(fl)
        # ~ posx = popts[:,5]
        # ~ posy = popts[:,6]
        # ~ amp = popts[:,0]
        # ~ th = 0
        # ~ th2 = 14 # 50e3     
        # ~ sel = (popts[:,-1]==0)*(posx>0.1)*(posy>0.1)*(posx<4.9)*(posy<4.9)*(amp>exp(th))*(amp<exp(th2))
        # ~ if sel.sum()>1000 and max(amp)<exp(th2):
            # ~ y = mean(amp[sel])/mean(popts[sel,4])
            # ~ ys.append(y)
            # ~ if y>2.5:
                # ~ nfiles.append(fl)
        # ~ else:
            # ~ ys.append(0)
    nfiles = dfiles 
    #drift, hdrift = extractdrift(nfiles,ulim = 14, llim = 0.0 ,weighted=True)
    drift, hdrift = extractdrift_clean(nfiles,ulim = 14, llim = 0.0 ,weighted=True)

                
    #nframes = drift.shape[0]
    time = 60.0
    dd = hdrift[:4]*325*6000/time

    print(wdir.split("/")[-3],' in ',wdir.split("/")[-4],'has a drift of (%.3f +- %.3f,%.3f +- %.3f) nm/s \n' % (dd[0],dd[2],dd[1],dd[3]))

    #if plot:
        #t = arange(drift.shape[0])
        #plot(0+t*hdrift[0]-drift[:,0],'.-',alpha=0.5) 
        ## ~ plot(0+t*hdrift[0]-driftw[:,0],'.-',alpha=0.5)

        #plot(0+t*hdrift[1]-drift[:,1],'.-',alpha=0.5) 
        ## ~ plot(0+t*hdrift[1]-driftw[:,1],'.-',alpha=0.5)

        #figure()
        #plot(drift[range(0,drift.shape[0],1),0],drift[range(0,drift.shape[0],1),1],'C1.-',alpha=0.5)
        #plot(drift[range(0,drift.shape[0],100),0],drift[range(0,drift.shape[0],100),1],'C0.-',alpha=0.5)
        ##plot(driftw[range(0,6000,100),0],driftw[range(0,6000,100),1],'.-',alpha=0.5)

        #plot(0+t*hdrift[0],0+t*hdrift[1],'k--',alpha=0.5)

    for ni,name in enumerate(dfiles):
        popts = load(name)
        t = arange(popts.shape[0])
        popts[:,5] = popts[:,5]-t*hdrift[0]
        popts[:,6] = popts[:,6]-t*hdrift[1]

        nname = name.split(".")
        nname = nname[0]+"_DC."+nname[1]
        
        save(nname,popts)

if __name__ == "__main__":
    main()
