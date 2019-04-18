from numpy import *
from matplotlib.pylab import *
import os
import scipy.optimize as opt
import sys
from scipy.optimize import minimize
#import pickle
#from core_libs import *
#import h5py,os
from scipy import signal
#%matplotlib inline
from core import gaussian_fit, get_intervals_fixedD,fixing_intervals2

def running_mean(x, N):
    cumsumt = cumsum(insert(x, 0, 0)) 
    return( (cumsumt[N:] - cumsumt[:-N]) / float(N))

freqadq = 200.0      

if len(sys.argv)>1:
    wdir = sys.argv[1]
    if wdir[-1] != '/':
        wdir = wdir+'/'
else:
    wdir = "./"
    print("No argument, we do it for current folder")

outputdir = wdir+'analysis/'
if not os.path.isdir(outputdir):
    try:
        os.system("mkdir "+outputdir)
    except ValueError:
        print("I cannot create the folder for the output!")
        raise SystemExit(0)

dirt = wdir

basedir = dirt
files = os.listdir(basedir)

dfilesA = []
dfilesB = []

for f in files:
    if f[:6]=='roi_sA' and f[-3:] == 'npy': dfilesA.append(basedir+f)
    if f[:6]=='roi_sB' and f[-3:] == 'npy': dfilesB.append(basedir+f)
 
dfilesA.sort()
dfilesB.sort()

selection = array(loadtxt(basedir+"roi_coincidentROIS.dat"),dtype=int)

dfiles = [[dfilesA[iA],dfilesB[iB]] for iA,iB in selection]

def signal(fname,nmin=3,nmax = None):
    roi = load(fname)
    sh = roi.shape[1]
    roi.sort(axis=1)
    roim = roi[:,:nmin].mean(axis=1)
    if nmax is None:
        nmax = nmin
        S = roi[:,nmin:].sum(axis=1)-roim*(sh-nmin)
    else:
        S = roi[:,-nmax:].sum(axis=1)-roim*nmax
    return(S)
        

for fa,fb in dfiles:
    FA = signal(fa); 
    FB = signal(fb); 
    # Total signal
    St = FA+FB

    namesp = (fa.split(sep=".")[0]).split(sep="/")[-1]

    Sshape = FA.shape

    # ~ namefigure = namesp
    # ~ figure(namefigure,figsize=(16,6))
    t = arange(0,Sshape[0])/freqadq  # in seconds
    # Regularise
    meanS = mean(St)
    stS = std(St)
    Sn = (St-meanS)/stS
    
    hlim = max(max(abs(Sn)),4)
    par = gaussian_fit(Sn,bins=linspace(-hlim,hlim,101))
    # Probably need to change mu2 -> mu1+e
    # par: p, mu1, s1, mu2, s2
    parR = 1.0*par
    parR[1:] = parR[1:]*stS
    parR[array([1,3])] += meanS
    on = 1
    threshold = parR[1]+parR[2]*2.0
    # ~ figure(namefigure,figsize=(16,6))
    # ~ plot(t,St,'C0.',alpha=0.5)
    # ~ plot(t,t*0+threshold,'r--')
    # ~ plot(t[St>threshold],St[St>threshold],'C1.',alpha=0.5)
    
    # ~ #draw(namefigure) 
    # ~ #show(block=False)               
    on = 1#input("Do you want to continue, it is not bimodal: 0, 1")
    # ~ close(namefigure)
        
        
        
    intervals = nan
    if on==1:
        intervals = get_intervals_fixedD(Sn,par,probl = 0.1)
    if intervals is nan: 
        print("nan!")
    else:
       s0,s1 = fixing_intervals2(intervals,length=7)
       if len(s0)<10: break
        #print(min(Sn),max(Sn))
        
        FApurged = []
        FBpurged = []
        Stindics = []
        # ~ s0A = s0
        # ~ s1A = s1
        difs = {k:[] for k in range(-3,8)}
        for i in range(s1.shape[0]):
            # indices of interval?
            indcs = arange(s1[i,0],s1[i,1])
            # times of interval
            tt = t[indcs]
            # their position in the sq wave
            tt0 = indcs%4
            # First time there is a beginning of a sq wave 
            sel0 = arange(len(tt0))[tt0 == 0]
            indc0 = indcs[sel0]        
            sel3 = arange(len(tt0))[tt0 == 3]
            indc3 = indcs[sel3]
            indcs = indcs[sel0[0]:sel3[-1]+1]
            
            #FA purged        
            FApurged.extend(FA[indcs].tolist())
            FBpurged.extend(FA[indcs].tolist())
            Stindics.extend(indcs.tolist())
            
        # ~ if len(FApurged)==0: break
        FApurged = array(FApurged)
        FBpurged = array(FBpurged)
        Stindics = array(Stindics)
        Stpurged = St[Stindics]
        
    Stt = Stpurged.reshape(Stpurged.shape[0]//4,4)
    # ~ Sttn = Stt[isfinite(Stt)]
    # ~ print(Sttn.shape)
    dS = (Stt[:,:2].sum(axis=1)-Stt[:,2:].sum(axis=1))
    # ~ dS = dS[isfinite(dS)]/mean(Sttn)
    Ssc = running_mean((dS)**2,12)     
    threshold = max(Ssc)*.6
    sel = arange(Ssc.shape[0])[Ssc>threshold]
    sel = sel[(sel>0)*(sel<Ssc.shape[0]-12)]
    csel = sel[1:]-sel[:-1]
    s0 = 0
    ssel = []
    for i,si in enumerate(csel):
        if si>1:
            ssel.append([sel[s0],sel[i]])
            s0 = i+1
    
    bsel = []
    for s in ssel:
        bsel.append([s[0],s[1]+12])
    
    bselm = []
    bs0a = bsel[0][0]
    bs1a = bsel[0][1]
    for bs in bsel[1:]:
        bs0 = bs[0]
        bs1 = bs[1]
        if bs1a< bs0:
            bselm.append([bs0a,bs1a])
            bs0a = bs0
        bs1a = bs1
    bselm.append([bs0a,bs1a])
    
    dSsel = []
    Stsel = []
    for bs in bselm:
        dSsel.extend(dS[bs[0]:bs[1]])
        Stsel.extend(Stt[(bs[0]*4):(bs[1]*4)])



