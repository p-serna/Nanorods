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
    St = FA#+FB

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
        # ~ intpre = array(intervals)
        #~ print("Intervals shape=",intpre.shape)
        s0,s1 = fixing_intervals2(intervals,length=7)
        #~ print("And now shape=",s1.shape)
        #~ s0,s1 = fixing_intervals(s0,s1,length=7)
        #~ print("And now! shape=",s1.shape)
        #~ s0,s1 = fixing_intervals(s0,s1,length=7)
        #~ print("And now!2 shape=",s1.shape)
        if len(s0)<10: break
        #print(min(Sn),max(Sn))
        
        # what is this for? I think it can be removed
        # ~ xbins = linspace(min(Sn)*0.99,max(Sn)*1.01,51)
        # ~ xbinb = (xbins[1:]+xbins[:-1])/2.0
        # ~ h0 = histogram(Sn,xbins)[0]
        

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
            
            #Need reviewing!!!!            
            # ~ for k in [0,1]:
                # ~ sel1 = sel0+3
                # ~ sel1 = sel1[sel1<len(indcs)]
                # ~ indc1 = indcs[sel1]
                # ~ for k1,k2 in zip(indc0,indc1):
                    # ~ Sn2 = (St[k2]+St[k2-1-k])/2.0
                    # ~ Sn1 = (St[k1]+St[k1+1+k])/2.0
                    # ~ difs[k].append((Sn1-Sn2)/meanS)
                    # ~ Sn2 = (Stb[k2]+Stb[k2-1-k])/2.0
                    # ~ Sn1 = (Stb[k1]+Stb[k1+1+k])/2.0
                    # ~ difs[k+2].append(Sn1-Sn2)
                    # ~ Sn2 = (St[k2]/(St[k2]+Stb[k2])+St[k2-1-k]/(St[k2-1-k]+Stb[k2-1-k]))/2.0
                    # ~ Sn1 = (St[k1]/(St[k1]+Stb[k1])+Stb[k1+1+k]/(St[k1+1+k]+Stb[k2-1-k]))/2.0
                    # ~ difs[k+4].append(Sn1-Sn2)
                    
            #FA purged        
            FApurged.extend(FA[indcs].tolist())
            FBpurged.extend(FA[indcs].tolist())
            Stindics.extend(indcs.tolist())
            
        if len(stnpurged)==0: break
        stnpurged = array(stnpurged)
        stnindics = array(stnindics)
        
        # ~ namefigt = namefigure+"Tseries"
        # ~ figure(namefigt,figsize=(16,6))
        # ~ plot(t,Sn,'.',alpha=0.3)
        # ~ plot(t[stnindics],stnpurged,'.',alpha=0.6)
        # ~ xlabel("t/s")
        # ~ ylabel("F (arbitrary units)")
        # ~ savefig(namefigt+".png")
        # ~ close(namefigt)
        
        # I LEFT IT HERE, BUT WE NEED TO REVIEW PREVIOUS SECTION!!!!
        
        
        
        
        
        on = False
        on2 = False
        xlabels = ["$(F(0) - F(\Delta t))/\\langle F\\rangle$","$(F_b(0)-F_r(0))- (F_b(\Delta t)-F_r(\Delta t)) $ (fluorescence)","$(F_b(0)-F_r(0))/F(0)- (F_b(\Delta t)-F_r(\Delta t))/F(\Delta t) $"]
        for kt in range(3):
            namefigt = namefigure+"Disdiff"+str(kt)
            figure(namefigt)
            #meanS
            dmin = min(array([min(difs[i]) for i in range(2*kt,2*(kt+1))]))
            dmax = max(array([max(difs[i]) for i in range(2*kt,2*(kt+1))]))
            xbins = linspace(dmin*0.99,dmax*1.01,21)
            xbinb = (xbins[:-1]+xbins[1:])/2.0
            ylimt = -1
            xm = []
            for i in range(2*kt,2*(kt+1)):
                ht = histogram(difs[i],xbins)[0]
                plot(xbinb,ht,'.-',label=5*(2-i%2))
                ylimt = max(ylimt,max(ht))
                xm.append([mean(difs[i]),sqrt(var(difs[i])),len(difs[i])])
                #vlines(mean(difs[i]),0,1,'C'+str(i))
            ylimt = ylimt
            xlimt = min(xbins)
            x1,x2,x3 = (xm[0][0],xm[0][1]/sqrt(xm[0][2]),xm[0][1])
            x1b,x2b,x3b = (xm[1][0],xm[1][1]/sqrt(xm[1][2]),xm[1][1])
            text_l = '''mean$\pm$ errmean, sd = 
                                Blue)\t %.2e$\pm$ %.2e,\t %.2e
                                Orange)\t %.2e$\pm$ %.2e,\t %.2e
                        ''' %(x1,x2,x3,x1b,x2b,x3b)
            for ir in range(1,6):
                if abs(x1)/abs(x2)>ir : 
                    text_l = '$\\bf{*}$'+text_l
                    on = True
                if abs(x1b)/abs(x2b)>ir : 
                    text_l = text_l+'$\\bf{*}$'
                    on2 = True
            text_l = text_l+"\n"

                                        
            text(xlimt,ylimt*0.95,text_l)
            
            ylim(-.03,ylimt*1.03)
            legend(title = "$\Delta t$/ms")
            xlabel(xlabels[kt])
            ylabel("Distribution")
            savefig(namefigt+".png")
            close(namefigt)
            if on:
                if abs(x1)/abs(x2)>4 : 
                    relevance4.append(ish)
                elif abs(x1)/abs(x2)>3 : 
                    relevance3.append(ish)
                elif abs(x1)/abs(x2)>2 : 
                    relevance2.append(ish)
                elif abs(x1)/abs(x2)>1 : 
                    relevance1.append(ish)
            if on2:
                if abs(x1b)/abs(x2b)>4 : 
                    relevance4b.append(ish)
                elif abs(x1b)/abs(x2b)>3 : 
                    relevance3b.append(ish)
                elif abs(x1b)/abs(x2b)>2 : 
                    relevance2b.append(ish)
                elif abs(x1b)/abs(x2b)>1 : 
                    relevance1b.append(ish)                    
        
        for i in range(3): close()
        #print(-1./ps0[0],-1./ps1[0])
    print("File name:",n)
    print("Relevance 4 sigma: ",relevance4)
    print("Relevance 3 sigma: ",relevance3)
    print("Relevance 2 sigma: ",relevance2)
    print("Relevance 1 sigma: ",relevance1)
    print("Control 3 sigma: ",relevance3b)
    print("Control 2 sigma: ",relevance2b)
    print("Control 1 sigma: ",relevance1b)

