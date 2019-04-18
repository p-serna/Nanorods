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

%pylab

def burstsearch(dS,thf=0.6,ncycles = 12):
    ''' It returns a set of pairs of numbers with the index of the beggining 
    and last cycle in a row in a burst interval
    '''
    Ssc = (running_mean(dS,12))**2
    threshold = (max(Ssc)-min(Ssc))*thf+min(Ssc)

    sel = arange(Ssc.shape[0])[Ssc>threshold]
    #sel = sel[(sel>0)*(sel<Ssc.shape[0]-12)]
    csel = sel[1:]-sel[:-1]
    s0 = 0
    ssel = []
    for i,si in enumerate(csel):
        if si>1:
            ssel.append([sel[s0],sel[i]])
            s0 = i+1
    ssel.append([sel[s0],sel[-1]])

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

    return(bselm)
    
def overlapintervals(int0,int1):
    di0 = int0[:,1]-int0[:,0]
    di1 = int1[:,1]-int1[:,0]
    
    total0 = di0.sum()
    total1 = di1.sum()
    if total1 >total0:
        overlap = overlapintervals(int1,int0)
        return(overlap)

    int1l = len(int1)
    overlap = 0
    k = 0
    for s0,s1 in int0:
        while k<int1l:
            t0,t1 = int1[k]
            if t1< s0:
                k +=1
            elif t1<= s1:
                overlap+= t1-max(t0,s0)
                k += 1
                break
            elif t1>s1:
                if t0<s1:
                    overlap += t1-max(t0,s0)
                    
                break
                    
    return( overlap/total0)
                
blues = ["/mnt/data/Anastasia/traces/ROI_990_OnState_SigSig_blue.dat","/mnt/data/Anastasia/traces/ROI_180_OnState_SigSig_blue.dat","/mnt/data/Anastasia/traces/ROI_103_OnState_SigSig_blue.dat","/mnt/data/Anastasia/traces/ROI_990_OnState_ContrlSig_blue.dat","/mnt/data/Anastasia/traces/ROI_180_OnState_ContrlSig_blue.dat","/mnt/data/Anastasia/traces/ROI_103_OnState_ContrlSig_blue.dat"]
reds = ["/mnt/data/Anastasia/traces/ROI_990_OnState_SigSig_red.dat","/mnt/data/Anastasia/traces/ROI_180_OnState_SigSig_red.dat","/mnt/data/Anastasia/traces/ROI_103_OnState_SigSig_red.dat","/mnt/data/Anastasia/traces/ROI_990_OnState_ContrlSig_red.dat","/mnt/data/Anastasia/traces/ROI_180_OnState_ContrlSig_red.dat","/mnt/data/Anastasia/traces/ROI_103_OnState_ContrlSig_red.dat"]
blues.sort()
reds.sort()

thf = 0.6

for red,blue in zip(reds,blues):
    FB = loadtxt(blue)
    FA = loadtxt(red)
    St = FA+FB
    SA = FA[isfinite(FA)]
    SB = FB[isfinite(FB)]


    # ~ Stpurged = 1.0*St
    Stt = St.reshape(St.shape[0]//4,4)
    if mean(SA)>=mean(SB):
        SXt = FA.reshape(FA.shape[0]//4,4)
    else:
        SXt = FB.reshape(FB.shape[0]//4,4)    
        

    SXt = SXt/Stt
    dR = (SXt[:,:2].sum(axis=1)-SXt[:,2:].sum(axis=1))
    dR = dR[isfinite(dR)]
    SXt = SXt[isfinite(SXt.sum(axis=1)),:]


    dS = (Stt[:,:2].sum(axis=1)-Stt[:,2:].sum(axis=1))
    dS = dS[isfinite(dS)]
    Stt = Stt[isfinite(Stt.sum(axis=1)),:]
    # ~ Sttn = Stt[isfinite(Stt)]
    print(Sttn.shape)
    dS = dS[isfinite(dS)]/mean(Stt)

    Ssc = (running_mean(dS,12))**2
    dRsc = (running_mean(dR,12))**2

    #hist(Ssc,bins=31)

    bselm = burstsearch(dS,thf=0.6,ncycles = 12)


    dSsel = []
    Sttsel = []
    dRsel = []
    for bs in bselm:
        dSsel.extend(dS[bs[0]:bs[1]])
        Sttsel.extend(Stt[bs[0]:bs[1],:])
        dRsel.extend(dR[bs[0]:bs[1]])

    dSsel = array(dSsel)
    Sttsel = array(Sttsel)
    dRselF = array(dRsel)
    
    figure(red[:-8]+"F")

    Sscore = []
    S2score = []
    Snscore  = []
    ts0 = 0
    for bs in bselm:
        ts = ts0+arange((bs[1]-bs[0])*4)
        Stts = Stt[bs[0]:bs[1],:]
        Stts = column_stack((Stts[:,:2].mean(axis=1),Stts[:,2:].mean(axis=1)))
        Sscore.append((Stts[:,0]-Stts[:,1]).sum()/Stts.mean())
        S2score.append(((Stts[:,0]-Stts[:,1])**2).sum()/Stts.mean()**2)
        Snscore.append(len(Stts[:,0]))
        Stts = Stts.flatten()
        idx = (ts-ts0)%4
        sel2 = (idx==0)+(idx==2)
        plot(ts[sel2],Stts,'C0.-',alpha=0.5)
        sel = idx[sel2]//2==1
         
        plot(ts[sel2][sel],Stts[sel],'C1.',alpha=0.5)
        vlines(ts0,min(Stts),max(Stts),linestyles='--',alpha=0.5)
        ts0 = ts[-1]+1
    xlabel("t(frames)")
    ylabel("F(cps)")        
    Sscore = array(Sscore)
    S2score = array(S2score)
    Snscore = array(Snscore)

    #figure()
    #plot(arange(len(Stsel))%4+randn(len(Stsel))*0.05,Stsel/mean(Sttn),'.')
    selt = arange(len(Stsel))%4
    
    burstIdF = array(bselm) 

    tempup = Sttsel[:,:2].flatten()/Sttsel.mean()
    temd = Sttsel[:,2:].flatten()/Sttsel.mean()
    
    
    # ~ figure(red[:-8]+"a")
    # ~ title("dF/F")
    # ~ dSsel = array(dSsel)
    # ~ ts = arange(len(dSsel))
    # ~ selu = dSsel>0
    # ~ plot(ts,ts*0+0.0,'k--')
    # ~ plot(ts,dSsel,'.',alpha=0.5)   
    # ~ plot(ts[selu],dSsel[selu],'C1.',alpha=0.5)
    
    
    print("Particle "+red[:-8])
    print("Selected dF/F +/- error")
    print(mean(dSsel),std(dSsel)/sqrt(len(dSsel)))
    print("Total dF/F +/- error")
    print(mean(dS),std(dS)/sqrt(len(dS)))
    
    # ~ print("Selected dF_0/F +/- error, dF_60/F +/- error")
    # ~ print([mean(tempup),std(tempup)/sqrt(len(tempup)),mean(temd),std(temd)/sqrt(len(temd))])
    # ~ print("Selected dF_0/F-dF60/F +/- error")
    # ~ print([mean(tempup)-mean(temd),sqrt(var(tempup)+var(temd))/sqrt(len(tempup)+len(temd))])

    print("Selected (from dF) dR +/- error")
    print(mean(dRsel),std(dRsel)/sqrt(len(dRsel)))

    print("Unselected dR +/- error")
    print(mean(dR),std(dR)/sqrt(len(dR)))

    bselm = burstsearch(dR,thf=0.6,ncycles = 12)
    burstIdR = array(bselm)

    dRsel = []
    Rtsel = []
    for bs in bselm:
        dRsel.extend(dR[bs[0]:bs[1]])
        Rtsel.extend(SXt[bs[0]:bs[1],:])
    Rtt = array(Rtsel)
        
    figure(red[:-8]+"R")
    
    Rscore = []
    R2score = []
    Rnscore = []
    ts0 = 0
    for bs in bselm:
        ts = ts0+arange((bs[1]-bs[0])*4)
        Stts = SXt[bs[0]:bs[1],:]
        Stts = column_stack((Stts[:,:2].mean(axis=1),Stts[:,2:].mean(axis=1)))
        Rscore.append((Stts[:,0]-Stts[:,1]).sum()/Stts.mean())
        R2score.append(((Stts[:,0]-Stts[:,1])**2).sum()/Stts.mean()**2)
        Rnscore.append(len(Stts[:,0]))
        Stts = Stts.flatten()
        idx = (ts-ts0)%4
        sel2 = (idx==0)+(idx==2)
        plot(ts[sel2],Stts,'C0.-',alpha=0.5)
        sel = idx[sel2]//2==1
         
        plot(ts[sel2][sel],Stts[sel],'C1.',alpha=0.5)
        vlines(ts0,min(Stts),max(Stts),linestyles='--',alpha=0.5)
        ts0 = ts[-1]+1
    xlabel("t(frames)")
    ylabel("X/(R+B) (X=R if <R> > <B> or X=B otherwise)")        
    Rscore = array(Rscore)
    R2score = array(R2score)
    Rnscore = array(Rnscore)
    
    # ~ figure(red[:-8]+"b")
    # ~ title("dR")
    # ~ dRsel = array(dRsel)
    # ~ ts = arange(len(dRsel))
    # ~ selu = dRsel>0
    # ~ plot(ts,ts*0+0.0,'k--')
    # ~ plot(ts,dRsel,'C0.',alpha=0.5)
    # ~ plot(ts[selu],dRsel[selu],'C1.',alpha=0.5)
    #plot(ts[seld],dRsel[seld],'C1.',alpha=0.5,label="down")
    #plot(ts[selu],dRsel[selu]-dRsel[seld],'.-',alpha=0.5)

    print("Selected (from dR) dR +/- error")
    print(mean(dRsel),std(dRsel)/sqrt(len(dRsel)))


    
    b = burstIdF*1.0
    db = b[:,1]-b[:,0]
    print("\n Burst statistics dF/F")
    print("Number of intervals = ",b.shape[0])
    print("Av. Length = ",mean(db),"+/-",std(db)/sqrt(db.shape[0]))
    print("Total length = ",db.sum())
    print("Burst score = ",mean(Sscore),"+-",std(Sscore)/sqrt(len(Sscore)))
    st2 = sqrt(S2score/Snscore-(Sscore/Snscore)**2)/sqrt(len(Snscore))
    mt2 = sum((Sscore/Snscore)/st2**2)/(sum(1./st2**2))
    et2 = sqrt(1.0/sum(1./st2**2))
    print("Burst score per time (units of frame)= ",mt2,"+- ",et2)    

    b = burstIdR*1.0
    db = b[:,1]-b[:,0]
    print("\n Burst statistics dR")
    print("Number of intervals = ",b.shape[0])
    print("Av. Length = ",mean(db),"+/-",std(db)/sqrt(db.shape[0]))
    print("Total length = ",db.sum())
    print("Burst score = ",mean(Rscore),"+-",std(Rscore)/sqrt(len(Rscore)))
    st2 = sqrt(R2score/Rnscore-(Rscore/Rnscore)**2)/sqrt(len(Rnscore))
    mt2 = sum((Rscore/Rnscore)/st2**2)/(sum(1./st2**2))
    et2 = sqrt(1.0/sum(1./st2**2))
    print("Burst score per time (units of frame)= ",mt2,"+- ",et2) 
    print("\n The overlap between the 2 burst searchs is %.2f %%" % round(100*overlapintervals(burstIdF,burstIdR),2))
    print("\n")
    print(12*"#"+"\n\n")
