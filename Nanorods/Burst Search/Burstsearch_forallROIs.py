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
import pickle

#%pylab
def running_mean(x, N, padding = "valid"):
    if padding =="same":
        N2 = N//2
        Nf = float(N)
        cumsumt = cumsum(concatenate((zeros(1),x,zeros(N-1))))
        runmean = (cumsumt[N:] - cumsumt[:-N]) / Nf
        runmean[-N+1:] = runmean[-N+1:]*Nf/(arange(Nf-1,0,-1))
    elif padding =="valid":
        cumsumt = cumsum(insert(x, 0, 0))
        runmean = (cumsumt[N:] - cumsumt[:-N]) / float(N)
    return(runmean)
    
def burstsearch(dS,thf=0.6,ncycles = 12, merged = True,returnstat = False):
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
    if merged:
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
    else:
        bselm = bsel
        
    if returnstat:
        Ssc.sort()
        Smedian = Ssc[len(Ssc)//2]
        return(bselm,(Ssc.max(),Ssc.mean(),Smedian))
    else:
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
                    
    return( overlap/total1, total0-total1)

def burstscore(b,dS, Fav = 1, idx= None, plot = False):
    Sscore = zeros((b.shape[0],4))
    if idx is None:
        idx  = arange(dS.shape[0])*2
    for i,bs in enumerate(b):
        # ~ print(bs[0],bs[1],dS.shape,idx.shape)
        dSt = dS[bs[0]:bs[1]] 
        t0,tf = idx[bs[0]:bs[1]][array([0,-1])]
        tf += 2
        # Burst score, S= sum dF/F, and burst score per time 
        # (S2 and Sn are there to calculate error bars later for S/time) 
        Sscore[i,:] = ((dSt).sum()/Fav,(dSt*dSt).sum()/Fav**2, dSt.shape[0],tf-t0)
        
    return(Sscore)

def preparedata(FA,FB,framespercycle=4,nframes=6000,removeNans = False):
    framespercycle2 = framespercycle//2
    St = FA+FB
    Fav = St[isfinite(St)].mean()
    # Reshape of total trace in 4 columns, its shape is now  (ncycles,4) 
    Stt = St.reshape(St.shape[0]//framespercycle,framespercycle)
    
    # ratio R of one channel over F. We select the channel with larger average
    if mean(FA[isfinite(FA)])>=mean(FB[isfinite(FB)]):
        SXt = FA.reshape(FA.shape[0]//framespercycle,framespercycle)
    else:
        SXt = FB.reshape(FB.shape[0]//framespercycle,framespercycle)    
    
    SXt = SXt/Stt
    SX = 1.0*SXt
    # Definition of dR and dF, for each interval: R_0+R_1-R_2-R_3   
    dRa = (SXt[:,:framespercycle2].mean(axis=1)-SXt[:,framespercycle2:].mean(axis=1))
    dFa = (Stt[:,:framespercycle2].mean(axis=1)-Stt[:,framespercycle2:].mean(axis=1))
    dRb = -(SXt[:-1,framespercycle2:].mean(axis=1)-SXt[1:,:framespercycle2].mean(axis=1))
    dFb = -(Stt[:-1,framespercycle2:].mean(axis=1)-Stt[1:,:framespercycle2].mean(axis=1))

    dR = concatenate((array(list(zip(dRa,dRb))).flatten(),dRa[-1:]))
    dF = concatenate((array(list(zip(dFa,dFb))).flatten(),dFa[-1:]))
    dF = dF/Fav

    # We remove intervals where there was at least a frame classified as blinking (dark state)    
    idxF = arange(0,nframes-2,2)
    sel = isfinite(dF)
    idxF = arange(0,nframes-2,2)[sel]
    if removeNans:
        sel = isfinite(dR)
        idxR = arange(0,nframes-2,2)[sel]
        dR = dR[sel]
        SXt = SXt[isfinite(SXt.sum(axis=1)),:]
        dF = dF[sel]
        Stt = Stt[isfinite(Stt.sum(axis=1)),:]
        # dF will be defined as dF/<F> with <F> average of the whole trace.
            
    return(St,SX,Stt,SXt, dR,dF,idxF,Fav)


def burstanalysis(FB,FA, verbose = True, plot = False,pname = "",thf = 0.6, ncycles = 12,framespercycle = 4 ):    
    # Printing statistics
    stats = {}

    nframes = FA.shape[0]
    # Total trace, F
    # ~ St = FA+FB
    # ~ Fav = St[isfinite(St)].mean()
    # ~ SA = FA[isfinite(FA)]
    # ~ SB = FB[isfinite(FB)]

    _,_,_,_,dR,dF,idxF,Fav = preparedata(FA,FB,framespercycle,nframes)

    # Burst search : it provides merged intervals [ [idx00,idxf0], [idx01,idxf1], ...]
    bselm, statsF_Brstsrch = burstsearch(dF,thf,ncycles,returnstat=True)
    burstIdF = array(bselm) 
    stats["burstIdF"] = burstIdF
    stats["burstIdF stats"] = statsF_Brstsrch

    # Now we take the parts of the trace selected by the burst search
    dFsel = []
    dRsel = []
    for bs in bselm:
        dFsel.extend(dF[bs[0]:bs[1]])
        dRsel.extend(dR[bs[0]:bs[1]])

    dFsel = array(dFsel)
    dRselF = array(dRsel)
    
    if plot:
        figure(pname+"F")

    # Calculating the S score
    Fscore = burstscore(burstIdF,dF,1.0,idx= idxF)
    
    # We stroe the score per time  
    stats["burstIdF score"] = Fscore

    if verbose:
        print("Particle "+pname)
        stats["Selected dF/F"] = (mean(dFsel),std(dFsel)/sqrt(len(dFsel))) 
        print("Selected dF/F = %.2e +/- %.2e" % stats["Selected dF/F"])
        stats["Unselected dF/F"] = (mean(dF),std(dF)/sqrt(len(dF)))
        print("Total dF/F = %.2e +/- %.2e" % stats["Unselected dF/F"])
        stats["SelecteddF dR"] = (mean(dRsel),std(dRsel)/sqrt(len(dRsel)))
        print("Selected (from dF) dR = %.2e +/- %.2e" % stats["SelecteddF dR"])
        stats["Unselected dR"] = (mean(dR),std(dR)/sqrt(len(dR)))
        print("Unselected dR = %.2e +/- %.2e" % stats["Unselected dR"])

    # Burst search for R
    bselm, statsR_Brstsrch = burstsearch(dR,thf,ncycles,returnstat=True)
    burstIdR = array(bselm)
    stats["burstIdR"] = burstIdR
    stats["burstIdR stats"] = statsR_Brstsrch
    
    dRsel = []
    for bs in bselm:
        dRsel.extend(dR[bs[0]:bs[1]])
        
    if plot:
        figure(pname+"R")
    
    # ~ print(burstIdR,dR.shape)

    #Score for the ratio R
    Rscore =   burstscore(burstIdR,dR,1,idx= idxF)

    stats["burstIdR score"] = Rscore

    # Last statistics
    if verbose:
        stats["Selected dR"] = (mean(dRsel),std(dRsel)/sqrt(len(dRsel)))
        print("Selected (from dR) dR = %.2e +/- %.2e" % stats["Selected dR"])

        b = burstIdF*1.0
        db = b[:,1]-b[:,0]
        Sscore = 1.0*Fscore
        print("\n Burst statistics dF/F")
        stats["F number"] = b.shape[0]
        print("Number of intervals = ",b.shape[0])
        stats["F Av length"] = (mean(db),std(db)/sqrt(db.shape[0]))
        print("Av. Length  (cycles) = %.2f +- %.2f" % (mean(db),std(db)/sqrt(db.shape[0])))
        stats["F Total length"] = db.sum()
        print("Total length  (cycles) = ",db.sum())
        stats["F Burst score"] = (mean(Sscore[:,0]),std(Sscore[:,0])/sqrt(Sscore.shape[0]))
        print("Burst score =  %.4e +- %.4e" % (mean(Sscore[:,0]),std(Sscore[:,0])/sqrt(Sscore.shape[0])))
        st2 = sqrt(Sscore[:,1]/Sscore[:,2]-(Sscore[:,0]/Sscore[:,2])**2)/sqrt(Sscore.shape[0])
        mt2 = sum((Sscore[:,0]/Sscore[:,2])/st2**2)/(sum(1./st2**2))
        et2 = sqrt(1.0/sum(1./st2**2))
        stats["F Burst score per cycle"] = (mt2,et2)
        print("Burst score per cycle=  %.4e +- %.4e" % (mt2,et2))    
        st2 = sqrt(Sscore[:,1]/Sscore[:,3]-(Sscore[:,0]/Sscore[:,3])**2)/sqrt(Sscore.shape[0])
        mt2 = sum((Sscore[:,0]/Sscore[:,3])/st2**2)/(sum(1./st2**2))
        et2 = sqrt(1.0/sum(1./st2**2))
        stats["F Burst score per frame"] = (mt2,et2)
        print("Burst score per frame= %.2e +- %.2e" % (mt2,et2))    
        stats["F skipped time"] = mean((Sscore[:,3]-2*Sscore[:,2]) /Sscore[:,3])
        print("fraction skipped time = %.4f" %  mean((Sscore[:,3]-2*Sscore[:,2]) /Sscore[:,3]))

        print("\n Burst statistics dR")
        b = burstIdR*1.0
        db = b[:,1]-b[:,0]
        Sscore = 1.0*Rscore
        stats["R number"] = b.shape[0]
        print("Number of intervals = ",b.shape[0])
        stats["R Av length "] = (mean(db),std(db)/sqrt(db.shape[0]))
        print("Av. Length  (cycles)=  %.4f +- %.4f" % (mean(db),std(db)/sqrt(db.shape[0])))
        stats["R Total length"] = db.sum()
        print("Total length (cycles) = ",db.sum())
        stats["R Burst score"] = (mean(Sscore[:,0]),std(Sscore[:,0])/sqrt(Sscore.shape[0]))
        print("Burst score =  %.4e +- %.4e" % (mean(Sscore[:,0]),std(Sscore[:,0])/sqrt(Sscore.shape[0])))
        st2 = sqrt(Sscore[:,1]/Sscore[:,2]-(Sscore[:,0]/Sscore[:,2])**2)/sqrt(Sscore.shape[0])
        mt2 = sum((Sscore[:,0]/Sscore[:,2])/st2**2)/(sum(1./st2**2))
        et2 = sqrt(1.0/sum(1./st2**2))
        stats["R Burst score per cycle"] = (mt2,et2)
        print("Burst score per cycle=  %.4e +- %.4e" % (mt2,et2))    
        st2 = sqrt(Sscore[:,1]/Sscore[:,3]-(Sscore[:,0]/Sscore[:,3])**2)/sqrt(Sscore.shape[0])
        mt2 = sum((Sscore[:,0]/Sscore[:,3])/st2**2)/(sum(1./st2**2))
        et2 = sqrt(1.0/sum(1./st2**2))
        stats["R Burst score per frame"] = (mt2,et2)
        print("Burst score per frame= %.2e +- %.2e" % (mt2,et2) )   
        stats["R skipped time"] = mean((Sscore[:,3]-2*Sscore[:,2]) /Sscore[:,3])
        print("fraction skipped time = %.4f" % mean( (Sscore[:,3]-2*Sscore[:,2]) /Sscore[:,3] ))
        
        stats["overlap"] = overlapintervals(burstIdF,burstIdR)
        print("\n The overlap between the 2 burst searchs is %.2f %%" % round(100*stats["overlap"][0],2))
        print("\n")
        print(12*"#"+"\n\n")
    
    return(stats)


def shadeBurstIntervals(burst,ts):
    '''burst = stats["burstIdR"]
    '''
    #print("Burst R:",burst)
    for b in burst:
        tsb = ts[b[0]:b[1],:]
        nsel = arange(tsb.shape[0]-1)[(tsb[1:]-tsb[:-1]-0.01>0.0005)]
        s0 = 0
        if nsel.shape[0]>0:
            s0 = 0
            for s1 in nsel:
                currentAxis.add_patch(Rectangle((tsb[s0], ymin), tsb[s1]-tsb[s0], ymax,alpha=0.7,color='deeppink'))
                s0 = s1+1
        else:
            currentAxis.add_patch(Rectangle((tsb[0], ymin), tsb[-1]-tsb[0], ymax,alpha=0.7,color='deeppink'))







def main(mainfile = "allROIsSignal_OnState_wave.dat",ctrlfile = None,outputfile = "fullstats.dat"):
    #ctrlfile = "allROIsSignal_OnState_contrl.dat"
    wave = loadtxt(mainfile)
    print("Wave shape",wave.shape)

    if ctrlfile is not None:
        ctrl = loadtxt(ctrlfile)
        print("Ctrl shape",ctrl.shape)

    
    nrois2, nframes = wave.shape
    nrois = nrois2//2

    thf = 0.6
    
    fullstats = {}

    for i in range(nrois):
        try:
            FB = 1.0*wave[i,:] # Blue channel
            FA = 1.0*wave[i+nrois,:] # red channel
            if isfinite(FA).sum()>0:
                stats = burstanalysis(FB,FA, verbose = True, plot = False,pname=str(i).zfill(3),thf = 0.6, ncycles = 12)
            else:
                stats = nan
            if ctrlfile is not None:
                FB = 1.0*ctrl[i,:] # Blue channel
                FA = 1.0*ctrl[i+nrois,:] # red channel
                if isfinite(FA).sum()>0:
                    statsC = burstanalysis(FB,FA, verbose = True, plot = False,pname=str(i).zfill(3)+" ctrl",thf = 0.6, ncycles = 12)
                else:
                    statsC = nan
                
                fullstats[i] = [stats,statsC]
            else:
                fullstats[i] = [stats]
        except  Exception as e:
            fullstats[i] = [nan,nan]
            print("For spine i=",i," ",e)
            
        if i>0 and i%100==0:      
            with open("fullstatstemp.dat", 'wb') as filehandler: 
                pickle.dump(fullstats, filehandler)

    with open(outputfile, 'wb') as filehandler: 
        pickle.dump(fullstats, filehandler)

def plottingF(mainfile = "allROIsSignal_OnState_wave.dat",ctrlfile = None,outputfolder = "imgs/",framespercycle = 4):
    wave = loadtxt(mainfile)
    print("Wave shape",wave.shape)
    nrois2, nframes = wave.shape
    nrois = nrois2//2
    for i in range(nrois):
        #print(i)
        FB = 1.0*wave[i,:] # Blue channel
        FA = 1.0*wave[i+nrois,:] # red channel
        if isfinite(FA).sum()>0:
            stats = burstanalysis(FB,FA, verbose = True, plot = False,pname=str(i).zfill(3),thf = 0.6, ncycles = 12,framespercycle = framespercycle)
        else:
            stats = nan
        
        framespercycle2 = framespercycle//2
        t = arange(0,FA.shape[0])*10e-3
        Sast = FA.reshape(FA.shape[0]//framespercycle,framespercycle)
        sel = isfinite(Sast.sum(axis=1))
        #print(sel.sum())
        Sp = FA+FB
        Stt = Sp.reshape(Sp.shape[0]//framespercycle,framespercycle)
        tst = t.reshape(t.shape[0]//framespercycle,framespercycle)
        ts = tst[sel,:]
        tse = column_stack((tst[:,:framespercycle2].mean(axis=1),tst[:,framespercycle2:].mean(axis=1)))
        tse = tse[sel,].flatten()
        Sps = column_stack((Stt[:,:framespercycle2].mean(axis=1),Stt[:,framespercycle2:].mean(axis=1)))
        Sps = Sps[sel,:].flatten()
        selnan = isnan(Sast.sum(axis=1))
        Stt[selnan,:] = array([nan]*framespercycle)
        tst[selnan,:] = array([nan]*framespercycle)
        tst = column_stack((tst[:,:framespercycle2].mean(axis=1),tst[:,framespercycle2:].mean(axis=1)))
        Stt = column_stack((Stt[:,:framespercycle2].mean(axis=1),Stt[:,framespercycle2:].mean(axis=1)))

        def thousands(x, pos):
            'The two args are the value and tick position'
            return '%1.0fk' % (x*1e-3)
        
        
        if Sps.shape[0]>0:
            formatter = FuncFormatter(thousands)


            fig = figure(figsize=(14,10))
            grid = plt.GridSpec(3, 3, wspace=0.4, hspace=0.2)
            plt.subplot(grid[:2, 0:])
            currentAxis = gca()
            plot(tst.flatten(),Stt.flatten(),'C0-',alpha=0.7)
            plot(tse,Sps,'C0.',alpha=0.7)
            #print(burst)
            ymin = min(Sps)*.98
            ymax = max(Sps)*1.02
            currentAxis.yaxis.set_major_formatter(formatter)


            #try:

            burst = stats["burstIdF"]
            for b in burst:
                tsb = ts[b[0]:b[1],:].flatten()
                nsel = arange(tsb.shape[0]-1)[(tsb[1:]-tsb[:-1]-0.01>0.0005)]
                if nsel.shape[0]>0:
                    #print(nsel)
                    s0 = 0
                    for s1 in nsel:
                        currentAxis.add_patch(Rectangle((tsb[s0], ymin), tsb[s1]-tsb[s0], ymax,alpha=0.7,color='gray'))
                        s0 = s1+1
                else:
                    currentAxis.add_patch(Rectangle((tsb[0], ymin), tsb[-1]-tsb[0], ymax,alpha=0.7,color='gray'))
            
            burst = stats["burstIdR"]
            #print("Burst R:",burst)
            for b in burst:
                tsb = ts[b[0]:b[1],:].flatten()
                nsel = arange(tsb.shape[0]-1)[(tsb[1:]-tsb[:-1]-0.01>0.0005)]
                s0 = 0
                if nsel.shape[0]>0:
                    s0 = 0
                    for s1 in nsel:
                        currentAxis.add_patch(Rectangle((tsb[s0], ymin), tsb[s1]-tsb[s0], ymax,alpha=0.7,color='deeppink'))
                        s0 = s1+1
                else:
                    currentAxis.add_patch(Rectangle((tsb[0], ymin), tsb[-1]-tsb[0], ymax,alpha=0.7,color='deeppink'))

                    
            xlabel("t (s)",fontsize = 16)
            ylabel("F (cps)",fontsize = 16)
            
            burst = stats["burstIdF"]
            bscore = stats["burstIdF score"]
            ab = bscore.argsort()
            if len(ab)>=3:
                for ib in range(3):
                    b = burst[ab][ib]
                    #print(db)
                    plt.subplot(grid[2, ib])
                    currentAxis = gca()
                    currentAxis.yaxis.set_major_formatter(formatter)
                    plot(tst.flatten(),Stt.flatten(),'C0-',alpha=0.7)
                    plot(tse,Sps,'C0.',alpha=0.7)
                    sqwup = (arange(Sps.shape[0])%2)==0
                    plot(tse[sqwup],Sps[sqwup],'C3.',alpha=0.7)
                    #print(Sps.shape)
                    tsb = ts[b[0]:b[1],:].flatten()
                    xlim(tsb[0],tsb[-1])
                    bottom, top = ylim()
                    text(tsb[0]+(tsb[-1]-tsb[0])*.02,top-(top-bottom)*.02,'Burst score: %.2e' % bscore[ab][ib], ha='left', va='top')
                    xlabel('t (s)')
                    tss = linspace(tsb[0],tsb[-1],200)
                    sawchain = 20*(floor(6*(tss/120.0e-3))%6)-20
                    plot(tss,bottom+(top-bottom)*0.005+(1+sign(sin(2*pi*tss*25.0)))*(top-bottom)*0.07/2.0,'k-')

                    
            elif len(ab)==2:
                for ib in range(2):
                    b = burst[ab][ib]
                    #print(db)
                    if ib==0:
                        plt.subplot(grid[2, :2])
                    else:
                        plt.subplot(grid[2,2])
                    currentAxis = gca()
                    currentAxis.yaxis.set_major_formatter(formatter)
                    plot(tst.flatten(),Stt.flatten(),'C0-',alpha=0.7)
                    plot(tse,Sps,'C0.',alpha=0.7)
                    sqwup = (arange(Sps.shape[0])%2)==0
                    plot(tse[sqwup],Sps[sqwup],'C3.',alpha=0.7)
                    #print(Sps.shape)
                    tsb = ts[b[0]:b[1],:].flatten()
                    xlim(tsb[0],tsb[-1])
                    xlabel('t (s)')
                    bottom, top = ylim()
                    text(tsb[0]+(tsb[-1]-tsb[0])*.02,top-(top-bottom)*.02,'Burst score: %.2e' % bscore[ab][ib], ha='left', va='top')
                    xlabel('t (s)')
                    tss = linspace(tsb[0],tsb[-1],200)
                    plot(tss,bottom+(top-bottom)*0.005+(1+sign(sin(2*pi*tss*25.0)))*(top-bottom)*0.07/2.0,'k-')
            elif len(ab) == 1:
                for ib in range(1):
                    b = burst[ab][ib]
                    #print(db)
                    plt.subplot(grid[2,:])
                    currentAxis = gca()
                    currentAxis.yaxis.set_major_formatter(formatter)
                    plot(tst.flatten(),Stt.flatten(),'C0-',alpha=0.7)
                    plot(tse,Sps,'C0.',alpha=0.7)
                    sqwup = (arange(Sps.shape[0])%2)==0
                    plot(tse[sqwup],Sps[sqwup],'C3.',alpha=0.7)
                    #print(Sps.shape)
                    tsb = ts[b[0]:b[1],:].flatten()
                    xlim(tsb[0],tsb[-1])
                    xlabel('t (s)')
                    bottom, top = ylim()
                    text(tsb[0]+(tsb[-1]-tsb[0])*.02,top-(top-bottom)*.02,'Burst score: %.2e' % bscore[ab][ib], ha='left', va='top')
                    xlabel('t (s)')
                    tss = linspace(tsb[0],tsb[-1],200)
                    plot(tss,bottom+(top-bottom)*0.005+(1+sign(sin(2*pi*tss*25.0)))*(top-bottom)*0.07/2.0,'k-')

            #except:
            #    pass
            savefig(outputfolder+"BurstF"+str(i).zfill(4)+".png")
            close()


def plottingR(mainfile = "allROIsSignal_OnState_wave.dat",ctrlfile = None,outputfolder = "imgs/", framespercycle=4):
    wave = loadtxt(mainfile)
    print("Wave shape",wave.shape)
    nrois2, nframes = wave.shape
    nrois = nrois2//2
    for i in range(nrois):
        #print(i)
        FB = 1.0*wave[i,:] # Blue channel
        FA = 1.0*wave[i+nrois,:] # red channel
        if isfinite(FA).sum()>0:
            stats = burstanalysis(FB,FA, verbose = True, plot = False,pname=str(i).zfill(3),thf = 0.6, ncycles = 12, framespercycle= framespercycle)
        else:
            stats = nan
        
        if mean(FA[isfinite(FA)])>=mean(FB[isfinite(FB)]):
            SXt = FA.reshape(FA.shape[0]//framespercycle,framespercycle)
        else:
            SXt = FB.reshape(FB.shape[0]//framespercycle,framespercycle)    
        
        framespercycle2 = framespercycle//2
        t = arange(0,FA.shape[0])*10e-3
        Sast = FA.reshape(FA.shape[0]//framespercycle,framespercycle)
        sel = isfinite(Sast.sum(axis=1))
        #print(sel.sum())
        Sp = FA+FB
        Stt = Sp.reshape(Sp.shape[0]//framespercycle,framespercycle)
        SXt = SXt/Stt
        
        tst = t.reshape(t.shape[0]//framespercycle,framespercycle)
        ts = tst[sel,:]
        tse = column_stack((tst[:,:framespercycle2].mean(axis=1),tst[:,framespercycle2:].mean(axis=1)))
        tse = tse[sel,].flatten()
        
        Sps = column_stack((SXt[:,:framespercycle2].mean(axis=1),SXt[:,framespercycle2:].mean(axis=1)))
        Sps = Sps[sel,:].flatten()

        selnan = isnan(SXt.sum(axis=1))
        SXt[selnan,:] = array([nan]*framespercycle)
        tst[selnan,:] = array([nan]*framespercycle)
        tst = column_stack((tst[:,:framespercycle2].mean(axis=1),tst[:,framespercycle2:].mean(axis=1)))
        Stt = column_stack((SXt[:,:framespercycle2].mean(axis=1),SXt[:,framespercycle2:].mean(axis=1)))
       
        
        if Sps.shape[0]>0:
            fig = figure(figsize=(14,10))
            grid = plt.GridSpec(3, 3, wspace=0.2, hspace=0.3)
            plt.subplot(grid[:2, 0:])
            currentAxis = gca()
            plot(tst.flatten(),Stt.flatten(),'C0-',alpha=0.7)
            plot(tse,Sps,'C0.',alpha=0.7)
            #print(burst)
            miSps = Sps.min()
            maSps = Sps.max()
            ymin = miSps-(maSps-miSps)*0.02
            ymax = maSps+(maSps-miSps)*0.02
            #currentAxis.yaxis.set_major_formatter(formatter)

            burst = stats["burstIdF"]
            for b in burst:
                tsb = ts[b[0]:b[1],:].flatten()
                nsel = arange(tsb.shape[0]-1)[(tsb[1:]-tsb[:-1]-0.01>0.0005)]
                if nsel.shape[0]>0:
                    #print(nsel)
                    s0 = 0
                    for s1 in nsel:
                        currentAxis.add_patch(Rectangle((tsb[s0], ymin), tsb[s1]-tsb[s0], ymax,alpha=0.7,color='gray'))
                        s0 = s1+1
                else:
                    currentAxis.add_patch(Rectangle((tsb[0], ymin), tsb[-1]-tsb[0], ymax,alpha=0.7,color='gray'))
            
            burst = stats["burstIdR"]
            #print("Burst R:",burst)
            for b in burst:
                tsb = ts[b[0]:b[1],:].flatten()
                nsel = arange(tsb.shape[0]-1)[(tsb[1:]-tsb[:-1]-0.01>0.0005)]
                s0 = 0
                if nsel.shape[0]>0:
                    s0 = 0
                    for s1 in nsel:
                        currentAxis.add_patch(Rectangle((tsb[s0], ymin), tsb[s1]-tsb[s0], ymax,alpha=0.7,color='deeppink'))
                        s0 = s1+1
                else:
                    currentAxis.add_patch(Rectangle((tsb[0], ymin), tsb[-1]-tsb[0], ymax,alpha=0.7,color='deeppink'))
                    
            xlabel("t (s)",fontsize = 16)
            ylabel("R",fontsize = 16)
            
            burst = stats["burstIdR"]
            bscore = stats["burstIdR score"]
            ab = bscore.argsort()
            if len(ab)>=3:
                for ib in range(3):
                    b = burst[ab][ib]
                    #print(db)
                    plt.subplot(grid[2, ib])
                    currentAxis = gca()
                    #currentAxis.yaxis.set_major_formatter(formatter)
                    plot(tst.flatten(),Stt.flatten(),'C0-',alpha=0.7)
                    plot(tse,Sps,'C0.',alpha=0.7)
                    sqwup = (arange(Sps.shape[0])%2)==0
                    plot(tse[sqwup],Sps[sqwup],'C3.',alpha=0.7)
                    #print(Sps.shape)
                    tsb = ts[b[0]:b[1],:].flatten()
                    xlim(tsb[0],tsb[-1])
                    bottom, top = ylim()
                    text(tsb[0]+(tsb[-1]-tsb[0])*.02,top-(top-bottom)*.02,'Burst score: %.2e' % bscore[ab][ib], ha='left', va='top')
                    xlabel('t (s)')
                    tss = linspace(tsb[0],tsb[-1],2000)
                    plot(tss,bottom+(top-bottom)*0.005+(1+sign(sin(2*pi*tss*25.0)))*(top-bottom)*0.07/2.0,'k-')
            elif len(ab)==2:
                for ib in range(2):
                    b = burst[ab][ib]
                    #print(db)
                    if ib==0:
                        plt.subplot(grid[2, :2])
                    else:
                        plt.subplot(grid[2,2])
                    currentAxis = gca()
                    #currentAxis.yaxis.set_major_formatter(formatter)
                    plot(tst.flatten(),Stt.flatten(),'C0-',alpha=0.7)
                    plot(tse,Sps,'C0.',alpha=0.7)
                    sqwup = (arange(Sps.shape[0])%2)==0
                    plot(tse[sqwup],Sps[sqwup],'C3.',alpha=0.7)
                    #print(Sps.shape)
                    tsb = ts[b[0]:b[1],:].flatten()
                    xlim(tsb[0],tsb[-1])
                    xlabel('t (s)')
                    bottom, top = ylim()
                    text(tsb[0]+(tsb[-1]-tsb[0])*.02,top-(top-bottom)*.02,'Burst score: %.2e' % bscore[ab][ib], ha='left', va='top')
                    xlabel('t (s)')
                    tss = linspace(tsb[0],tsb[-1],2000)
                    plot(tss,bottom+(top-bottom)*0.005+(1+sign(sin(2*pi*tss*25.0)))*(top-bottom)*0.07/2.0,'k-')
            elif len(ab) == 1:
                for ib in range(1):
                    b = burst[ab][ib]
                    #print(db)
                    plt.subplot(grid[2,:])
                    currentAxis = gca()
                    #currentAxis.yaxis.set_major_formatter(formatter)
                    plot(tst.flatten(),Stt.flatten(),'C0-',alpha=0.7)
                    plot(tse,Sps,'C0.',alpha=0.7)
                    sqwup = (arange(Sps.shape[0])%2)==0
                    plot(tse[sqwup],Sps[sqwup],'C3.',alpha=0.7)
                    #print(Sps.shape)
                    tsb = ts[b[0]:b[1],:].flatten()
                    xlim(tsb[0],tsb[-1])
                    xlabel('t (s)')
                    bottom, top = ylim()
                    text(tsb[0]+(tsb[-1]-tsb[0])*.02,top-(top-bottom)*.02,'Burst score: %.2e' % bscore[ab][ib], ha='left', va='top')
                    xlabel('t (s)')
                    tss = linspace(tsb[0],tsb[-1],2000)
                    plot(tss,bottom+(top-bottom)*0.005+(1+sign(sin(2*pi*tss*25.0)))*(top-bottom)*0.07/2.0,'k-')

            savefig(outputfolder+"BurstR"+str(i).zfill(4)+".png")
            close()

if __name__ == '__main__':
    main(mainfile = "allROIsSignal_OnState_wave.dat",ctrlfile = "allROIsSignal_OnState_contrl.dat")
        
