function xavg = running_mean(x, N)
    N2 = floor(N/2);
    cumsumt = cumsum([zeros(1,N2) x zeros(1,N2)]); 
    xavg = (cumsumt(N+1:end) - cumsumt(1:end-N)) / N;
    xavg(1:N2) = xavg(1:N2)*N./(N2:N-1);
    xavg(end-N2+1:end) = xavg(end-N2+1:end)*N./(N-1:-1:N2);
    end

ncycles = 12;
function burst =  burstsearch(dS,thf,ncycles,merged)
    # First average of dS over ncycles
    Ssc = (running_mean(dS,ncycles)).^2;

    # We set a threshold to keep only the tail of the distribution
    threshold = (max(Ssc)-min(Ssc))*thf+min(Ssc);
    sel = (1:size(Ssc)(2))(Ssc>threshold);
    
    # This is the way to select separate intervals:
    # If diference in numbering is more than 1 then the intervals are disconnected
    csel = sel(2:end)-sel(1:end-1);
    s0 = 1;
    ssel = [];
    for i=1:size(csel)(2);
        si =csel(i);
        if si>1
            ssel= [ssel;[sel(s0),sel(i)]];
            s0 = i+1;
        endif
    endfor
    ssel= [ssel;[sel(s0),sel(end)]];

    # Instead of the numbering on dS, we want to have the numbering in the cycles,
    # which is the index of the point selected + 12
    bsel = [];
    for i = 1:size(ssel)(1)
        bsel= [bsel;[ssel(i,1),ssel(i,2)+12]];
    endfor
    
    # This is to merge intervals that overlap, if merged option is 1.
    bselm = [];
    if merged == 1
        bs0a = bsel(1,1);
        bs1a = bsel(1,2);
        bsize = size(bsel)(1);
        for ibs = 2:bsize
            bs0 = bsel(ibs,1);
            bs1 = bsel(ibs,2);
            if bs1a< bs0
                bselm =[ bselm; [bs0a,bs1a]];
                bs0a = bs0;
            endif
            bs1a = bs1;
        endfor
        bselm= [bselm; [bs0a,bs1a]];
    else
        bselm = bsel;
    endif
    
    burst = bselm;
end
    
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


def burstanalysis(FB,FA, verbose = True, plot = False,pname = "",thf = 0.6, ncycles = 12):    
    # Printing statistics
    stats = {}

    # Total trace, F
    St = FA+FB
    Fav = St[isfinite(St)].mean()
    # ~ SA = FA[isfinite(FA)]
    # ~ SB = FB[isfinite(FB)]

    # Reshape of total trace in 4 columns, its shape is now  (ncycles,4) 
    Stt = St.reshape(St.shape[0]//4,4)
    # ~ # Stt2: we take averages for every 2 frames
    # ~ Stt2 = column_stack((Stt[:,:2].mean(axis=1),Stt[:,2:].mean(axis=1)))
    
    # ratio R of one channel over F. We select the channel with larger average
    if mean(FA[isfinite(FA)])>=mean(FB[isfinite(FB)]):
        SXt = FA.reshape(FA.shape[0]//4,4)
    else:
        SXt = FB.reshape(FB.shape[0]//4,4)    
    
    SXt = SXt/Stt
    # ~ # Now we take averages 0,1 and 2,3 for each interval
    # ~ SXt = column_stack((SXt[:,:2].mean(axis=1),SXt[:,2:].mean(axis=1)))
    
    # Definition of dR and dF, for each interval: R_0+R_1-R_2-R_3   
    dR = (SXt[:,:2].mean(axis=1)-SXt[:,2:].mean(axis=1))
    dF = (Stt[:,:2].mean(axis=1)-Stt[:,2:].mean(axis=1))

    # We remove intervals where there was at least a frame classified as blinking (dark state)    
    dR = dR[isfinite(dR)]
    SXt = SXt[isfinite(SXt.sum(axis=1)),:]
    dF = dF[isfinite(dF)]
    Stt = Stt[isfinite(Stt.sum(axis=1)),:]

    # dF will be defined as dF/<F> with <F> average of the whole trace.
    dF = dF[isfinite(dF)]/Fav

    # Burst search : it provides merged intervals [ [idx00,idxf0], [idx01,idxf1], ...]
    bselm = burstsearch(dF,thf,ncycles)
    burstIdF = array(bselm) 
    stats["burstIdF"] = burstIdF
    # Now we take the parts of the trace selected by the burst search
    dFsel = []
    Sttsel = []
    dRsel = []
    for bs in bselm:
        dFsel.extend(dF[bs[0]:bs[1]])
        Sttsel.extend(Stt[bs[0]:bs[1],:])
        dRsel.extend(dR[bs[0]:bs[1]])

    dFsel = array(dFsel)
    Sttsel = array(Sttsel)
    dRselF = array(dRsel)
    
    if plot:
        figure(pname+"F")

    # Calculating the S score
    Sscore = []
    S2score = []
    Snscore  = []
    ts0 = 0
    for bs in bselm:
        # time frame of the interval from 0 to length*4
        ts = ts0+arange((bs[1]-bs[0])*4)
        # We select the cycles from bs[0] to bs[1]
        Stts = Stt[bs[0]:bs[1],:]
        Stts = column_stack((Stts[:,:2].mean(axis=1),Stts[:,2:].mean(axis=1)))
        
        # Burst score, S= sum dF/F, and burst score per time 
        # (S2 and Sn are there to calculate error bars later for S/time) 
        Sscore.append((Stts[:,0]-Stts[:,1]).sum()/Fav)
        S2score.append(((Stts[:,0]-Stts[:,1])**2).sum()/Fav**2)
        Snscore.append(len(Stts[:,0]))
        
        # Plotting the trace
        if plot:
            Stts = Stts.flatten()
            idx = (ts-ts0)%4
            sel2 = (idx==0)+(idx==2)
            plot(ts[sel2],Stts,'C0.-',alpha=0.5)
            sel = idx[sel2]//2==1
            plot(ts[sel2][sel],Stts[sel],'C3.',alpha=0.5)
            vlines(ts0,min(Stts),max(Stts),linestyles='--',alpha=0.5)
        # To keep track of the time
        ts0 = ts[-1]+1
    
    if plot:
        xlabel("t(frames)")
        ylabel("F(cps)")        
    
    # We make them arrays    
    Sscore = array(Sscore)
    S2score = array(S2score)
    Snscore = array(Snscore)
    stats["burstIdF score"] = Sscore/Snscore


    # ~ tempup = Sttsel[:,:2].flatten()/Sttsel.mean()
    # ~ temd = Sttsel[:,2:].flatten()/Sttsel.mean()
    
    # ~ figure(red[:-8]+"a")
    # ~ title("dF/F")
    # ~ dSsel = array(dSsel)
    # ~ ts = arange(len(dSsel))
    # ~ selu = dSsel>0
    # ~ plot(ts,ts*0+0.0,'k--')
    # ~ plot(ts,dSsel,'.',alpha=0.5)   
    # ~ plot(ts[selu],dSsel[selu],'C1.',alpha=0.5)
    
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
    bselm = burstsearch(dR,thf,ncycles)
    burstIdR = array(bselm)
    stats["burstIdR"] = burstIdR

    dRsel = []
    Rtsel = []
    for bs in bselm:
        dRsel.extend(dR[bs[0]:bs[1]])
        Rtsel.extend(SXt[bs[0]:bs[1],:])
    Rtt = array(Rtsel)
        
    if plot:
        figure(pname+"R")
    
    #Score for the ratio R
    Rscore = []
    R2score = []
    Rnscore = []
    ts0 = 0
    for bs in bselm:
        ts = ts0+arange((bs[1]-bs[0])*4)
        
        Stts = SXt[bs[0]:bs[1],:]
        Stts = column_stack((Stts[:,:2].mean(axis=1),Stts[:,2:].mean(axis=1)))
        
        # Burst score, S= sum dR (in this case without norm factor, it could explode), 
        #               and also burst score per time 
        # (S2 and Sn are there to calculate error bars later for S/time) 
        Rscore.append((Stts[:,0]-Stts[:,1]).sum())
        R2score.append(((Stts[:,0]-Stts[:,1])**2).sum())
        Rnscore.append(len(Stts[:,0]))
        if plot:
            Stts = Stts.flatten()
            idx = (ts-ts0)%4
            sel2 = (idx==0)+(idx==2)
            plot(ts[sel2],Stts,'C0.-',alpha=0.5)
            sel = idx[sel2]//2==1             
            plot(ts[sel2][sel],Stts[sel],'C3.',alpha=0.5)
            vlines(ts0,min(Stts),max(Stts),linestyles='--',alpha=0.5)
        ts0 = ts[-1]+1
    if plot:
        xlabel("t(frames)")
        ylabel("X/(R+B) (X=R if <R> > <B> or X=B otherwise)")        

    Rscore = array(Rscore)
    R2score = array(R2score)
    Rnscore = array(Rnscore)
    stats["burstIdR score"] = Rscore/Rnscore
    
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

    # Last statistics
    if verbose:
        stats["Selected dR"] = (mean(dRsel),std(dRsel)/sqrt(len(dRsel)))
        print("Selected (from dR) dR = %.2e +/- %.2e" % stats["Selected dR"])

        b = burstIdF*1.0
        db = b[:,1]-b[:,0]
        print("\n Burst statistics dF/F")
        stats["F number"] = b.shape[0]
        print("Number of intervals = ",b.shape[0])
        stats["F Av length"] = (mean(db),std(db)/sqrt(db.shape[0]))
        print("Av. Length = ",mean(db),"+/-",std(db)/sqrt(db.shape[0]))
        stats["F Total length"] = db.sum()
        print("Total length = ",db.sum())
        stats["F Burst score"] = (mean(Sscore),std(Sscore)/sqrt(len(Sscore)))
        print("Burst score = ",mean(Sscore),"+-",std(Sscore)/sqrt(len(Sscore)))
        st2 = sqrt(S2score/Snscore-(Sscore/Snscore)**2)/sqrt(len(Snscore))
        mt2 = sum((Sscore/Snscore)/st2**2)/(sum(1./st2**2))
        et2 = sqrt(1.0/sum(1./st2**2))
        stats["F Burst score per time"] = (mt2,et2)
        print("Burst score per time (units of frame)= ",mt2,"+- ",et2)    

        b = burstIdR*1.0
        db = b[:,1]-b[:,0]
        stats["R number"] = b.shape[0]
        print("Number of intervals = ",b.shape[0])
        stats["R Av length"] = (mean(db),std(db)/sqrt(db.shape[0]))
        print("Av. Length = ",mean(db),"+/-",std(db)/sqrt(db.shape[0]))
        stats["R Total length"] = db.sum()
        print("Total length = ",db.sum())
        stats["R Burst score"] = (mean(Sscore),std(Sscore)/sqrt(len(Sscore)))
        print("Burst score = ",mean(Sscore),"+-",std(Sscore)/sqrt(len(Sscore)))
        st2 = sqrt(S2score/Snscore-(Sscore/Snscore)**2)/sqrt(len(Snscore))
        mt2 = sum((Sscore/Snscore)/st2**2)/(sum(1./st2**2))
        et2 = sqrt(1.0/sum(1./st2**2))
        stats["R Burst score per time"] = (mt2,et2)
        print("Burst score per time (units of frame)= ",mt2,"+- ",et2) 
        
        stats["overlap"] = overlapintervals(burstIdF,burstIdR)
        print("\n The overlap between the 2 burst searchs is %.2f %%" % round(100*stats["overlap"][0],2))
        print("\n")
        print(12*"#"+"\n\n")
    
    return(stats)

def main(mainfile = "allROIsSignal_OnState_wave.dat",ctrlfile = None,outputfile = "fullstats.dat"):
    #ctrlfile = "allROIsSignal_OnState_contrl.dat"
    wave = loadtxt(mainfile)
    print("Wave shape",wave.shape)

    if ctrlfile is not None:
        ctrl = loadtxt(ctrfile)
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
        except:
            fullstats[i] = [nan,nan]
            
        if i>0 and i%100==0:      
            with open("fullstatstemp.dat", 'wb') as filehandler: 
                pickle.dump(fullstats, filehandler)

    with open(outputfile, 'wb') as filehandler: 
        pickle.dump(fullstats, filehandler)

if __name__ == '__main__':
    main()
        
