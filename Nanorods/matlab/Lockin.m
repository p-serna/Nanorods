# Sine-wave with frequency freq and dephase angle at time t.
function f = fsin(t, freq, dephase)
    f = sin(2.0*pi*freq*t+phi); 
    end

# St: trace
# t: time array
# freq: target freq
# tint: integration time
# phi: dephasing angle... leaving at zero unless you want to have it locked to a specific angle

function [tm,r,phi] =  product(St,t,freq,tint,phi)
    sxsine1= St*fsin(t,freq,phi)
    sxsine2= St*fsin(t,freq,phi+pi/2.0)
    sel = t<max(t)-tint
    tm = t[sel]
    x = tm*0
    y = tm*0
    for i,tu in enumerate(tm):
        sel = (t>tu)*(t<tu+T)
        st = sxsine1[sel]
        tt = t[sel]
        x[i] = integrate.simps(st,tt)
        st = sxsine2[sel]
        y[i] = integrate.simps(st,tt)
    if polar:
        return (tm,sqrt(x**2+y**2),arctan2(x,y))
    else:
        return (tm,x,y)
    
    
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
        bsel= [bsel;[ssel(i,1),ssel(i,2)+ncycles]];
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


# How to use them: Once you have dF/F trace
# burst = burstsearch(dF,0.6,12,1) will provide the intervals
# with 0.6 a threshold parameter, and 12 number of cycles it averages.
