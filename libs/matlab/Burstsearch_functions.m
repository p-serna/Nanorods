# running mean with padding = "valid".
function xavg = running_mean(x, N)
    cumsumt = cumsum([0 x]); 
    xavg = (cumsumt(N+1:end) - cumsumt(1:end-N)) / N;
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
