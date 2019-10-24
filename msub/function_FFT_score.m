function  score_FFT = FFT_score(signal, length_trh, int, Tacq)
NROI=size(signal,2);

for roi=1:NROI
    if strcmp(class(signal),'cell')
        sig=signal{roi};
    elseif NROI==1
        sig=signal(:,roi)';
    else
        sig=signal(:,roi);
    end
    
    Nframes=size(sig,2);
        
    if Nframes>length_trh
        Y = fft(sig);
        P2 = abs(Y/length(sig));
        P1 = P2(1:length(sig)/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = (1/Tacq)*(0:(length(sig)/2))/length(sig);
        [~, ind_start]=min(abs(f-0.99/(4*Tacq)));
        [~, ind_stop]=min(abs(f-1.01/(4*Tacq)));
        [~, inx_start]=min(abs(f-(1/(4*Tacq)-int/2)));
        [~, inx_stop]=min(abs(f-(1/(4*Tacq)+int/2)));
        m=1;
        for i=ind_start:ind_stop
            temp1(m)=(P1(i)-mean(P1(inx_start:inx_stop)))/std(P1(inx_start:inx_stop));
            m=m+1;
        end
        score_FFT(roi)=max(temp1);
    else
        score_FFT(roi)=NaN;
        
    end
    clearvars temp1
end
