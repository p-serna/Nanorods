clearvars Signal_wave Signal_contrl NROI Tacq Nframes Mask_CMOS SumMovie_wave SumMovie_contrl  pts_blue pts_red ROI_contrl
firstrun=1;
Keep_Only_Double=1;
for exp=1:length(Exp_Data)
    if firstrun
        Signal_wave_blue=Exp_Data(exp).Signal_wave(:,1:Exp_Data(exp).NROI);
        Signal_wave_red=Exp_Data(exp).Signal_wave(:,Exp_Data(exp).NROI+1:end);   
        Signal_contrl_blue=Exp_Data(exp).Signal_contrl(:,1:Exp_Data(exp).NROI);
        Signal_contrl_red=Exp_Data(exp).Signal_contrl(:,Exp_Data(exp).NROI+1:end);
        NROI_all=Exp_Data(exp).NROI;
        Tacq_all=Exp_Data(exp).Tacq;
        Nframes_all=Exp_Data(exp).Nframes;
        Mask_CMOS=Exp_Data(exp).Mask_CMOS;
        SumMovie_wave=Exp_Data(exp).SumMovie_wave;
        SumMovie_contrl=Exp_Data(exp).SumMovie_contrl; 
        firstrun=0;
    else
        Signal_wave_blue=cat(2, Signal_wave_blue, Exp_Data(exp).Signal_wave(:,1:Exp_Data(exp).NROI));
        Signal_wave_red=cat(2, Signal_wave_red, Exp_Data(exp).Signal_wave(:,Exp_Data(exp).NROI+1:end));
        Signal_contrl_blue=cat(2, Signal_contrl_blue, Exp_Data(exp).Signal_contrl(:,1:Exp_Data(exp).NROI));
        Signal_contrl_red=cat(2, Signal_contrl_red, Exp_Data(exp).Signal_contrl(:,Exp_Data(exp).NROI+1:end));
        NROI_all=cat(1, NROI_all, Exp_Data(exp).NROI);
        Tacq_all=cat(1, Tacq_all, Exp_Data(exp).Tacq);
        Nframes_all=cat(1, Nframes_all, Exp_Data(exp).Nframes);
        Mask_CMOS=cat(3, Mask_CMOS, Exp_Data(exp).Mask_CMOS);
        SumMovie_wave=cat(3, SumMovie_wave, Exp_Data(exp).SumMovie_wave);
        SumMovie_contrl=cat(3, SumMovie_contrl, Exp_Data(exp).SumMovie_contrl);

    end
    pts_blue{exp}=Exp_Data(exp).pts_blue;
    pts_red{exp}=Exp_Data(exp).pts_red;
    ROI_contrl{exp}=Exp_Data(exp).ROI_contrl;
end

clearvars firstrun
Signal_wave=[Signal_wave_blue Signal_wave_red];
Signal_contrl=[Signal_contrl_blue Signal_contrl_red];
NROI=sum(NROI_all);
Nframes=size(Signal_wave,1);
Tacq=mean(Tacq_all);

Signal_2ch_wave=Signal_wave(:,1:NROI)+Signal_wave(:,NROI+1:end);
Signal_2ch_contrl=Signal_contrl(:,1:NROI)+Signal_contrl(:,NROI+1:end);

clearvars Signal_wave_blue Signal_wave_red Signal_contrl_blue Signal_contrl_red
%% calculation of correlation coefficient between blue and red traces
clearvars CorCoef_wave CorCoef_contrl
for roi=1:NROI
    temp=corrcoef(Signal_wave(:,roi),Signal_wave(:,NROI+roi));
    CorCoef_wave(1,roi)=temp(2);
    temp=corrcoef(Signal_contrl(:,roi),Signal_contrl(:,NROI+roi));
    CorCoef_contrl(1,roi)=temp(2);
end
clearvars temp roi

  
    
%% Thresholding ON/OFF states
clearvars threshold_2ch_wave ROI_present_contrl
for exp=1:length(ROI_contrl)
    temp=zeros(1, NROI_all(exp));
    temp(ROI_contrl{exp})=1;
    if exp==1
        ROI_present_contrl=temp;
    else
        ROI_present_contrl=cat(2, ROI_present_contrl, temp);
    end
end
clearvars temp exp 

%calculation of the threshold based on gaussian fit.
for j=1:NROI
    if (Keep_Only_Double==1 && ROI_present_contrl(j)==1) || Keep_Only_Double==0;
    [A, edges]=histcounts(Signal_2ch_wave(:,j),40);
    [xData, yData] = prepareCurveData(edges(1:end-1), A );
    ft = fittype( 'gauss2' );
    sig_mean=mean(Signal_2ch_wave(:,j));
    sig_std=std(Signal_2ch_wave(:,j));
    seed=[max(A) sig_mean-sig_std sig_std/2 max(A) sig_mean+sig_std sig_std/2];
    [fitresult, gof] = fit( xData, yData, ft, 'StartPoint', seed);
    gauss=coeffvalues(fitresult);
    if gauss(2)>gauss(5)
        if gauss(5)+gauss(6)*2<gauss(2)
            threshold_2ch_wave(j)=gauss(5)+gauss(6)*2;
        else
            threshold_2ch_wave(j)=mean(gauss([2 5]));
        end
    else
        if gauss(2)+gauss(3)*2<gauss(5)
            threshold_2ch_wave(j)=gauss(2)+gauss(3)*2;
        else
            threshold_2ch_wave(j)=mean(gauss([2 5]));
        end
    end
        else threshold_2ch_wave(j)=max(Signal_2ch_wave(:,j))+1;
    end
end
clearvars seed b edges xData yData gauss fitresult gof ft opts.Display opts opts.Lower opts.StartPoint A j sig_mean sig_std Keep_Only_Double
%%

clearvars threshold_2ch_contrl 

for j=1:NROI
    if ROI_present_contrl(j)==1;
        [A, edges]=histcounts(Signal_2ch_contrl(:,j),40);
        [xData, yData] = prepareCurveData(edges(1:end-1), A );
        ft = fittype( 'gauss2' );
        sig_mean=mean(Signal_2ch_contrl(:,j));
        sig_std=std(Signal_2ch_contrl(:,j));
        seed=[max(A) sig_mean-sig_std sig_std/2 max(A) sig_mean+sig_std sig_std/2];
        [fitresult, gof] = fit( xData, yData, ft, 'StartPoint', seed);
        gauss=coeffvalues(fitresult);
        if gauss(2)>gauss(5)
            if gauss(5)+gauss(6)*2<gauss(2)
                threshold_2ch_contrl(j)=gauss(5)+gauss(6)*2;
            else
                threshold_2ch_contrl(j)=mean(gauss([2 5]));
            end
        else
            if gauss(2)+gauss(3)*2<gauss(5)
                threshold_2ch_contrl(j)=gauss(2)+gauss(3)*2;
            else
                threshold_2ch_contrl(j)=mean(gauss([2 5]));
            end
        end
    else threshold_2ch_contrl(j)=max(Signal_2ch_contrl(:,j))+1;
    end
end

clearvars seed b edges xData yData gauss fitresult gof ft opts.Display opts opts.Lower opts.StartPoint A j sig_mean sig_std


%% Defining OnState periods
OnState_wave=bsxfun(@ge,Signal_2ch_wave,threshold_2ch_wave); %compares column-wise signal_sum with threshold
Signal_2ch_OnState_wave=Signal_2ch_wave.*OnState_wave;
Signal_2ch_OnState_wave(Signal_2ch_OnState_wave==0)=nan;

OnState_wave=repmat(OnState_wave,1,2); 
Signal_OnState_wave=Signal_wave.*OnState_wave;
Signal_OnState_wave(Signal_OnState_wave==0)=nan;

% calculation of correlation coefficient between blue and red traces OnState 
for j=1:NROI
    temp_blue=Signal_OnState_wave(:,j);
    temp_blue(isnan(temp_blue)) = [];
    temp_red=Signal_OnState_wave(:,NROI+j);
    temp_red(isnan(temp_red)) = [];
    temp=corrcoef(temp_blue,temp_red);
    if size(temp)>1
        CorCoef_wave(2,j)=temp(2);
    else CorCoef_wave(2,j)=nan;
    end
end

OnState_contrl=bsxfun(@ge,Signal_2ch_contrl,threshold_2ch_contrl); %compares column-wise signal_sum with threshold
Signal_2ch_OnState_contrl=Signal_2ch_contrl.*OnState_contrl;
Signal_2ch_OnState_contrl(Signal_2ch_OnState_contrl==0)=nan;

OnState_contrl=repmat(OnState_contrl,1,2); 
Signal_OnState_contrl=Signal_contrl.*OnState_contrl;
Signal_OnState_contrl(Signal_OnState_contrl==0)=nan;

% calculation of correlation coefficient between blue and red traces OnState 
for j=1:NROI
    temp_blue=Signal_OnState_contrl(:,j);
    temp_blue(isnan(temp_blue)) = [];
    temp_red=Signal_OnState_contrl(:,NROI+j);
    temp_red(isnan(temp_red)) = [];
    temp=corrcoef(temp_blue,temp_red);
    if size(temp)>1
        CorCoef_contrl(2,j)=temp(2);
    else CorCoef_contrl(2,j)=nan;
    end
end

clearvars temp temp_blue temp_red OnState OnState_OffPhase n j

%% calculate blue/sum ratio

Ratio_wave=zeros(Nframes,NROI);
for i=1:NROI
    Ratio_wave(:,i)=Signal_OnState_wave(:,i)./Signal_2ch_OnState_wave(:,i);
end
clearvars i

split_wave=nanmean(Ratio_wave);
Sig_FFT_wave=zeros(Nframes,NROI);
Sig_FFT_wave(:,split_wave>0.5)=Signal_OnState_wave(:,split_wave>0.5);
Sig_FFT_wave(:,split_wave<=0.5)=Signal_OnState_wave(:,NROI+find(split_wave<=0.5));

Ratio_contrl=zeros(Nframes(1),NROI);
for i=1:NROI
    Ratio_contrl(:,i)=Signal_OnState_contrl(:,i)./Signal_2ch_OnState_contrl(:,i);
end
clearvars i

split_contrl=nanmean(Ratio_contrl);
Sig_FFT_contrl=zeros(Nframes,NROI);
Sig_FFT_contrl(:,split_contrl>0.5)=Signal_OnState_contrl(:,split_contrl>0.5);
Sig_FFT_contrl(:,split_contrl<=0.5)=Signal_OnState_contrl(:,NROI+find(split_contrl<=0.5));

%% FFT score calculation
clearvars score_FFT_wave
int=20; % interval for STD calculation in Hz
trh=1; %length of the minimum ON interval (in periods of the wave, x4 in frames)
length_trh=100; % threshold for the length of the OnState trace
normalize=0;
score_trh=5;

[Signal_aligned_wave, ~]=align_ONState(Sig_FFT_wave, trh, normalize);
for roi=1:NROI
    score_FFT_wave(roi)= FFT_score(Signal_aligned_wave{roi}', length_trh, int, Tacq);  
end

ROI_selected_FFT_wave=find(score_FFT_wave>score_trh);


[Signal_aligned_contrl, ~]=align_ONState(Sig_FFT_contrl, trh, normalize);
for roi=1:NROI
    score_FFT_contrl(roi)= FFT_score(Signal_aligned_contrl{roi}', length_trh, int, Tacq);  
end

ROI_selected_FFT_contrl=find(score_FFT_contrl>score_trh);
clearvars int trh length_trh normalize score_trh