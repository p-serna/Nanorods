%% --- user defined variables --- %
root='D:\LocalData\aludwig\Data\MATLAB\scripts for Bar_Ilan\sample_rec_vsNR\image files';

folder_wave='wave'; 
folder_ctrl='control';

file_wave='cell1_2.tif'; %file with voltage modulation
file_ctrl='cell1_2.tif'; %file without voltage modulation

file_mask='cell1_BeRST.tif'; %file with BeRST recording


ROIsize=5;
calib_file='calib.tif'; %calibration image for optosplit; 
%in the same folder there must be a calib.csv file with parwise coordinates of
%ROIs on the left and right side of the image.

%% extract some info from the files, read calib file and calculate ROI displasment between two channels (dx and dy)

imcalib=imread(fullfile(root, calib_file));
pts_calib = ReadROIs(fullfile(root, [calib_file(1:end-3) 'csv'])); %coordinates of ROIs for the calibration image
pts_calib = CenterROIs(imcalib,pts_calib,9);
dx=round(mean(pts_calib(2:2:end,1)-pts_calib(1:2:end,1)));
dy=round(mean(pts_calib(2:2:end,2)-pts_calib(1:2:end,2)));

info = imfinfo(fullfile(root, folder_wave, file_wave));
info = info(1, :);
NFramesStr = regexp(info.ImageDescription, 'images=(\d*)', 'tokens');
Nframes = str2double(NFramesStr{1});
TacqStr = regexp(info.ImageDescription, 'finterval=(\d?\D?\d*)','tokens');
Tacq= str2double(TacqStr{1}); %frame rate


% Generating stimulation trace
Stim=repmat([0 0 -60 -60],1,Nframes/4); %stimulation trace

clearvars imcalib pts_calib TacqStr NFramesStr info f

%% Make mask

f=3; %averaging factor for the mask

%EMCCD-CMOS calib results
[ax ay bx by]=EMCCD_CMOS_calib([root, '\EMCCD-CMOS calib']); %use EMCCD and CMOS calibration images to extract coordinates displasment

% extract some info from the BeRST file
info = imfinfo(fullfile(root, file_mask));
info = info(1, :);
NFramesStr = regexp(info.ImageDescription, 'images=(\d*)', 'tokens');
Nframes_mask = str2double(NFramesStr{1});
TacqStr = regexp(info.ImageDescription, 'finterval=(\d?\D?\d*)','tokens');
Tacq_mask = str2double(TacqStr{1});


clearvars imcalib pts_calib TacqStr pts_calib_EMCCD pts_calib_CMOS info NFramesStr

% Generating stimulation trace
Stim_mask=repmat([-60 0], 1, Nframes_mask/2)'; %stimulation for the BeRST recording

% read mask movie 
Movie_mask=ReadLongMovie(fullfile(root, file_mask), Nframes_mask); %read BeRST recording


SizeX=size(Movie_mask, 1);
SizeY=size(Movie_mask, 2);

% calculation of correlation coefficient matrix between mask and a stimulation trace

clearvars CorCoef_mask 
for x=1:SizeX
    for y=1:SizeY
         temp=corrcoef(Movie_mask(x,y,:),Stim_mask);
         CorCoef_mask(x,y)=temp(2);
    end
end

clearvars temp x y;


% pseudocolor plot CorCoef_mask and creation of the Mask_EMCCD
C=CorCoef_mask;

T=1; %if T=0, pseudocolor CorCoef image is plotted, if T=1 the mask is plotted

thld=mean(CorCoef_mask(:))+2.2*std(CorCoef_mask(:)); % mask thresholding
X=[1:size(C,2)];
Y=[1:size(C,1)];

if T
    C(C<thld)=0;
    C(C>=thld)=1;
    Mask_EMCCD=C;
    figure;
    imshow(Mask_EMCCD);
else
figure('position', [400, 50, 800, 800]);
pcolor(X,Y,C)    
end

clearvars X Y C thld T 

% Filter and Expand the Mask
temp1=movmean (Mask_EMCCD,5,1,'omitnan');
temp2=movmean (temp1,5,2,'omitnan');
temp3=Mask_EMCCD;
temp3(temp2<=0.15)=nan;
temp3(temp3==0)=nan;

temp1=movmean (temp3,f,1,'omitnan');
Mask_EMCCD_exp=movmean (temp1,f,2,'omitnan');

figure;
imshow(Mask_EMCCD_exp);

clearvars temp1 temp2 temp3

% Translation of Mask_EMCCD to Mask_CMOS
Mask_CMOS=calib(Mask_EMCCD_exp, ax, ay, bx, by); %creation of CMOS mask based of EMCCD mask and coordinates displacement values
figure
imshow(Mask_CMOS)

clearvars Mask_EMCCD Mask_EMCCD_exp Nframes_mask Movie_mask ax ay bx by CorCoef_mask SumMovie_mask Tacq_mask SizeX SizeY Stim_mask

%% WAVE movie  (recording with voltage modulation)
%% reading wave movie

Movie=ReadLongMovie(fullfile(root, folder_wave,file_wave), Nframes);
SumMovie_wave=sum(Movie,3);


%% Generating ROIs within CMOS mask

clearvars pts
ROIsize=5; %size of ROI in pixels (5x5)
n=1;
SE = strel('square',5);
J = imdilate(SumMovie_wave,SE); 
Mask_CMOS(Mask_CMOS==0)=nan;
for x=2*ROIsize+1:size(Mask_CMOS,1)-2*ROIsize
    for y=2*ROIsize:size(Mask_CMOS,2)-2*ROIsize
        if ~isnan(Mask_CMOS(x,y))&& J(x, y)==SumMovie_wave(x, y)
            pts(n,:)=[y x]; %coordinates of bright clusters inside the mask; only left (blue) part of the splitted image is considered
            n=n+1;
        end
    end
end

clearvars n x y SE J
%%  check if ROI after translation to the right (red) side are still inside the image boundaries
x=pts(:,1)+dx; % translation to the right side
y=pts(:,2)+dy;

a=x+floor(ROIsize/2)<size(SumMovie_wave,2); 
x=x.*a;
y=y.*a;

a=y+floor(ROIsize/2)<size(SumMovie_wave,1);
x=x.*a;
y=y.*a;

a=y-floor(ROIsize/2)>0;
x=x.*a;
y=y.*a;

x(x==0)=[]; % removal of ROIs that are outside the boundaries
y(y==0)=[];
pts=[x-dx y-dy];


clearvars x y n s J SE a

%% reading ROIs
pts_blue=pts(:,1:2);
pts_blue_center=CenterROIs(SumMovie_wave,pts_blue,ROIsize);

% transtlate blue ROIs to red side of the camera and re-center
pts_red = pts_blue_center+repmat([dx dy],size(pts_blue_center,1),1);
pts_red_center = pts_red;%CenterROIs(SumMovie_wave,pts_red,ROIsize);

ROIs_Movie = MakeROIMoive(Movie,[pts_blue_center;pts_red_center],ROIsize); %reading values for pixels inside ROIs

NROI_auto=size(ROIs_Movie,4)/2;

clearvars pts_blue pts_red;

%% Determining ROI signal
% Creation of a matrix with all ROI data. Dimensions of Linear_ROIs_Moive: (ROI Pixels,Frame,ROI#)
Linear_ROIs_Movie=reshape(ROIs_Movie,ROIsize^2,1,Nframes,NROI_auto*2); %linearizing ROI movie from  5x5x(frame number)x(ROIs number) to 25xx(frame number)x(ROIs number)
Linear_ROIs_Movie=permute(Linear_ROIs_Movie,[1,3,4,2]);

% Sorting ROI Pixels in Linear_ROIs_Moive from min to max
for j=1:NROI_auto*2 % loop over all ROI
    for i=1:Nframes % loop over all frames
        Linear_ROIs_Movie(:,i,j)=sortrows(Linear_ROIs_Movie(:,i,j));
    end
end

Linear_ROIs_Movie_Back=Linear_ROIs_Movie(1:ROIsize,:,:); %background pixels (5 smallest brightness)
Linear_ROIs_Movie_Signal=Linear_ROIs_Movie(ROIsize+1:end,:,:); %signal pixels (20 highest brightness)


Signal_auto=squeeze(sum(Linear_ROIs_Movie_Signal)-mean(Linear_ROIs_Movie_Back)*(ROIsize^2-ROIsize))/Tacq; %calculation of the signal as sum of signal pixels ...
%minus mean of background pixels multiplied by 20; first all blue ROIs, then corresponding red ROIs 

Signal_auto_2ch=Signal_auto(:,1:NROI_auto)+Signal_auto(:,NROI_auto+1:end); %sum of blue and red channels
clearvars  Linear_ROIs_Movie_Back Linear_ROIs_Movie_Signal j i

%% classifier (selection of blinking traces)

temp1=mean(Signal_auto_2ch);
temp2=mean((Signal_auto_2ch-temp1).^2);
temp3=mean((Signal_auto_2ch-temp1).^3)./temp2.^(3.0/2.0);
temp4=mean((Signal_auto_2ch-temp1).^4)./temp2.^(4.0/2.0);
classifier=(temp3.^2+1)./temp4;

clearvars temp1 temp2 temp3 temp4
% selecting ROIs

clearvars ROI_selected
ROI_selected=find(classifier>0.45); %selection of ROI with classifier above the threshold

if sum(ROI_selected)==0
    disp('no ROI selected')
    return;
end

NROI=size(ROI_selected,2);
Signal_wave=Signal_auto(:,[ROI_selected NROI_auto+ROI_selected]); %only selected ROI's signal
Signal_wave_2ch=Signal_auto_2ch(:,ROI_selected);


pts_blue=pts_blue_center(ROI_selected,:);
pts_red=pts_red_center(ROI_selected,:);

clearvars Signal_auto Signal_auto_2ch

%% Plot ROIs 
RB=1;
norm=1;
PlotROIs(SumMovie_wave,[pts_blue;pts_red],ROIsize,RB, norm)
clearvars RB norm;

%% CONTROL movie %movie without voltage modulation
%%  reading control file

Movie=ReadLongMovie(fullfile(root, folder_ctrl,file_ctrl), Nframes);

SumMovie_ctrl=sum(Movie,3);
pts_blue_center=CenterROIs(SumMovie_ctrl,pts_blue,ROIsize);
pts_red_center=CenterROIs(SumMovie_ctrl,pts_red,ROIsize);

ROIs_Movie = MakeROIMoive(Movie,[pts_blue_center;pts_red_center],ROIsize); %reading ROI movie, the same ROI as in the recording with voltage modulation are used

%% Determining ROI signal in control movie
% Creation of a matrix with all ROI data. Dimensions of Linear_ROIs_Moive: (ROI Pixels,Frame,ROI#)
Linear_ROIs_Movie=reshape(ROIs_Movie,ROIsize^2,1,Nframes,NROI*2);
Linear_ROIs_Movie=permute(Linear_ROIs_Movie,[1,3,4,2]);

% Sorting ROI Pixels in Linear_ROIs_Moive from min to max
for j=1:NROI*2 % loop over all ROI
    for i=1:Nframes % loop over all frames
        Linear_ROIs_Movie(:,i,j)=sortrows(Linear_ROIs_Movie(:,i,j));
    end
end

Linear_ROIs_Movie_Back=Linear_ROIs_Movie(1:ROIsize,:,:);
Linear_ROIs_Movie_Signal=Linear_ROIs_Movie(ROIsize+1:end,:,:);


Signal_ctrl=squeeze(sum(Linear_ROIs_Movie_Signal)-mean(Linear_ROIs_Movie_Back)*(ROIsize^2-ROIsize))/Tacq;
Signal_ctrl_2ch=Signal_ctrl(:,1:NROI)+Signal_ctrl(:,NROI+1:end);
clearvars  Linear_ROIs_Movie_Back Linear_ROIs_Movie_Signal j i

%% classifier control

temp1=mean(Signal_ctrl_2ch);
temp2=mean((Signal_ctrl_2ch-temp1).^2);
temp3=mean((Signal_ctrl_2ch-temp1).^3)./temp2.^(3.0/2.0);
temp4=mean((Signal_ctrl_2ch-temp1).^4)./temp2.^(4.0/2.0);
classifier=(temp3.^2+1)./temp4;

clearvars temp1 temp2 temp3 temp4
% selecting ROIs

clearvars ROI_selected_contrl
ROI_ctrl=find(classifier>0.45); %finding ROIs that are blinking in control

if sum(ROI_ctrl)==0
    disp('no ROI selected')
    return;
end

clearvars j button x y m n limy limx temp fitobject gof STD
        

Signal_wave=Signal_wave(:,[ROI_ctrl NROI+ROI_ctrl]); %only ROI present in both wave and control
Signal_wave_2ch=Signal_wave_2ch(:,ROI_ctrl);
Signal_ctrl=Signal_ctrl(:,[ROI_ctrl NROI+ROI_ctrl]); %only ROI present in both wave and control
Signal_ctrl_2ch=Signal_ctrl_2ch(:,ROI_ctrl);
NROI=size(ROI_ctrl,2);
%% Plot ROIs present in control
RB=1;
norm=1;
PlotROIs(SumMovie_wave,[pts_blue(ROI_ctrl,:);pts_red(ROI_ctrl,:)],ROIsize,RB, norm)
clearvars RB norm;

%% Thresholding ON/OFF states (blinking threshold)
%WAVE
clearvars threshold_wave_2ch

for j=1:NROI
    [A, edges]=histcounts(Signal_wave_2ch(:,j),40);
    [xData, yData] = prepareCurveData(edges(1:end-1), A );
    ft = fittype( 'gauss2' );
    sig_mean=mean(Signal_wave_2ch(:,j));
    sig_std=std(Signal_wave_2ch(:,j));
    seed=[max(A) sig_mean-sig_std sig_std/2 max(A) sig_mean+sig_std sig_std/2];
    [fitresult, gof] = fit( xData, yData, ft, 'StartPoint', seed);
    gauss=coeffvalues(fitresult);
    if gauss(2)>gauss(5)
        if gauss(5)+gauss(6)*2<gauss(2)
            threshold_wave_2ch(j)=gauss(5)+gauss(6)*2;
        else
            threshold_wave_2ch(j)=mean(gauss([2 5]));
        end
    else
        if gauss(2)+gauss(3)*2<gauss(5)
            threshold_wave_2ch(j)=gauss(2)+gauss(3)*2;
        else
            threshold_wave_2ch(j)=mean(gauss([2 5]));
        end
    end
end
clearvars seed b edges xData yData gauss fitresult gof ft opts.Display opts opts.Lower opts.StartPoint A j sig_mean sig_std Keep_Only_Double

%CONTROL
clearvars threshold_ctrl_2ch

for j=1:NROI
    
    [A, edges]=histcounts(Signal_ctrl_2ch(:,j),40);
    [xData, yData] = prepareCurveData(edges(1:end-1), A );
    ft = fittype( 'gauss2' );
    sig_mean=mean(Signal_ctrl_2ch(:,j));
    sig_std=std(Signal_ctrl_2ch(:,j));
    seed=[max(A) sig_mean-sig_std sig_std/2 max(A) sig_mean+sig_std sig_std/2];
    [fitresult, gof] = fit( xData, yData, ft, 'StartPoint', seed);
    gauss=coeffvalues(fitresult);
    if gauss(2)>gauss(5)
        if gauss(5)+gauss(6)*2<gauss(2)
            threshold_ctrl_2ch(j)=gauss(5)+gauss(6)*2;
        else
            threshold_ctrl_2ch(j)=mean(gauss([2 5]));
        end
    else
        if gauss(2)+gauss(3)*2<gauss(5)
            threshold_ctrl_2ch(j)=gauss(2)+gauss(3)*2;
        else
            threshold_ctrl_2ch(j)=mean(gauss([2 5]));
        end
    end
    
end

%clearvars seed b edges xData yData gauss fitresult gof ft opts.Display opts opts.Lower opts.StartPoint A j sig_mean sig_std

%% Plotting signal_2ch and threshold. Sum of channels. Random ROI
sig=Signal_wave_2ch;
thresh=threshold_wave_2ch;
roi=randi(NROI);


figure('position', [200, 400, 1200, 600]);
plot(sig(:,roi),'.-') 
hold on

plot([0 Nframes],[thresh(roi) thresh(roi)],'--','LineWidth',3);
title(['time trace ROI', num2str(roi)]);
set(gca,'FontSize',30)
xlabel('frame #, 10ms/frame')
ylabel('cps')
axis('tight');

figure('position', [1500, 400, 800, 600]);
histogram(sig(:,roi));
temp=ylim;
hold on
plot([thresh(roi) thresh(roi)], temp,'--','LineWidth',3)


clearvars sig thresh temp roi
%% Defining OnState periods (non-blinking periods)
OnState_wave=bsxfun(@ge,Signal_wave_2ch,threshold_wave_2ch); %compares column-wise signal_sum with threshold
Signal_wave_2ch_OnState=Signal_wave_2ch.*OnState_wave;
Signal_wave_2ch_OnState(Signal_wave_2ch_OnState==0)=nan;

OnState_wave=repmat(OnState_wave,1,2); 
Signal_wave_OnState=Signal_wave.*OnState_wave;
Signal_wave_OnState(Signal_wave_OnState==0)=nan;


OnState_ctrl=bsxfun(@ge,Signal_ctrl_2ch,threshold_ctrl_2ch); %compares column-wise signal_sum with threshold
Signal_ctrl_2ch_OnState=Signal_ctrl_2ch.*OnState_ctrl;
Signal_ctrl_2ch_OnState(Signal_ctrl_2ch_OnState==0)=nan;

OnState_ctrl=repmat(OnState_ctrl,1,2); 
Signal_ctrl_OnState=Signal_ctrl.*OnState_ctrl;
Signal_ctrl_OnState(Signal_ctrl_OnState==0)=nan;


clearvars OnState_wave OnState_ctrl

%% calculate blue/sum ratio (ratio between the two channels)

Ratio_wave=zeros(Nframes,NROI);
for i=1:NROI
    Ratio_wave(:,i)=Signal_wave_OnState(:,i)./Signal_wave_2ch_OnState(:,i);
end
clearvars i

split_wave=nanmean(Ratio_wave); %selecting the best channel (highest signal of the two) for FFT calculation
Sig_wave_FFT=zeros(Nframes,NROI);
Sig_wave_FFT(:,split_wave>0.5)=Signal_wave_OnState(:,split_wave>0.5);
Sig_wave_FFT(:,split_wave<=0.5)=Signal_wave_OnState(:,NROI+find(split_wave<=0.5));

Ratio_ctrl=zeros(Nframes(1),NROI);
for i=1:NROI
    Ratio_ctrl(:,i)=Signal_ctrl_OnState(:,i)./Signal_ctrl_2ch_OnState(:,i);
end
clearvars i

split_ctrl=nanmean(Ratio_ctrl);
Sig_ctrl_FFT=zeros(Nframes,NROI);
Sig_ctrl_FFT(:,split_ctrl>0.5)=Signal_ctrl_OnState(:,split_ctrl>0.5);
Sig_ctrl_FFT(:,split_ctrl<=0.5)=Signal_ctrl_OnState(:,NROI+find(split_ctrl<=0.5));

clearvars split_wave split_ctrl
%% FFT score calculation
clearvars score_wave_FFT score_ctrl_FFT
int=20; % interval for STD calculation in Hz
trh=1; %length of the minimum ON interval (in periods of the wave, x4 in frames)
length_trh=100; % threshold for the length of the OnState trace
normalize=0;
score_trh=5;

[Signal_wave_aligned, ~, ~]=align_ONState(Sig_wave_FFT, trh, normalize);
for roi=1:NROI
    score_wave_FFT(roi)= FFT_score(Signal_wave_aligned{roi}', length_trh, int, Tacq);  
end

ROI_selected_wave_FFT=find(score_wave_FFT>score_trh);


[Signal_ctrl_aligned, ~, ~]=align_ONState(Sig_ctrl_FFT, trh, normalize);
for roi=1:NROI
    score_ctrl_FFT(roi)= FFT_score(Signal_ctrl_aligned{roi}', length_trh, int, Tacq);  
end

ROI_selected_ctrl_FFT=find(score_ctrl_FFT>score_trh);

clearvars int trh length_trh normalize score_trh

