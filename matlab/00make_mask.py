#%% --- user defined variables --- %
#folder='C:\Users\ludwig\Documents\MATLAB\19_02_05_pd1_02_div4_NR_BeRST';
#file_wave='cell8_BeRST.tif'; 
#f=3; %averaging factor for the mask

#%EMCCD-CMOS calib results
#[ax ay bx by]=EMCCD_CMOS_calib([folder, '\EMCCD_CMOS_calib']);
import sys
from tifmethods import readtifInfo, readBigTifFile, readtifImage

def createMask(file_wave,calibpar = None, folder = None, avgf = 3):
    '''file_wave: name of the file with pixels very correlated with voltage wave
       calibpar: parameters of calibration ax,ay,bx,by
       avgf: averaging factor for the mask 
    '''
    if folder is None:
        if file_wave.find('/') >=0:
            fsp = file_wave.split('/')
            folder = '/'.join(fsp[:-1])+'/'
            file_wave = fsp[-1]
        if file_wave.find('\\') >=0:
            fsp = file_wave.split('/')
            folder = '\\'.join(fsp[:-1])+'\\'
            file_wave = fsp[-1]
    if calibpar is None:
        calibpar = EMCCD_CMOS_calib(folder+'EMCCD_CMOS_calib')

    # We extract basic info: number of frames and interval (is this the average?)
    info = readtifInfo(folder+'image files/'+file_wave)

    NframesStr = [d for d in info[270].split('\n') if d.find('frames')>=0][0].split('=')[-1]
    # Why float and not integer
    Nframes_wave = float(NFramesStr)

    TacqStr = [d for d in info[270].split('\n') if d.find('finterval')>=0][0].split('=')[-1]
    Tacq_wave = float(TacqStr)

    # Generate stimulation trace
    Stim_wave = (((arange(Nframes_wave).reshape(Nframes_wave//2,2)+1)%2)*(-60.0)).flatten()

    # read movie wave
    Movie_wave=readBigTifFile(folder+'image files/'+file_wave), Nframes_wave);

SumMovie_wave=sum(Movie_wave,3);

SizeX=size(SumMovie_wave, 1);
SizeY=size(SumMovie_wave, 2);


% %Plot images
% temp=SumMovie_wave;
% img=temp/max(max(temp));
% img1=imadjust(img);
% figure;
% imshow(img1)

clearvars img img1 temp


% calculation of correlation coefficient matrix between NR wave movie and a stimulation trace
SizeX=size(Movie_wave,1);
SizeY=size(Movie_wave,2);
clearvars CorCoef_wave 
for x=1:SizeX
    for y=1:SizeY
         temp=corrcoef(Movie_wave(x,y,:),Stim_wave);
         CorCoef_wave(x,y)=temp(2);
    end
end

clearvars temp x y;
% 


%% pseudocolor plot CorCoef_wave and creation of the Mask_EMCCD
C=CorCoef_wave;
T=1;
thld=mean(CorCoef_wave(:))+2.2*std(CorCoef_wave(:));%0.12;
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

%% 

% Filter and Expand the Mask
temp1=movmean (Mask_EMCCD,5,1,'omitnan');
temp2=movmean (temp1,5,2,'omitnan');
temp3=Mask_EMCCD;
temp3(temp2<=0.15)=nan;
temp3(temp3==0)=nan;

temp1=movmean (temp3,f,1,'omitnan');
Mask_EMCCD_exp=movmean (temp1,f,2,'omitnan');

% figure;
% imshow(Mask_EMCCD);
% figure;
% imshow(temp3);
figure;
imshow(Mask_EMCCD_exp);

clearvars temp1 temp2 temp3
% 
 
%% 

% Translation of Mask_EMCCD to Mask_CMOS
Mask_CMOS=calib(Mask_EMCCD_exp, ax, ay, bx, by);
%Mask_CMOS=Mask_CMOS(256:767,1:512);
figure
imshow(Mask_CMOS)



% Saving 

temp=Mask_EMCCD;
temp1=temp/max(max(temp));
imwrite(temp1,fullfile(folder, [file_wave(1:end-10), '_Mask_EMCCD.tif']),'tif');

temp=Mask_CMOS;
temp1=temp/max(max(temp));
imwrite(temp1,fullfile(folder, [file_wave(1:end-10), '_Mask_CMOS.tif']),'tif');

clearvars temp temp1
% 

clearvars  ans Movie_wave SizeX SizeY Nframes_wave Tacq_wave NFramesStr Stim_wave info SumMovie_wave Mask_EMCCD_exp Mask_EMCCD CorCoef_wave 
save(fullfile(folder, 'mat files', [file_wave(1:end-10), '_mask']))
